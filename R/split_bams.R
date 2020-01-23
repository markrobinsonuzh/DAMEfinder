
#' Detect allele-specific methylation from a bam file
#'
#' The function takes a bam (from bismark) and vcf file for each sample. For
#' each SNP contained in the vcfile it calculates the proportion of methylated
#' reads for each CpG site at each allele. At the end it returns (saves to
#' working directory) a GRanges list, where each GRanges contains all the CpG
#' sites overlapping the reads containing a specific SNP.
#'
#' @param bamFiles List of bam files.
#' @param vcfFiles List of vcf files.
#' @param sampleNames Names of files in the list.
#' @param referenceFile fasta file used to generate the bam files. Or 
#'   \code{DNAStringSet} with DNA sequence.
#' @param coverage Minimum number of reads covering a CpG site on each allele.
#'   Default = 2.
#' @param cores Number of cores to use. See package {parallel} for description
#'   of core. Default = 1.
#' @param verbose default = TRUE
#'
#' @return A list of GRanges for each sample. Each list is saved in a separate
#'   .rds file.
#' @examples
#' DATA_PATH_DIR <- system.file("extdata", ".", package = "DAMEfinder")
#' get_data_path <- function(file_name) file.path(DATA_PATH_DIR, file_name)
#' bamFiles <- get_data_path("NORM1_chr19_trim.bam")
#' vcfFiles <- get_data_path("NORM1.chr19.trim.vcf")
#' sampleNames <- "NORM1"
#' #referenceFile <- get_data_path("19.fa")
#'
#' #GRanges_list <- extract_bams(bamFiles, vcfFiles, sampleNames, referenceFile)
#'
#' @importFrom BiocGenerics start
#' @importFrom BiocGenerics end
#' @importFrom S4Vectors mcols
#' @importFrom S4Vectors mcols<-
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom GenomicAlignments readGAlignmentPairs
#'
#' @export
extract_bams <- function(bamFiles, vcfFiles, sampleNames, referenceFile,
                       coverage = 4, cores = 1, verbose = TRUE){

if(verbose) message("Reading reference file", appendLF = TRUE)
if(typeof(referenceFile) == "character"){
  fa <- open(Rsamtools::FaFile(referenceFile,
                               index = paste0(referenceFile, ".fai")))
  
}

#get names and indeces per sample to apply to
sample_list <- 1:length(sampleNames)
names(sample_list) <- sampleNames

lapply(sample_list, function(samp){

  message(sprintf("Running sample %s", sampleNames[samp]))

  if(verbose) message("Reading VCF file", appendLF = TRUE)
  vcf <- vcfR::getFIX(vcfR::read.vcfR(vcfFiles[samp], verbose = verbose))
  bam.file <- bamFiles[samp]

  #message(sprintf("Processing chromosome %s",chrom))
  if(verbose) message("Extracting methylation per SNP", appendLF = TRUE)
  snp.table <- parallel::mcmapply(function(t, u, v){

    #Applying to a GRanges takes too long so I create one each time
    snp <- GRanges(seqnames = t, IRanges(start = as.integer(u), width = 1))

    #Ignore non-standard chromosomes
    if(length(grep(t, c(as.character(1:22), "X", "Y"))) == 0){
      if(verbose) message("Bad chrom", appendLF = TRUE)
      return(NULL)
      }

    #message(sprintf("Processing chromosome %s",t))
    #message(".")
    chrom <- t

    #Get the reads that align to the specified SNP
    # flags <- Rsamtools::scanBamFlag(
    #   isPaired = T,
    #   isProperPair = T,
    #   isUnmappedQuery = F,
    #   hasUnmappedMate = F,
    #   isNotPassingQualityControls = F,
    #   isDuplicate = F,
    #   isSecondaryAlignment = F,
    #   isSupplementaryAlignment = F)

    params <- Rsamtools::ScanBamParam(
      #flag = flags,
      tag = c("MD","XM","XR","XG"),
      #mapqFilter = 20,
      which = snp)

    alns.pairs <- readGAlignmentPairs(bam.file,
                                      use.names = TRUE,
                                      param = params)

    alns <- unlist(alns.pairs) #unpaired
    split <- splitReads(alns, v, snp)

    alt.reads <- split$alt.reads
    ref.reads <- split$ref.reads

    if(length(ref.reads) <= 1 | length(alt.reads) <= 1){
      if(verbose) message("0 or 1 read in one or both alleles", appendLF = TRUE)
      return(NULL)
      #next
    }

    ####get reference and CpG positions####-----------------------------------

    #Get limits for calling methylation and grabing reference sequence
    left <- min(start(alns))
    right <- max(end(alns))
    window <- GRanges(GenomeInfoDb::seqnames(snp), IRanges(left, right))

    if(typeof(referenceFile) == "S4"){
      dna <- referenceFile[window]
      
    } else if(typeof(referenceFile) == "character"){
    #open reference seq to get correct CpG locations within that reference chunk
    dna <- Rsamtools::scanFa(fa, param=window)}

    cgsite <- stringr::str_locate_all(dna, "CG")[[1]][,1] #also look at GpCs?

    if(length(cgsite) < 1){
      if(verbose) message("No CpG sites associated to this SNP", appendLF = TRUE)
      return(NULL)
      #next
    }
    mepos <- cgsite / Biostrings::nchar(dna) #location of CpG sites

    ##### Use MD tag from bam to extract methylation status for each site
    ##### detected above ####----------------------------
    conversion <- vapply(alns.pairs, function(x){
      #change C locations for G locations if reference context is different
      XGcondition <- mcols(x)$XG[1]
      if(XGcondition == "GA"){
        cgsite <- cgsite + 1
      }

      #fuse pair info in one row
      conversion <- matrix(NA,1, length(cgsite))

      read.start <- start(x) - left + 1 #start of read
      read.end <- end(x) - left + 1 #end

      for(pair in 1:2){
        something <- mcols(x[pair])$MD
        tag <- getMD(something)
        MDtag <- tag$MDtag
        nucl.num <- tag$nucl.num

        #Get the locations from the specific read:
        #Count along the numerical MDtag until you reach the actual CpG start
        #Fill in conversion table with meth state for each read
        #1 is unmethylated
        #2 is methylated
        #NA is not in read

        count <- 0
        this.mepos <- cgsite - read.start[pair] + 1

        for(p in seq_along(this.mepos)){
          if(is.na(conversion[1,p])){
            if(cgsite[p] > read.end[pair] || cgsite[p] < read.start[pair] ){
              conversion[1,p] <- NA
              next} else{
                for(i in seq_along(nucl.num)){
                  count <- count + nucl.num[i]
                  if(count >= this.mepos[p]){
                    count <- 0
                    break}
                }

                #if there is a mismatch, the MDtag shows what the base is in the
                #ref. So, if it shows a C or G, it means that the ref is C or G,
                #and the read has T or A, therefore it is unmethylated
                #(converted)
                if(MDtag[i] %in% c("C", "G")){
                  conversion[1,p] <- 1 } else{conversion[1,p] <- 2}
              }
          } else{next}
        }
      }
      return(conversion)
    }, double(length(cgsite)))

    #Return methylated reads and total coverage for each allele

    if(is.null(dim(conversion))){
      ref.conv <- conversion[names(conversion) %in% ref.reads ]
      alt.conv <- conversion[names(conversion) %in% alt.reads ]

      ref.cov <- sum(!is.na(ref.conv))
      alt.cov <- sum(!is.na(alt.conv))

      ref.meth <- sum(!is.na(ref.conv) & ref.conv == 2)
      alt.meth <- sum(!is.na(alt.conv) & alt.conv == 2)

    } else{
      ref.conv <- conversion[,colnames(conversion) %in% ref.reads ]
      alt.conv <- conversion[,colnames(conversion) %in% alt.reads ]

      ref.cov <- BiocGenerics::rowSums(!is.na(ref.conv))
      alt.cov <- BiocGenerics::rowSums(!is.na(alt.conv))

      ref.meth <- BiocGenerics::rowSums(!is.na(ref.conv) & ref.conv == 2)
      alt.meth <- BiocGenerics::rowSums(!is.na(alt.conv) & alt.conv == 2)
    }

    #Build GRanges
    GR <- GRanges(gsub("chr","",chrom), IRanges(start = left + cgsite,
                                                width = 1))
    mcols(GR)$cov.ref <- ref.cov
    mcols(GR)$cov.alt <- alt.cov
    mcols(GR)$meth.ref <- ref.meth
    mcols(GR)$meth.alt <- alt.meth
    mcols(GR)$snp <- paste0(GenomeInfoDb::seqnames(snp),".",start(snp))

    #Keep CpG sites where each allele has at least 'coverage' reads.
    filt <- BiocGenerics::rowSums(cbind(GR$cov.ref,
                                        GR$cov.alt) >= coverage) >= 2
    if(sum(filt) < 2){
      if(verbose) message("No CpG sites sufficiently covered at this SNP",
                          appendLF = TRUE)
      return(NULL)
    }

    gr <- GR[filt]
    return(gr)

  },t = vcf[,1], u = vcf[,2], v = vcf[,4], SIMPLIFY = F, USE.NAMES = F,
  mc.cores = cores)

  message(sprintf("Done with sample %s", sampleNames[samp]))
  return(snp.table)
}
)
}
