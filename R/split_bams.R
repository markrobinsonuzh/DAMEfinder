
#' Detect allele-specific methylation from a bam file
#'
#' The function takes a bam (from bismark) and vcf file for each sample.
#' For each SNP contained in the vcfile it calculates the proportion of
#' methylated reads for each CpG site at each allele. At the end it returns
#' (saves to working directory) a GRanges list, where each GRanges contains
#' all the CpG sites overlapping the reads containing a specific SNP.
#'
#' @param bam_files List of bam files.
#' @param vcf_files List of vcf files.
#' @param sample_names Names of files in the list.
#' @param reference_file fasta file used to generate the bam files.
#' @param coverage Minimum number of reads covering a CpG site on each allele. default = 2
#' @param cores Number of cores to use. See package {parallel} for description of core. default = 1
#'
#' @return A list of GRanges for each sample. Each list is saved in separate .rds files
#' @examples
#' @importFrom BiocGenerics start
#' @importFrom BiocGenerics end
#' @importFrom S4Vectors mcols
#'
#'
#' @export
split_bams <- function(bam_files, vcf_files, sample_names, reference_file, coverage = 2, cores = 1){

message("Reading reference file")
fa <- open(Rsamtools::FaFile(reference_file))

for(samp in 1:length(sample_names)){

  message(sprintf("Running sample %s", sample_names[samp]))
  chrom <- "1"

  message("Reading VCF file")
  vcf <- vcfR::getFIX(vcfR::read.vcfR(vcf_files[samp], verbose = FALSE))
  bam.file <- bam_files[samp]

  message(sprintf("Processing chromosome %s",chrom))

  snp.table <- parallel::mcmapply(function(t, u, v){

    #Appearently applying to a GRanges takes too long ? so I will create one each time
    snp <- GenomicRanges::GRanges(seqnames = t, IRanges::IRanges(start = as.integer(u), width = 1))

    if(length(grep(t, c(as.character(1:21), "X", "Y"))) == 0){
      message("Bad chrom")
      return(NULL)
      next
      }
    if(t != chrom){message(sprintf("Processing chromosome %s",t))}
    chrom <- t

    #Get the reads that align to the specified SNP
    alns.pairs <- GenomicAlignments::readGAlignmentPairs(bam.file,
                                      param = Rsamtools::ScanBamParam(
                                        tag= c("MD","XM","XR","XG"),
                                        which=snp),
                                      use.names = TRUE)

    alns <- unlist(alns.pairs) #unpaired
    split <- splitReads(alns, v, snp)

    alt.reads <- split$alt.reads
    ref.reads <- split$ref.reads

    if(length(ref.reads) <= 1 | length(alt.reads) <= 1){
      message("0 or 1 read in one or both alleles")
      return(NULL)
      next
    }

    ####get reference and CpG positions####---------------------------------------

    #Get limits for calling methylation and grabing reference sequence
    left <- min(start(alns))
    right <- max(end(alns))
    window <- GenomicRanges::GRanges(GenomeInfoDb::seqnames(snp), IRanges::IRanges(left, right))

    #open reference seq to get correct CpG locations within that reference chunk
    dna <- Rsamtools::scanFa(fa, param=window)

    cgsite <- stringr::str_locate_all(dna, "CG")[[1]][,1] #also look at GpCs?

    if(length(cgsite) < 1){
      message("No CpG sites associated to this SNP")
      return(NULL)
      next
    }
    mepos <- cgsite/nchar(dna) #location of CpG sites

    ##### Use MD tag from bam to extract methylation status for each site detected above ####----------------------------
    conversion <- vapply(alns.pairs, function(x){

      #change C locations for G locations if reference context is different
      if(mcols(x)$XG[1] == "GA"){
        cgsite <- cgsite + 1
      }

      #fuse pair info in one row
      conversion <- matrix(NA,1, length(cgsite))

      read.start <- start(x) - left + 1 #start of read
      read.end <- end(x) - left + 1 #end

      for(pair in 1:2){

        tag <- getMD(x[pair])
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

                #if there is a mismatch, the MDtag shows what the base is in the ref. So,
                #if it shows a C or G, it means that the ref is C or G, and the read has T or A,
                #therefore it is unmethylated (converted)
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
    GR <- GenomicAlignments::GRanges(gsub("chr","",chrom), IRanges::IRanges(start = left + cgsite, width = 1))
    mcols(GR)$cov.ref <- ref.cov
    mcols(GR)$cov.alt <- alt.cov
    mcols(GR)$meth.ref <- ref.meth
    mcols(GR)$meth.alt <- alt.meth
    mcols(GR)$snp <- paste0(GenomeInfoDb::seqnames(snp),".",start(snp))

    #Keep CpG sites where each allele has at least 'coverage' reads.
    filt <- BiocGenerics::rowSums(cbind(mcols(GR)$cov.ref,
                                        mcols(GR)$cov.alt) >= coverage) >= 2
    if(sum(filt) < 2){
      return(NULL)
      break
    }

    gr <- GR[filt]
    return(gr)

  },t = vcf[,1], u = vcf[,2], v = vcf[,4], SIMPLIFY = F, USE.NAMES = F, mc.cores = cores)

  message(sprintf("Saving results for sample %s", sample_names[samp]))
  saveRDS(snp.table, file = sprintf("snp.table.full.%s.rds", sample_names[samp]))

}
}
