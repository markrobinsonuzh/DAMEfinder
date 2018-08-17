#' Draw methylation circle plot
#'
#' Draws CpG site methylation status as points, in reads containing a specific SNP.
#' Generates one plot per bam file.
#'
#' @param snp GRanges object containing SNP location
#' @param cpgsite (optional) GRanges object containing a single CpG site location of interest
#' @param vcf.file vcf file
#' @param bam.file bismark bam file path
#' @param ref.file fasta reference file path
#' @param dame (optional) GRanges object containing a region to plot
#' @param letter.size Size of alleles drawn in plot
#' @param sample.name FIX: this is to save the vcf file to not generate it every time you run the function
#'
#' @return Plot
#' @examples
#' snp <- GRanges(19, IRanges(267039, width = 1))
#' cpgsite <- GRanges(19, IRanges(266998, width = 1))
#' bam.files <- list.files("Shared_taupo/steph/CRC.bismark.bams/", "_pe.dedupl_s.bam$")
#' bam.files <- paste0("Shared_taupo/steph/CRC.bismark.bams/", bam.files)
#' ref.file <- "data/annotation/Human/GRCH37/Bisulfite_Genome.release91/GRCh37.91.fa"
#' vcf.file <- "Shared_taupo/steph/CRC.vcfs/NORM1.het.snp.raw.vcf"
#'
#' methyl_circle_plot(snp = snp, vcf.file = vcf.file, bam.file = bam.files[5],
#' ref.file = ref.file, letter.size = 2.5, cpgsite = cpgsite)
#'
#' @importFrom S4Vectors mcols
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics start
#' @importFrom BiocGenerics end
#' @importFrom GenomeInfoDb seqnames
#' @import ggplot2
#'
#' @export

methyl_circle_plot <- function(snp, vcf.file, bam.file, ref.file, dame = NULL, letter.size = 2.5, sample.name = "sample1", cpgsite = NULL){

  message("Reading vcf file")
  if(!file.exists(paste0(sample.name,".RData"))){
    vcf <- vcfR::getFIX(vcfR::read.vcfR(vcf.file, verbose = T))
    save(vcf, file = paste0(sample.name, ".RData"))
  } else {load(paste0(sample.name, ".RData"))}

  snp.loc <- which(as.integer(vcf[,2]) == start(snp))
  ref <- vcf[,4][snp.loc]
  alt <- vcf[,5][snp.loc]

  message("Getting reads")
  alns.pairs <- GenomicAlignments::readGAlignmentPairs(bam.file,
                                    param = Rsamtools::ScanBamParam(what = c("flag","mapq","mrnm","mpos","isize","seq","qual"),
                                                         tag= c("MD","XM","XR","XG"),
                                                         which=snp),
                                    use.names = TRUE)

  alns <- unlist(alns.pairs) #unpaired

  split <- splitReads(alns, ref, snp)

  alt.reads <- split$alt.reads
  ref.reads <- split$ref.reads

  message("Sorting reads")
  alt.pairs <- alns.pairs[names(alns.pairs) %in% alt.reads] #don't unlist, paired
  ref.pairs <- alns.pairs[names(alns.pairs) %in% ref.reads]
  alns.pairs <- c(alt.pairs, ref.pairs) #I just do this to sort them

  ####get reference and CpG positions####--------------------------------------------------------

  #Get limits for plotting

  if(is.null(dame)){
    left <- min(start(alns))
    right <- max(end(alns))
    window <- GRanges(seqnames(snp), IRanges(left, right))
    } else {
      window <- dame
      left <- start(dame)
      right <- end(dame)
    }

  message("Reading reference")
  #open reference seq to get correct CpG locations within that reference chunk
  fa <- open(Rsamtools::FaFile(ref.file))
  #idx <- scanFaIndex(fa)
  dna <- Rsamtools::scanFa(fa, param=window)

  cgsite <- stingr::str_locate_all(dna, "CG")[[1]][,1] #also look at GpCs?

  if(length(cgsite) < 1){
    stop("No CpG sites associated to this SNP")
  }

  mepos <- cgsite/nchar(dna) #location of CpG sites

  ##### Use MD tag from bam to extract methylation status for each site detected above ####----------------------------
  message("Getting meth state per read-pair")
  conversion <- vapply(alns.pairs, function(x){

    #change C locations for G locations if reference context is different
    if(mcols(x)$XG[1] == "GA"){
      cgsite <- cgsite + 1
    }

    #fuse pair info in one row
    conversion <- matrix(NA,1, length(cgsite))

    #a <- mcols(x)$MD
    read.start <- start(x) - left + 1 #start of read
    read.end <- end(x) - left + 1 #end

    for(pair in 1:2){

      tag <- getMD(x[pair])
      MDtag <- tag$MDtag
      nucl.num <- tag$nucl.num

      #Get the locations from the specific read:
      #Count along the numerical MDtag until you reach the actual CpG start
      count <- 0
      this.mepos <- cgsite - read.start[pair] + 1

      #Fill in conversion table with meth state for each read
      #21 is unmethylated
      #19 is methylated
      #0 is not in read

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
          # if it shows a C or G, it means that the ref is C and the read has T or A, therefore it is unmethylated (converted)
          if(MDtag[i] %in% c("C", "G")){
            conversion[1,p] <- 21 } else{conversion[1,p] <- 19}
            }
          } else{next}
        }
      }
    return(conversion)
  }, double(length(cgsite)))

  #Remove reads without CpG sites
  rem <- BiocGenerics::colSums(is.na(conversion)) >= dim(conversion)[1]
  conversion <- conversion[,!rem]

  alt.reads <- alt.reads[alt.reads %in% colnames(conversion)]
  ref.reads <- ref.reads[ref.reads %in% colnames(conversion)]
  alns.pairs <- alns.pairs[names(alns.pairs) %in% colnames(conversion)]

  #### Get snp position in window ####------------------------------------------------------------

  letter <- c(rep(alt, length(alt.reads)), rep(ref, length(ref.reads)))

  snp.start <- start(snp) - left + 1

  #and cpg of interest
  if(!is.null(cpgsite)){
    cpg.start <- start(cpgsite) - left + 1
  }else{cpg.start = NA}

  #### plot ####---------------------------------------------------------------------------------

  message("Plotting")

  #data for points
  d <- data.frame(CpG=rep(cgsite,length(alns.pairs)),
                  read=rep(1:length(alns.pairs), each=length(cgsite)),
                  value=as.vector(conversion))

  #data for segments
  reads <- d$read

  xstart <- sapply(reads, function(x){
    newd <- d[which(d$read == x),]
    vals <- which(!is.na(newd$value))
    min(newd$CpG[vals], snp.start)
  })

  xend <- sapply(reads,function(x){
    newd <- d[which(d$read == x),]
    vals <- which(!is.na(newd$value))
    max(newd$CpG[vals], snp.start)
  })

  d2 <- unique(data.frame(xstart, xend, reads))
  d2$snp <- c(rep("a", length(alt.reads)), rep("r", length(ref.reads)))

  #To manually scale the colors
  cols <- c("#0E4071", "#d55e00", "#0E4071", "#d55e00")
  names(cols) <- c("a", "r", alt, ref)

  ggplot() +
    scale_shape_identity() +
    theme_void() +
    geom_segment(data=d2, aes(x=xstart, y=reads, xend=xend, yend=reads, colour = snp), size = 0.2) +
    geom_point(data=d, aes(x=CpG, y=read, shape=value), fill = "white", size=1) +
    geom_point(aes(x=snp.start, y=1:length(alns.pairs), shape = letter, colour = letter), size=letter.size) +
    geom_point(aes(x=cpg.start, y = 0), shape = 24, size = 3, fill = "green")  +
    scale_color_manual(values = cols) +
    guides(color=FALSE)
}

#' @describeIn methylCirclePlot Draws CpG site methylation status as points, in reads containing a specific CpG site.
#'
#' @examples
#' methyl_circle_plotCpG(cpgsite = cpgsite, bam.file = bam.files[5], ref.file = ref.file,
#' dame = NULL)
#'
methyl_circle_plotCpG <- function(cpgsite = cpgsite, bam.file = bam.file, ref.file = ref.file, dame = NULL){

  alns.pairs <- GenomicAlignments::readGAlignmentPairs(bam.file,
                                    param = Rsamtools::ScanBamParam(tag= c("MD","XM","XR","XG"),
                                                         which=cpgsite),
                                    use.names = TRUE)
  alns <- unlist(alns.pairs)
  ####get reference and CpG positions####--------------------------------------------------------

  #Get limits for plotting

  if(is.null(dame)){
    left <- min(start(alns))
    right <- max(end(alns))
    window <- GRanges(seqnames(cpgsite), IRanges(left, right))
  } else {
    window <- dame
    left <- start(dame)
    right <- end(dame)
  }

  message("Reading reference")
  #open reference seq to get correct CpG locations within that reference chunk
  fa <- open(Rsamtools::FaFile(ref.file))
  #idx <- scanFaIndex(fa)
  dna <- Rsamtools::scanFa(fa, param=window)

  cgsite <- stringr::str_locate_all(dna, "CG")[[1]][,1] #also look at GpCs?

  if(length(cgsite) < 1){
    stop("No CpG sites in these reads")
  }

  mepos <- cgsite/nchar(dna) #location of CpG sites

  ##### Use MD tag from bam to extract methylation status for each site detected above ####----------------------------
  message("Getting meth state per read-pair")
  conversion <- vapply(alns.pairs, function(x){

    #change C locations for G locations if reference context is different
    if(mcols(x)$XG[1] == "GA"){
      cgsite <- cgsite + 1
    }

    #fuse pair info in one row
    conversion <- matrix(NA,1, length(cgsite))

    #a <- mcols(x)$MD
    read.start <- start(x) - left + 1 #start of read
    read.end <- end(x) - left + 1 #end

    for(pair in 1:2){

      tag <- getMD(x[pair])
      MDtag <- tag$MDtag
      nucl.num <- tag$nucl.num

      #Get the locations from the specific read:
      #Count along the numerical MDtag until you reach the actual CpG start
      count <- 0
      this.mepos <- cgsite - read.start[pair] + 1

      #Fill in conversion table with meth state for each read
      #21 is unmethylated
      #19 is methylated
      #0 is not in read

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
              # if it shows a C or G, it means that the ref is C and the read has T or A, therefore it is unmethylated (converted)
              if(MDtag[i] %in% c("C", "G")){
                conversion[1,p] <- 21 } else{conversion[1,p] <- 19}
            }
        } else{next}
      }
    }
    return(conversion)
  }, double(length(cgsite)))

  #Remove reads without CpG sites
  rem <- colSums(is.na(conversion)) >= dim(conversion)[1]
  conversion <- conversion[,!rem]

  alns.pairs <- alns.pairs[names(alns.pairs) %in% colnames(conversion)]

  #### Get CpG of interest in window ####------------------------------------------------------------

  cpg.start <- start(cpgsite) - left + 1

  #### plot ####---------------------------------------------------------------------------------

  message("Plotting")

  #data for points
  d <- data.frame(CpG=rep(cgsite,length(alns.pairs)),
                  read=rep(1:length(alns.pairs), each=length(cgsite)),
                  value=as.vector(conversion))

  #data for segments
  reads <- d$read

  xstart <- sapply(reads, function(x){
    newd <- d[which(d$read == x),]
    vals <- which(!is.na(newd$value))
    min(newd$CpG[vals])
  })

  xend <- sapply(reads,function(x){
    newd <- d[which(d$read == x),]
    vals <- which(!is.na(newd$value))
    max(newd$CpG[vals])
  })

  d2 <- unique(data.frame(xstart, xend, reads))

  ggplot() +
    scale_shape_identity() +
    theme_void() +
    geom_segment(data=d2, aes(x=xstart, y=reads, xend=xend, yend=reads), colour = "grey", size = 0.5) +
    geom_point(data=d, aes(x=CpG, y=read, shape=value), fill = "white", size=1) +
    geom_point(aes(x=cpg.start, y = 0), shape = 24, size = 3, fill = "green")  +
    guides(color=FALSE)
}

