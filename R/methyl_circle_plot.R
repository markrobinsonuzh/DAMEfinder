#' Draw methylation circle plot
#'
#' Draws CpG site methylation status as points, in reads containing a specific
#' SNP. Generates one plot per bam file.
#'
#' @param snp GRanges object containing SNP location.
#' @param cpgsite (optional) GRanges object containing a single CpG site
#'   location of interest.
#' @param vcfFile vcf file.
#' @param bamFile bismark bam file path.
#' @param refFile fasta reference file path. Or \code{DNAStringSet} with DNA
#'  sequence.
#' @param build genome build used. default = "hg19"
#' @param dame (optional) GRanges object containing a region to plot.
#' @param letterSize Size of alleles drawn in plot. Default = 2.5.
#' @param pointSize Size of methylation circles. Default = 3.
#' @param sampleName FIX?: this is to save the vcf file to not generate it every
#'   time you run the function.
#' @param sampleReads Whether a subset of reads should be plotted.
#'   Default = FALSE.
#' @param numReads Number of reads to plot per allele, if sampleReads is TRUE.
#'   Default = 20
#'
#' @return Plot
#' @examples
#' DATA_PATH_DIR <- system.file('extdata', '.', package = 'DAMEfinder')
#'
#' get_data_path <- function(file_name) file.path(DATA_PATH_DIR, file_name)
#' bam_files <- get_data_path('NORM1_chr19_trim.bam')
#' vcf_files <- get_data_path('NORM1.chr19.trim.vcf')
#' sample_names <- 'NORM1'
#' #reference_file
#' suppressPackageStartupMessages({library(BSgenome.Hsapiens.UCSC.hg19)})
#' genome <- BSgenome.Hsapiens.UCSC.hg19
#' seqnames(genome) <- gsub("chr","",seqnames(genome))
#' dna <- DNAStringSet(genome[[19]], use.names = TRUE)
#' names(dna) <- 19
#'
#' snp <- GenomicRanges::GRanges(19, IRanges::IRanges(292082, width = 1))
#' methyl_circle_plot(snp = snp,
#'  vcfFile = vcf_files,
#'  bamFile = bam_files,
#'  refFile = dna,
#'  sampleName = sample_names)
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics start
#' @importFrom BiocGenerics end
#' @importFrom GenomeInfoDb seqnames
#' @import ggplot2
#'
#' @export

methyl_circle_plot <- function(snp, vcfFile, bamFile, refFile, build = "hg19",
    dame = NULL, letterSize = 2.5, pointSize = 3, sampleName = "sample1", 
    cpgsite = NULL, sampleReads = FALSE, numReads = 20) {
    
    message("Reading vcf file")
    if (!file.exists(paste0(sampleName, ".RData"))) {
      param <- VariantAnnotation::ScanVcfParam(fixed="ALT", info=NA, geno=NA)  
      vcf <- VariantAnnotation::readVcf(vcfFile, genome = build, 
                                        param = param)
        save(vcf, file = paste0(sampleName, ".RData"))
    } else {
        load(paste0(sampleName, ".RData"))
    }
    
    snp.loc <- which(start(vcf) == start(snp) &  
                       as.character(seqnames(vcf)) == levels(seqnames(snp)))
    if (length(snp.loc) == 0) {
        message("Sample does not contain SNP")
        return(NULL)
        stop()
    }
    ref <- as.character(VariantAnnotation::ref(vcf))[snp.loc]
    alt <- as.character(
      BiocGenerics::unlist(VariantAnnotation::alt(vcf)))[snp.loc]
    
    message("Getting reads")
    alns.pairs <- GenomicAlignments::readGAlignmentPairs(bamFile, 
        param = Rsamtools::ScanBamParam(tag = c("MD", "XM", "XR", 
            "XG"), which = snp), use.names = TRUE)
    
    alns <- unlist(alns.pairs)  #unpaired
    
    split <- splitReads(alns, ref, snp)
    
    alt.reads <- split$alt.reads
    ref.reads <- split$ref.reads
    
    if (sampleReads) {
        alt.reads <- sample(alt.reads, numReads)
        ref.reads <- sample(ref.reads, numReads)
    }
    
    message("Sorting reads")
    # don't unlist, paired
    alt.pairs <- alns.pairs[names(alns.pairs) %in% alt.reads]
    ref.pairs <- alns.pairs[names(alns.pairs) %in% ref.reads]
    alns.pairs <- c(alt.pairs, ref.pairs)  #I just do this to sort them
    
    #### get reference and CpG
    #### positions####------------------------------------
    
    # Get limits for plotting
    
    if (is.null(dame)) {
        left <- min(start(alns))
        right <- max(end(alns))
        window <- GRanges(seqnames(snp), IRanges(left, right))
    } else {
        window <- dame
        left <- start(dame)
        right <- end(dame)
    }
    
    message("Reading reference")
    # open reference seq to get correct CpG locations within that
    # reference chunk
    if (typeof(refFile) == "character") {
        fa <- open(Rsamtools::FaFile(refFile, index = paste0(refFile, 
            ".fai")))
        dna <- Rsamtools::scanFa(fa, param = window)
    } else {
        dna <- refFile[window]
    }
    
    
    cgsite <- stringr::str_locate_all(dna, "CG")[[1]][, 1]  #also look at GpCs?
    
    if (length(cgsite) < 1) {
        return(NULL)
        stop("No CpG sites associated to this SNP")
    }
    
    mepos <- cgsite/Biostrings::nchar(dna)  #location of CpG sites
    
    # Use MD tag from bam to extract methylation status for each
    # site detected above
    message("Getting meth state per read-pair")
    conversion <- vapply(alns.pairs, function(x) {
        
        # change C locations for G locations if reference context is
        # different
        if (S4Vectors::mcols(x)$XG[1] == "GA") {
            cgsite <- cgsite + 1
        }
        
        # fuse pair info in one row
        conversion <- matrix(NA, 1, length(cgsite))
        
        # a <- mcols(x)$MD
        read.start <- start(x) - left + 1  #start of read
        read.end <- end(x) - left + 1  #end
        
        for (pair in c(1, 2)) {
            
            something <- mcols(x[pair])$MD
            tag <- getMD(something)
            MDtag <- tag$MDtag
            nucl.num <- tag$nucl.num
            
            # Get the locations from the specific read: Count along the
            # numerical MDtag until you reach the actual CpG start
            count <- 0
            this.mepos <- cgsite - read.start[pair] + 1
            
            # Fill in conversion table with meth state for each read 21
            # is unmethylated 19 is methylated 0 is not in read
            
            for (p in seq_along(this.mepos)) {
                if (is.na(conversion[1, p])) {
                    if (cgsite[p] > read.end[pair] || cgsite[p] < 
                    read.start[pair]) {
                    conversion[1, p] <- NA
                    next
                    } else {
                    for (i in seq_along(nucl.num)) {
                        count <- count + nucl.num[i]
                        if (count >= this.mepos[p]) {
                        count <- 0
                        break
                        }
                    }
                    
                    # if there is a mismatch, the MDtag shows what the base is 
                    # in the ref. So, if it shows a C or G, it means that the 
                    # ref is C and the read has T or A, therefore it is 
                    # unmethylated (converted)
                    if (MDtag[i] %in% c("C", "G")) {
                        conversion[1, p] <- 21
                    } else {
                        conversion[1, p] <- 19
                    }
                        }
                } else {
                    next
                }
            }
        }
        return(conversion)
    }, double(length(cgsite)))
    
    # Remove reads without CpG sites
    rem <- BiocGenerics::colSums(is.na(conversion)) >= dim(conversion)[1]
    conversion <- conversion[, !rem]
    
    alt.reads <- alt.reads[alt.reads %in% colnames(conversion)]
    ref.reads <- ref.reads[ref.reads %in% colnames(conversion)]
    alns.pairs <- alns.pairs[names(alns.pairs) %in% colnames(conversion)]
    
    #### Get snp position in window
    #### ####----------------------------------------
    
    letter <- c(rep(alt, length(alt.reads)), rep(ref, length(ref.reads)))
    
    snp.start <- start(snp) - left + 1
    
    # and cpg of interest
    if (!is.null(cpgsite)) {
        cpg.start <- start(cpgsite) - left + 1
    } else {
        cpg.start = snp.start
    }
    
    #### plot
    #### ####-------------------------------------------------------------
    
    message("Plotting")
    
    # data for points
    d <- data.frame(CpG = rep(cgsite, length(alns.pairs)), 
        read = rep(seq(from = 1, 
        to = length(alns.pairs), by = 1), each = length(cgsite)), 
        value = as.vector(conversion), stringsAsFactors = FALSE)
    
    # data for segments
    reads <- d$read
    
    xstart <- vapply(reads, function(x) {
        newd <- d[which(d$read == x), ]
        vals <- which(!is.na(newd$value))
        min(newd$CpG[vals], snp.start)
    }, FUN.VALUE = double(1))
    
    xend <- vapply(reads, function(x) {
        newd <- d[which(d$read == x), ]
        vals <- which(!is.na(newd$value))
        max(newd$CpG[vals], snp.start)
    }, FUN.VALUE = double(1))
    
    d2 <- unique(data.frame(xstart, xend, reads, stringsAsFactors = FALSE))
    d2$snp <- c(rep("a", length(alt.reads)), rep("r", length(ref.reads)))
    
    # To manually scale the colors
    cols <- c("#0E4071", "#d55e00", "#0E4071", "#d55e00")
    names(cols) <- c("a", "r", alt, ref)
    
    mp <- ggplot() + 
        scale_shape_identity() + 
        theme_void() + 
        geom_segment(data = d2, aes(x = xstart, y = reads, xend = xend, 
            yend = reads, colour = snp), size = 0.2) + 
        geom_point(data = d, aes_(x = ~CpG, y = ~read, shape = ~value), 
            fill = "white", size = pointSize) + 
        geom_point(aes(x = snp.start, y = seq(from = 1, to = length(alns.pairs),
            by = 1), shape = letter, colour = letter), size = letterSize) + 
        geom_point(aes(x = cpg.start, y = 0), shape = 24, size = 3, 
            fill = "green") + 
        scale_color_manual(values = cols) + 
        guides(color = FALSE) + 
        ggtitle(sampleName)
    
    return(mp)
}

#' Draw methylation circle plot without SNP
#'
#' Draws CpG site methylation status as points, in reads containing a specific
#' CpG site. Generates one plot per bam file.
#'
#' @param cpgsite GRanges object containing a single CpG site location of
#'   interest
#' @param bamFile bismark bam file path
#' @param refFile fasta reference file path
#' @param pointSize Size of methylation circles. Default = 3.
#' @param dame (optional) GRanges object containing a region to plot
#' @param order Whether reads should be sorted by methylation status. Default=
#' False.
#' @param sampleName Plot title.
#' @param sampleReads Whether a subset of reads should be plotted.
#'   Default = FALSE.
#' @param numReads Number of reads to plot, if sampleReads is TRUE. Default = 20
#' @return Plot
#' @examples
#' DATA_PATH_DIR <- system.file('extdata', '.', package = 'DAMEfinder')
#' get_data_path <- function(file_name) file.path(DATA_PATH_DIR, file_name)
#' bam_files <- get_data_path('NORM1_chr19_trim.bam')
#' sample_names <- 'NORM1'
#' #reference_file
#' suppressPackageStartupMessages({library(BSgenome.Hsapiens.UCSC.hg19)})
#' genome <- BSgenome.Hsapiens.UCSC.hg19
#' seqnames(genome) <- gsub("chr","",seqnames(genome))
#' dna <- DNAStringSet(genome[[19]], use.names = TRUE)
#' names(dna) <- 19
#'
#' cpg <- GenomicRanges::GRanges(19, IRanges::IRanges(292082, width = 1))
#' methyl_circle_plotCpG(cpgsite = cpg,
#'  bamFile = bam_files,
#'  refFile = dna)
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics start
#' @importFrom BiocGenerics end
#' @importFrom GenomeInfoDb seqnames
#' @import ggplot2
#'
#' @export
methyl_circle_plotCpG <- function(cpgsite = cpgsite, bamFile = bamFile, 
    pointSize = 3, refFile = refFile, dame = NULL, order = FALSE, 
    sampleName = NULL, sampleReads = FALSE, numReads = 20) {
    
    alns.pairs <- GenomicAlignments::readGAlignmentPairs(bamFile, 
        param = Rsamtools::ScanBamParam(tag = c("MD", "XM", "XR", 
            "XG"), which = cpgsite), use.names = TRUE)
    
    
    if (sampleReads) {
        ran.names <- sample(names(alns.pairs), numReads)
        alns.pairs <- alns.pairs[names(alns.pairs) %in% ran.names]
    }
    
    alns <- unlist(alns.pairs)
    
    #### get reference and CpG
    #### positions####-------------------------------------
    
    # Get limits for plotting
    
    if (is.null(dame)) {
        left <- min(start(alns))
        right <- max(end(alns))
        window <- GRanges(seqnames(cpgsite), IRanges(left, right))
    } else {
        window <- dame
        left <- start(dame)
        right <- end(dame)
    }
    
    message("Reading reference")
    # open reference seq to get correct CpG locations within that
    # reference chunk
    if (typeof(refFile) == "character") {
        fa <- open(Rsamtools::FaFile(refFile, index = paste0(refFile, 
            ".fai")))
        dna <- Rsamtools::scanFa(fa, param = window)
    } else {
        dna <- refFile[window]
    }
    
    cgsite <- stringr::str_locate_all(dna, "CG")[[1]][, 1]  #also look at GpCs?
    
    if (length(cgsite) < 1) {
        stop("No CpG sites in these reads")
    }
    
    mepos <- cgsite/Biostrings::nchar(dna)  #location of CpG sites
    
    ##### Use MD tag from bam to extract methylation status for each
    ##### site detected above ####----------------------------
    message("Getting meth state per read-pair")
    conversion <- vapply(alns.pairs, function(x) {
        
        # change C locations for G locations if reference context is
        # different
        if (S4Vectors::mcols(x)$XG[1] == "GA") {
            cgsite <- cgsite + 1
        }
        
        # fuse pair info in one row
        conversion <- matrix(NA, 1, length(cgsite))
        
        # a <- mcols(x)$MD
        read.start <- start(x) - left + 1  #start of read
        read.end <- end(x) - left + 1  #end
        
        for (pair in c(1, 2)) {
            
            something <- mcols(x[pair])$MD
            tag <- getMD(something)
            MDtag <- tag$MDtag
            nucl.num <- tag$nucl.num
            
            # Get the locations from the specific read: Count along the
            # numerical MDtag until you reach the actual CpG start
            count <- 0
            this.mepos <- cgsite - read.start[pair] + 1
            
            # Fill in conversion table with meth state for each read 21
            # is unmethylated 19 is methylated 0 is not in read
            
            for (p in seq_along(this.mepos)) {
                if (is.na(conversion[1, p])) {
                    if (cgsite[p] > read.end[pair] || cgsite[p] < 
                    read.start[pair]) {
                    conversion[1, p] <- NA
                    next
                    } else {
                    for (i in seq_along(nucl.num)) {
                        count <- count + nucl.num[i]
                        if (count >= this.mepos[p]) {
                            count <- 0
                            break
                        }
                    }
                    
                    # if there is a mismatch, the MDtag shows what the base is 
                    # in the ref. So, if it shows a C or G, it means that the 
                    # ref is C and the read has T or A, therefore it is 
                    # unmethylated (converted)
                    if (MDtag[i] %in% c("C", "G")) {
                        conversion[1, p] <- 21
                    } else {
                        conversion[1, p] <- 19
                    }
                    }
                } else {
                    next
                }
            }
        }
        return(conversion)
    }, double(length(cgsite)))
    
    # Remove reads without CpG sites
    rem <- colSums(is.na(conversion)) >= dim(conversion)[1]
    conversion <- conversion[, !rem]
    
    alns.pairs <- alns.pairs[names(alns.pairs) %in% colnames(conversion)]
    
    #### Get CpG of interest in window
    #### ####------------------------------------
    
    cpg.start <- start(cpgsite) - left + 1
    
    #### plot
    #### ####-------------------------------------------------------------
    
    message("Plotting")
    
    # data for points
    d <- data.frame(CpG = rep(cgsite, length(alns.pairs)), 
        read = rep(seq(from = 1, 
        to = length(alns.pairs), by = 1), each = length(cgsite)), 
        value = as.vector(conversion))
    
    # reorder reads to plot
    if (isTRUE(order)) {
        trans <- ifelse(d$value == 19, 1, 0)
        und <- unique(d$read)
        sums_per_read <- vapply(und, function(i) sum(trans[d$read == 
            i], na.rm = TRUE))
        und_ord <- und[order(sums_per_read)]
        
        d$order <- 0
        for (i in und) {
            d$order[d$read == i] <- which(und_ord == i)
        }
        d$read <- d$order
    }
    
    # data for segments
    reads <- d$read
    
    xstart <- vapply(reads, function(x) {
        newd <- d[which(d$read == x), ]
        vals <- which(!is.na(newd$value))
        min(newd$CpG[vals])
    }, FUN.VALUE = double(1))
    
    xend <- vapply(reads, function(x) {
        newd <- d[which(d$read == x), ]
        vals <- which(!is.na(newd$value))
        max(newd$CpG[vals])
    }, FUN.VALUE = double(1))
    
    d2 <- unique(data.frame(xstart, xend, reads))
    
    ggplot() + 
        scale_shape_identity() + 
        theme_void() + 
        geom_segment(data = d2, aes(x = xstart, y = reads, xend = xend,
            yend = reads), colour = "grey", size = 0.5) + 
        geom_point(data = d, aes_(x = ~CpG, y = ~read, shape = ~value), 
            fill = "white", size = pointSize) + 
        geom_point(aes(x = cpg.start, y = 0), shape = 24, size = 3, 
            fill = "green") + 
        guides(color = FALSE) + 
        ggtitle(sampleName)
    
}

