## ----setup, include=FALSE--------------------------------------------------
knitr::opts_chunk$set(cache = 0, fig.path = "figure/", fig.width = 6, fig.height = 7)

## ----style, echo = FALSE, results = 'asis'---------------------------------
BiocStyle::markdown()

## ---- echo=FALSE, out.width='250%', fig.align='center'---------------------
knitr::include_graphics('DAMEfinder_workflow.png')

## --------------------------------------------------------------------------
library(DAMEfinder)

bam_files <- c(system.file("extdata", "NORM1_chr19.bam", package = "DAMEfinder"),
               system.file("extdata", "CRC1_chr19.bam", package = "DAMEfinder"))

vcf_files <- c(system.file("extdata", "NORM1.chr19.moretrim.vcf", package = "DAMEfinder"),
               system.file("extdata", "CRC1.chr19.moretrim.vcf", package = "DAMEfinder"))

sample_names <- c("NORM1", "CRC1")

reference_file <- system.file("extdata", "19.fa", package = "DAMEfinder")

#Split reads and extract methylation according to allele
split_bams(bam_files, vcf_files, sample_names, reference_file)

#Read in one of the generated files
snp.table <- readRDS("snp.table.NORM1.rds")

#CpG sites for first SNP in VCF file
snp.table[[1]]

#CpG sites for second SNP in VCF file
snp.table[[2]]

#And so on...

## --------------------------------------------------------------------------

rds_files <- list("snp.table.NORM1.rds", "snp.table.CRC1.rds")
derASM <- calc_derivedasm(rds_files)

derASM

## --------------------------------------------------------------------------
suppressPackageStartupMessages(library(SummarizedExperiment))

head(assays(derASM)[["der.ASM"]])


## ---- eval=FALSE-----------------------------------------------------------
#  
#  # Make a design matrix that specifies the two conditions per sample
#  cols <- rep(1,ncol(derASM))
#  n <- grep("NORM", colnames(derASM))
#  cols[n] <- 0
#  
#  mod <- matrix(data=c(rep(1,ncol(derASM)), cols), ncol = 2)
#  mod
#  
#  # Get t-Statistics
#  tStatistics <- get_tstats(derASM, mod, method = "ls")
#  

## ---- eval=FALSE-----------------------------------------------------------
#  
#  dames <- find_dames(tStatistics, Q=0.9, maxGap=20)
#  
#  head(dames)
#  

## --------------------------------------------------------------------------
tuple_files <- c(system.file("extdata", "NORM1_chr19_qs.CG.2.tsv.gz", package = "DAMEfinder"),
                 system.file("extdata", "CRC1_chr19_qs.CG.2.tsv.gz", package = "DAMEfinder"))


tuple_list <- read_tuples(files = tuple_files, sample_names)

head(tuple_list$NORM1)

## --------------------------------------------------------------------------

ASM_mat <- calc_asm(sample_list = tuple_list)
ASM_mat

## ---- eval=FALSE-----------------------------------------------------------
#  #sample_list_no_snps <- parallel::mclapply(sample_list, remove_snps, mc.cores=6)

## ---- eval=FALSE-----------------------------------------------------------
#  
#  tStatistics <- get_tstats(ASm_mat)
#  
#  dames <- find_dames(tStatistics, Q=0.9, maxGap=20)
#  
#  head(dames)

## ---- fig1-----------------------------------------------------------------
library(GenomicRanges)
snp <- GRanges(19, IRanges(267039, width = 1))

snp

bam.file <- system.file("extdata", "NORM1_chr19.bam", package = "DAMEfinder")

vcf.file <- system.file("extdata", "NORM1.chr19.moretrim.vcf", package = "DAMEfinder")

ref.file <- system.file("extdata", "19.fa", package = "DAMEfinder")

methyl_circle_plot(snp = snp, vcf.file = vcf.file, bam.file = bam.file, ref.file = ref.file, sample.name = "sample1")

## ---- fig2-----------------------------------------------------------------

cpgsite <- GRanges(19, IRanges(266998, width = 1))

methyl_circle_plot(snp = snp, vcf.file = vcf.file, bam.file = bam.file, ref.file = ref.file, sample.name = "sample1", cpgsite = cpgsite)

## ---- fig3-----------------------------------------------------------------

cpgsite <- GRanges(19, IRanges(266998, width = 1))

methyl_circle_plotCpG(cpgsite = cpgsite, bam.file = bam.file, ref.file = ref.file)

