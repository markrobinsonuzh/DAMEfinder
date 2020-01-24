## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(cache = 0, fig.width = 6, 
                      fig.height = 7)

## ----style, echo = FALSE, results = 'asis'------------------------------------
BiocStyle::markdown()

## -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(DAMEfinder)
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg19)
  })

bam_files <- c(system.file("extdata", "NORM1_chr19_trim.bam", 
                           package = "DAMEfinder"),
               system.file("extdata", "CRC1_chr19_trim.bam", 
                           package = "DAMEfinder"))

vcf_files <- c(system.file("extdata", "NORM1.chr19.trim.vcf", 
                           package = "DAMEfinder"),
               system.file("extdata", "CRC1.chr19.trim.vcf", 
                           package = "DAMEfinder"))

sample_names <- c("NORM1", "CRC1")

#Use another reference file for demonstration, and fix the seqnames
genome <- BSgenome.Hsapiens.UCSC.hg19
seqnames(genome) <- gsub("chr","",seqnames(genome))
reference_file <- DNAStringSet(genome[[19]], use.names = TRUE)
names(reference_file) <- 19

#Extract reads and extract methylation according to allele
snp.list <- extract_bams(bam_files, vcf_files, sample_names, reference_file,
                       coverage = 2)

#CpG sites for first SNP in VCF file from sample NORM1
snp.list$NORM1[[1]]

#CpG sites for first SNP in VCF file from sample CRC1
snp.list$CRC1[[1]]


## -----------------------------------------------------------------------------

derASM <- calc_derivedasm(snp.list)

derASM
assays(derASM)

## -----------------------------------------------------------------------------
x <- assay(derASM, "der.ASM")
head(x)

## -----------------------------------------------------------------------------

data(extractbams_output)

#The data loaded is an output from `split_bams()`, therefore we run 
#`calc_derivedasm` to get the SummarizedExperiment
derASM <- calc_derivedasm(extractbams_output, cores = 1, verbose = FALSE)

#We remove all CpG sites with any NA values, but not 0s
filt <- rowSums(!is.na(assay(derASM, "der.ASM"))) == 8 
derASM <- derASM[filt,]

#set the design matrix
grp <- factor(c(rep("CRC",4),rep("NORM",4)), levels = c("NORM", "CRC"))
mod <- model.matrix(~grp)
mod

#Run default
dames <- find_dames(derASM, mod)

head(dames)

#Run empirical method
dames <- find_dames(derASM, mod, pvalAssign = "empirical")

head(dames)


## -----------------------------------------------------------------------------
tuple_files <- c(system.file("extdata", "NORM1_chr19.qs.CG.2.tsv.gz", 
                             package = "DAMEfinder"),
                 system.file("extdata", "CRC1_chr19.qs.CG.2.tsv.gz", 
                             package = "DAMEfinder"))

sample_names <- c("NORM1", "CRC1")

tuple_list <- read_tuples(tuple_files, sample_names)

head(tuple_list$NORM1)

## -----------------------------------------------------------------------------

ASM_mat <- calc_asm(tuple_list)
ASM_mat

## -----------------------------------------------------------------------------
#load package data
data(readtuples_output)

#run calc_asm and filter object
ASMscore <- calc_asm(readtuples_output)
filt <- rowSums(!is.na(assay(ASMscore, "asm"))) == 5 #filt to avoid warnings
ASMscore <- ASMscore[filt,]

#make design matrix (or specify a contrast)
grp <- factor(c(rep("CRC",3),rep("NORM",2)), levels = c("NORM", "CRC"))
mod <- model.matrix(~grp)

#run default and increase maxGap to get longer, more sparse regions
dames <- find_dames(ASMscore, mod, maxGap = 300)

head(dames)

#run alternative mode
dames <- find_dames(ASMscore, mod,  maxGap = 300, pvalAssign = "empirical")

head(dames)


## ---- dametrack---------------------------------------------------------------

#Here I will use the tuple-ASM SummExp
colData(ASMscore)$group <- grp
colData(ASMscore)$samples <- colnames(ASMscore)

#Set a DAME as a GRanges. I choose a one from the tables we obtained above
dame <- GRanges(19,IRanges(323736,324622))

dame_track(dame = dame,
           ASM = ASMscore)

## ---- dt2---------------------------------------------------------------------
dame_track(dame = dame,
           ASM = ASMscore,
           window = 2)


## ---- dt3---------------------------------------------------------------------

dame <- GRanges(19,IRanges(387966,387983))

grp <- factor(c(rep("CRC",4),rep("NORM",4)), levels = c("NORM", "CRC"))
colData(derASM)$group <- grp

dame_track(dame = dame,
           derASM = derASM)


## ---- dt4---------------------------------------------------------------------
dame_track(dame = dame,
           derASM = derASM,
           plotSNP = TRUE)


## ---- dt5---------------------------------------------------------------------
dame_track(dame = dame,
           derASM = derASM,
           ASM = ASMscore)


## ---- dt6---------------------------------------------------------------------
dame_track_mean(dame = dame,
           derASM = derASM,
           ASM = ASMscore)


## ---- fig1--------------------------------------------------------------------
#put SNP in GRanges (you can find the SNP with the dame_track function)
snp <- GRanges(19, IRanges(267039, width = 1)) #always set the width if your 
#GRanges has 1 site

snp

bam.file <- system.file("extdata", "CRC1_chr19_trim.bam", 
                        package = "DAMEfinder")

vcf.file <- system.file("extdata", "CRC1.chr19.trim.vcf", 
                        package = "DAMEfinder")

methyl_circle_plot(snp = snp, vcfFile = vcf.file, bamFile = bam.file, 
                   refFile = reference_file)

## ---- fig2--------------------------------------------------------------------

cpgsite <- GRanges(19, IRanges(266998, width = 1))

methyl_circle_plot(snp = snp, vcfFile = vcf.file, bamFile = bam.file, 
                   refFile = reference_file, cpgsite = cpgsite)

## ---- fig3--------------------------------------------------------------------

cpgsite <- GRanges(19, IRanges(266998, width = 1))

methyl_circle_plotCpG(cpgsite = cpgsite, bamFile = bam.file, 
                      refFile = reference_file)

## ---- fig4--------------------------------------------------------------------

#a random region
dame <- GRanges(19, IRanges(266998,267100))

methyl_circle_plot(snp = snp, vcfFile = vcf.file, bamFile = bam.file, 
                   refFile = reference_file, dame = dame)


## ---- fig5--------------------------------------------------------------------

grp <- factor(c(rep("CRC",3),rep("NORM",2)), levels = c("NORM", "CRC"))
methyl_MDS_plot(ASMscore, group = grp)


## ---- fig6, eval=FALSE--------------------------------------------------------
#  
#  #path to chrom.sizes file
#  chromsizes <- system.file("extdata", "hg19.chrom.sizes.mod",
#                            package = "DAMEfinder")
#  
#  #For a single sample
#  make_bigwig(samples = "NORM1", scoreObj = derASM, folder = "",
#              chromsizes.file = chromsizes)
#  
#  #For all samples in object
#  vapply(colnames(derASM), make_bigwig, scoreObj = derASM, folder = "",
#         chromsizesFile = chromsizes)
#  

