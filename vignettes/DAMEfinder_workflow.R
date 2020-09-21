## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(cache = 0, fig.width = 6, 
                      fig.height = 7)


## ----style, echo = FALSE, results = 'asis'------------------------------------
BiocStyle::markdown()


## ----eval=FALSE---------------------------------------------------------------
## if (!requireNamespace("BiocManager", quietly = TRUE))
##     install.packages("BiocManager")
## 
## BiocManager::install("DAMEfinder")


## #Check quality of reads

## fastqc -t 2  sample1_R1.fastq.gz sample1_R2.fastq.gz

## 
## #Trim reads to remove bad quality regions and adapter sequence

## trim_galore --paired sample1_R1.fastq.gz sample2_R2.fastq.gz


## #Build bisulfite reference

## bismark_genome_preparation <path_to_genome_folder>

## 
## #run Bismark

## bismark -B sample1 --genome <path_to_genome_folder>

##     -1 sample1_R1_val_1.fq.gz

##     -2 sample1_R2_val_2.fq.gz

## 
## #deduplicate (optional)

## deduplicate_bismark -p --bam sample1_pe.bam

## 
## #sort and index files

## samtools sort -m 20G -O bam -T _tmp

##     -o sample1_pe.dedupl_s.bam sample1_pe.deduplicated.bam

## samtools index file1_pe.dedupl_s.bam


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



## 
## # Sort bam file by query name

## samtools sort -n -@ 10 -m 20G -O bam -T _tmp

##     -o sample1_pe_sorted.bam sample1_pe.deduplicated.bam

## 
## # Run methtuple

## methtuple --sc --gzip -m 2 sample1_pe_sorted.bam


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



## -----------------------------------------------------------------------------
utils::sessionInfo()

