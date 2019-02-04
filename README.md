<!--- [![Build Status](https://travis-ci.com/csoneson/rnaseqworkflow.svg?branch=master)](https://travis-ci.com/csoneson/rnaseqworkflow) -->

# DAMEfinder

**DAMEfinder** (**D**ifferential **A**llele-specific **ME**thylation **finder**) is an R-package that detects allele-specific methylation (ASM) in a cohort of samples, and detects regions of differential ASM within groups of interest, based on **Bisulfite-sequencing** files.

DAMEfinder runs in two modes: **SNP-based** (exhaustive-mode) and **tuple-based** (fast-mode), which converge when calculating differential methylation.

![](vignettes/DAMEfinder_workflow.png)

Please refer to the vignette for more details on running the pipeline.

## What mode should I choose?

It depends on what you want to do and how much time you have. Either way you have to align your reads with [Bismark](https://github.com/FelixKrueger/Bismark) (apologies to other-aligner users).

### SNP-based

To run the **SNP-based** mode you need processed `bam` files *AND* a VCF file for each of your samples with heterozygous SNPs. I know this is typically not the case, so you could alternatively extract heterozygous SNPs using [BisSNP](https://github.com/dnaase/Bis-tools/tree/master/Bis-SNP) (which I have used), or [biscuit](https://github.com/zwdzwd/biscuit) (which I will test at some point).

I call this the "exhaustive-mode" because it extracts an ASM score for every CpG site in the reads containing each SNP from the VCF file. Based on this score DAMEs are detected.

From a biological point of view, you might want to run this mode if you are interested in loss or gain of allele-specificity linked to somatic heterozygous SNPs. More specifically, you could detect genes that exhibit loss of imprinting (e.g. [in colorectal cancer](http://cancerres.aacrjournals.org/content/62/22/6442.long)).

### tuple-based

To run the **tuple-based** mode you have to run [methtuple](https://github.com/PeteHaitch/methtuple) first. The methtuple output is the only thing needed for this mode. 

I call this the fast-mode because you don't need SNP information. The assumption is that intermediate levels of methylation represent ASM along the genome. For example, we have shown that the ASM score can distinguish females from males in the X chromosome. Using SNP information this wouldn't be possible.


## How do I install it?

Since DAMEfinder is not (yet) on Bioconductor, you have to install all dependencies before:

```{r}
## Install `BiocManager` if needed
if (!("BiocManager" %in% installed.packages()[, "Package"])) {
  install.packages("BiocManager")
}

## List dependencies
pkg <- c("BiocGenerics", "GenomeInfoDb", "GenomicRanges", "IRanges", 
"SummarizedExperiment", "limma", "bumphunter", "readr", 
"Rsamtools", "ggplot2",)

## Check if dependencies are already installed
pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]

## If some dependency is missing, install it
if (length(pkg) > 0) {
	BiocManager::install(pkg, dependencies = TRUE, ask = FALSE)
}
```

Now you can intall DAMEfinder

```{r}
BiocManager::install("markrobinsonuzh/DAMEfinder")
```
