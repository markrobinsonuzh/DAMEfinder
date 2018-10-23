# DAMEfinder
Finds DAMEs - Differential Allelicly MEthylated regions

To install, you'll need the following Bioconductor packages:
SummarizedExperiment, limma, bumphunter

These can be installed with the following commands:

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite(c("SummarizedExperiment","limma","bumphunter","readr","devtools", "Rsamtools"))
```

And, DAMEfinder can be installed (assuming you already have devtools) using:
```{r}
devtools::install_github("markrobinsonuzh/DAMEfinder", 
                         auth_token="f75f9e682dbbd16c118c9fc7963467bf8d1c60b6", 
                         dependencies = FALSE)
```

##run methtuple
python methtuple --sc --gzip -m 2 **.bam
python methtuple --sc --gzip -m 2 /home/Shared_sherborne/steph/B02_CRC_bismark_dedupl.bam
