# DAMEfinder
Finds DAMEs - Differential Allelicly MEthylated regions

To install, you'll need the following Bioconductor packages:
SummarizedExperiment, limma, bumphunter

These can be installed with the following commands:

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite(c("SummarizedExperiment","limma","bumphunter"))
```

And, DAMEfinder can be installing (assuming you already have devtools) using:
```{r}
devtools::install_github("markrobinsonuzh/DAMEfinder")
```
