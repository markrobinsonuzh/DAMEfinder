## ----setup, include=FALSE--------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----style, echo = FALSE, results = 'asis'---------------------------------
BiocStyle::markdown()

## ----eval=FALSE------------------------------------------------------------
#  library(DAMEfinder)
#  
#  file.copy(system.file(package="DAMEfinder", "extdata"), ".", recursive=TRUE)
#  meta_data <- read.table("extdata/metadata.txt", header = TRUE)
#  file_list <- list.files(path = "extdata/", pattern = ".tsv.gz")
#  m <- match(file_list, meta_data$methtuplefile)
#  file_list <- paste0("extdata/", file_list)
#  names(file_list) <- meta_data$shortname[m]
#  
#  file_list
#  

## ----eval=FALSE------------------------------------------------------------
#  
#  sample_list <- read_tuples(files = file_list, sample_names = names(file_list))
#  
#  head(sample_list$"adenoma1")
#  

## ----eval=FALSE------------------------------------------------------------
#  
#  sample_list_no_snps <- parallel::mclapply(sample_list, remove_snps, mc.cores=6)

## ----eval=FALSE------------------------------------------------------------
#  
#  ASM_score_matrix <- calc_asm(sample_list = sample_list)
#  

## ----eval=FALSE------------------------------------------------------------
#  
#  # Make a design matrix that specifies the two conditions per sample
#  cols <- rep(1,ncol(ASM_score_matrix))
#  n <- grep("normal", colnames(ASM_score_matrix))
#  cols[n] <- 0
#  
#  mod <- matrix(data=c(rep(1,ncol(ASM_score_matrix)), cols), ncol = 2)
#  mod
#  
#  # Get t-Statistics
#  tStatistics <- get_tstats(ASM_score_matrix, mod, method = "ls")
#  

## ----eval=FALSE------------------------------------------------------------
#  
#  dames <- find_dames(tStatistics)
#  
#  head(dames)
#  

