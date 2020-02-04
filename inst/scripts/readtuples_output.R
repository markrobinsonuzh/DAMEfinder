library(SummarizedExperiment)


DATA_PATH_DIR <- "inst/extdata/moredata"

get_data_path <- function(file_name) file.path(DATA_PATH_DIR, file_name)

tuple_files <- list.files(DATA_PATH_DIR, ".tsv.gz")
tuple_files <- get_data_path(tuple_files)

readtuples_output <- read_tuples(tuple_files, 
                                 c("CRC1", "CRC2", "CRC3", "NORM1", "NORM3"))


usethis::use_data(readtuples_output, overwrite = TRUE, compress = 'xz')
