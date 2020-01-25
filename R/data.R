#' extract_bams() output.
#'
#' 4 Patients from a previous study (Parker et al, 2018.) with colorectal cancer
#' were sequenced and the normal and cancerous tissue of each patient was
#' obtained. The data includes a subset of chromosome 19.
#' 
#' @format A large list with 8 elements. Each element is a list of
#'   \code{GRanges} for each sample. Each GRanges in the list includes the
#'   location of the CpG sites contained in the reads for each SNP. The GRanges
#'   metadata table contains: 
#'   \describe{
#' \item{\code{cov.ref}}{Number of reads of "reference" allele in that SNP}
#' \item{\code{cov.alt}}{Number of reads of "alternative" allele in that SNP}
#' \item{\code{meth.ref}}{Number of methylated reads of "reference" allele in
#' that SNP}
#' \item{\code{cov.ref}}{Number of methylated reads of "alternative" allele in
#' that SNP}
#' \item{\code{snp}}{The SNP containing the reads}
#' }
#' 
#' For further details, 
#' see \url{https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6949/}
#' sample names in in ArrayExpress do not necessarily match names given here!
"extractbams_output"

#' read_tuples() output.
#'
#' 3 Patients from a previous study (Parker et al, 2018.) with colorectal cancer
#' were sequenced and the normal and cancerous tissue of each patient was
#' obtained. The data includes a subset of chromosome 19. Here one normal sample
#' is not included.
#' 
#' @format A large list with 5 elements. Each element is a \code{tibble} with
#'   the coordinates of the pairs of CpG sites (tuples). Rest of the tibble
#'   contains:
#'   \describe{
#' \item{\code{MM}}{Number of reads with both CpG sites methylated}
#' \item{\code{MU}}{Number of reads with first CpG site methylated}
#' \item{\code{UM}}{Number of reads with second CpG site methylated}
#' \item{\code{UU}}{Number of reads with both CpG sites unmethylated}
#' \item{\code{cov}}{Coverage, total reads at tuple}
#' \item{\code{inter_dist}}{Distance in bp between CpG sites}
#' }
#' 
#' For further details, 
#' see \url{https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6949/}
#' sample names in in ArrayExpress do not necessarily match names given here!
"readtuples_output"

