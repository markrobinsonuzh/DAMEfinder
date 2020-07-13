#' DAMEfinder: Method to detect allele-specific methylation (ASM), and
#' differential ASM from Bisulfite sequencing data in R.
#'
#' The package allows the user to extract an ASM score in two ways: either from
#' a \code{bismark} bam file(s) and VCF file(s), or from the output from
#' \code{methtuple}. Either way the final output is a list of regions with
#' diferential allele-specific methylated between groups of samples of interest.
#' The package also provides functions to visualize ASM at the read level or the
#' score level
#' 
#' @section DAMEfinder functions: \code{calc_asm} extracts ASM for pairs of CpG
#'   sites from a methtuple file, \code{calc_derivedasm} extracts ASM at each
#'   CpG site linked to a SNP from the VCF file. Both functions generate a
#'   \code{RangedSummarizedExperiment}, which is the input for the main function
#'   \code{find_dames}, that generates a \code{data.frame} with regions
#'   exhibiting differential ASM between a number of samples.
#' 
#' @author Stephany Orjuela \email{sorjuelal@@gmail.com}
#' @author Dania Machlab
#' @author Mark D Robinson \email{mark.robinson@@imls.uzh.ch}
#' @name DAMEfinder
#' @docType package
NULL