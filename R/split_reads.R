#' Divide read names by allele
#'
#' Takes a GenomicAlignments object and returns a list of read
#' names dividied by allele
#'
#' @param alns GenomicAlignments object.
#' @param snp GRanges object containing SNP location.
#' @param v Nucleotide of reference (or alternative) allele
#'
#' @return A named list of vectors, each vector containing read names for each allele
#' @examples
#' @importFrom BiocGenerics start
#' @importFrom BiocGenerics end
#'
#'
#' @export
splitReads <- function(alns, v, snp){

  #Use MD tag from bam to check wich reads have the ref or alt allele
  alleles <- vapply(alns, function(x){

    #extract matches and mismatches from MD tag
    tag <- getMD(x)
    MDtag <- tag$MDtag
    nucl.num <- tag$nucl.num

    #Get right location of snp in sequence of read
    snp.start <- start(snp) - start(x) + 1

    #Get the snp (location) from the specific read
    count <- 0

    for(i in seq_along(nucl.num)){
      count <- count + nucl.num[i]
      if(count >= snp.start){break}
    }
    return(MDtag[i])
  }, character(1))

  #Get read names from alt and ref reads
  ref.reads <- sort(names(which(alleles != v)))
  alt.reads <- sort(names(which(alleles == v)))

  return(list(ref.reads = unique(ref.reads), alt.reads = unique(alt.reads)))
}


#' MDtag parser
#'
#' Takes a GenomicAlignments object containing the MDtag, and transforms it
#' into a vector of characters and numbers
#'
#' @param y GenomicAlignments object containing an MDtag
#'
#' @return A named list of vectors, each vector a parsed version of MDtag:
#' - nucl.num: Numeric representation of MDtag.
#' - MDtag: a split version of MDtag
#'
#' @examples
#'
#' @importFrom S4Vectors mcols
#'
#'
#' @export
getMD <- function(y){

  a <- mcols(y)$MD

  #extract matches and mismatches from MD tag
  numbers <- as.integer(stringr::str_extract_all(a, "[0-9]{1,}")[[1]])
  nucl <- stringr::str_extract_all(a, "[A-Z^]{1,2}")[[1]]
  MDtag <- c(numbers, nucl)[order(c(seq_along(numbers)*2 - 1, seq_along(nucl)*2))]

  if(length(MDtag) == 1){
    nucl.num <- as.integer(MDtag)
  } else{
    nucl.num <- numeric(length(numbers)+length(nucl))
    nucl.num[seq(2,length(nucl.num), 2)] <- 1
    nucl.num[nucl.num != 1] <- numbers
  }
  return(list(MDtag = MDtag, nucl.num = nucl.num))
}



