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
#' @importFrom S4Vectors mcols
#' @importFrom stringr str_extract_all
#' @importFrom stringr str_extract
#'
#' @export
splitReads <- function(alns, v, snp){

  fullMD <- mcols(alns)$MD
  fullstart <- start(alns)
  snp.loc <- start(snp)

  #Use MD tag from bam to check which reads have the ref or alt allele
  alleles <- mapply(function(a, x, snp.loc){

    #extract matches and mismatches from MD tag
    tag <- getMD(a)
    MDtag <- tag$MDtag
    nucl.num <- tag$nucl.num

    #Get right location of snp in sequence of read
    snp.start <- snp.loc - x + 1

    #Get the snp (location) from the specific read
    count <- 0

    for(i in seq_along(nucl.num)){
      count <- count + nucl.num[i]
      if(count >= snp.start){break}
    }
    return(MDtag[i])
  }, a = fullMD, x = fullstart, MoreArgs = list(snp.loc = snp.loc), USE.NAMES = F)

  #Get read names from alt and ref reads
  names(alleles) <- names(alns)
  ref.reads <- sort(names(which(alleles != v)))
  alt.reads <- sort(names(which(alleles == v)))

  return(list(ref.reads = unique(ref.reads), alt.reads = unique(alt.reads)))
}


#' MDtag parser
#'
#' Takes a GenomicAlignments object containing the MDtag, and transforms it
#' into a vector of characters and numbers
#'
#' @param a Vector of MDtags (single characters)
#'
#' @return A named list of vectors, each vector a parsed version of MDtag:
#' - nucl.num: Numeric representation of MDtag.
#' - MDtag: a split version of MDtag
#'
#' @examples
#'
#'
#'
#' @export
getMD <- function(a){

  #extract matches and mismatches from MD tag
  numbers <- as.integer(str_extract_all(a, "[0-9]{1,}")[[1]])
  nucl <- str_extract_all(a, "[A-Z^]{1,}")[[1]]
  MDtag <- c(numbers, nucl)[order(c(seq_along(numbers)*2 - 1, seq_along(nucl)*2))]
  second <- str_extract(MDtag, "[A-Z]{2,}")
  
  # Small fix for large insertions (letters next to each other)
  if(any(!is.na(second))){
    rem <- which(!is.na(second))
    temp <- str_extract_all(MDtag[rem], "[A-Z]{1}")[[1]]
    MDtag <- MDtag[-rem]
    MDtag <- append(MDtag, temp, after = rem-1)
  } 

  if(length(MDtag) == 1){
    nucl.num <- as.integer(MDtag)
  } else{
    nucl.num <- numeric(length(numbers)+length(nucl))
    nucl.num[seq(2,length(nucl.num), 2)] <- 1
    if(sum(nucl.num != 1) == length(numbers)){
    nucl.num[nucl.num != 1] <- numbers } else{
      message("I don't understand this error")
    }
  }
  return(list(MDtag = MDtag, nucl.num = nucl.num))
}



