#' Divide read names by allele
#'
#' Takes a GenomicAlignments object and returns a list of read names dividied by
#' allele.
#'
#' @param alns GenomicAlignments object.
#' @param snp GRanges object containing SNP location.
#' @param v Nucleotide of reference (or alternative) allele.
#'
#' @return A named list of vectors, each vector containing read names for each
#'   allele.
#' @importFrom BiocGenerics start
#' @importFrom BiocGenerics end
#' @importFrom S4Vectors mcols
#' @importFrom stringr str_extract_all
#' @importFrom stringr str_extract
#' @importFrom stringr str_locate_all
#' @keywords internal
splitReads <- function(alns, v, snp) {
    
    fullMD <- mcols(alns)$MD
    fullstart <- start(alns)
    snp.loc <- start(snp)
    
    # Use MD tag from bam to check which reads have the ref or
    # alt allele
    alleles <- mapply(function(a, x, snp.loc) {
        
        # extract matches and mismatches from MD tag
        tag <- getMD(a)
        MDtag <- tag$MDtag
        nucl.num <- tag$nucl.num
        
        # Get right location of snp in sequence of read
        snp.start <- snp.loc - x + 1
        
        # Get the snp (location) from the specific read
        count <- 0
        
        for (i in seq_along(nucl.num)) {
            count <- count + nucl.num[i]
            if (count >= snp.start) {
                break
            }
        }
        return(MDtag[i])
    }, a = fullMD, x = fullstart, MoreArgs = list(snp.loc = snp.loc), 
        USE.NAMES = FALSE)
    
    # Get read names from alt and ref reads
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
#' @keywords internal
getMD <- function(a) {
    
    # extract matches and mismatches from MD tag
    numbers <- str_extract_all(a, "[0-9]{1,}")[[1]]
    numbers.loc <- str_locate_all(a, "[0-9]{1,}")[[1]][, 1]
    
    # nucl <- str_extract_all(a, '[A-Z^]{1,}')[[1]]
    nucl <- str_extract_all(a, "[A-Z^]{1}")[[1]]
    nucl.loc <- str_locate_all(a, "[A-Z^]{1}")[[1]][, 1]
    
    MDtag <- character(nchar(a))
    MDtag[numbers.loc] <- numbers
    MDtag[nucl.loc] <- nucl
    
    # remove empty and ^ positions
    MDtag <- MDtag[(MDtag != "" & MDtag != "^")]
    
    # extract MDtag as numbers
    if (length(MDtag) == 1) {
        nucl.num <- as.integer(MDtag)
    } else {
        nucl.num <- numeric(length(MDtag))
        nucl.num[MDtag %in% c("A", "C", "G", "T", "N")] <- 1
        
        if (sum(nucl.num != 1) == length(numbers)) {
            nucl.num[nucl.num != 1] <- as.integer(numbers)
        } else {
            message("I don't understand this error")
        }
    }
    return(list(MDtag = MDtag, nucl.num = nucl.num))
}
