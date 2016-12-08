
# B # Calculate ASM score for a list of samples in the output format of the result of read_tuples
# This functions uses the following other functions: process, calcScore, calcWeight

#' Title
#'
#' @param sample_list
#' @param beta
#' @param a
#' @param return_matrix
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
calc_asm <- function(sample_list, beta=0.5, a=0.2, return_matrix=TRUE, ...) {

  # Add 1 to every count and calculate odds ratio
  ASM_sample_list <- lapply(sample_list, calc_logodds)
  

  # calculate ASM score per tuple per sample
  for (i in 1:length(ASM_sample_list)) {
    df <- ASM_sample_list[[i]]
    df$ASM_score <- calcScore(beta,a,df)
    ASM_sample_list[[i]] <- df
  }

  # If return_matrix is TRUE, calculate a matrix of scores across samples, otherwise return list of samples
  if (return_matrix) {

    # get key of unique tuples
    all_pos <- do.call("rbind", ASM_sample_list)
    all_pos <- all_pos[,c("chr","pos1","pos2")]
    key <- paste0(all_pos$chr,'.',all_pos$pos1, '.', all_pos$pos2)
    key <- unique(key)

    # get matrix of ASM scores across all samples
    score_matrix <- matrix(data=NA, nrow=length(key), ncol=length(ASM_sample_list))
    colnames(score_matrix) <- names(ASM_sample_list)
    rownames(score_matrix) <- key

    for (i in 1:length(ASM_sample_list)) {
      df <- ASM_sample_list[[i]]
      key_s <- paste0(df$chr, '.', df$pos1, '.', df$pos2)
      m <- match(key, key_s)
      ind <- m[which(!is.na(m))]
      ind_m <- which(!is.na(m))
      score_matrix[ind_m,i] <- df$ASM_score[ind]
    }
    return(score_matrix)

  }
  else {
    return(ASM_sample_list)
  }

}

