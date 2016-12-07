#' Title
#'
#' @param score_matrix
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
order_pos_by_median <- function(score_matrix, ...) {

  # set median of each tuple as the genomic position
  pos <- rownames(score_matrix)
  chr <- limma::strsplit2(pos, ".", fixed=T)[,1]
  pos1 <- as.numeric(limma::strsplit2(pos, ".", fixed=T)[,2])
  pos2 <- as.numeric(limma::strsplit2(pos, ".", fixed=T)[,3])
  midpt <- floor((pos2 - pos1)/2)
  median_pos <- pos1 + midpt

  # sort the score matrix by median position (important for regionFinder and bumphunting)
  pos_df <- data.frame(chr=chr, pos=median_pos)
  o <- order(pos_df[,"chr"], pos_df[,"pos"])
  score_matrix[o,]

}
