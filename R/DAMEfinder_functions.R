## DAMEfinder functions


# A # Read in list of methtuple files and filter out according to set min_cov and max_gap. 

# input: list of methtuple files, names of samples (optional)
# output: list of data frames of samples containing tuple information remaining after filtering, 
# ordered by chr, pos1, pos2 


read_tuples <- function(files, sample_names=seq.int(1:length(files)), min_coverage=10, max_gap=150, ... ) {
  # read in methtuple files
  methtuple_list <- parallel::mclapply(files, function(f){read.table(f, header=T)}, 
                                       mc.cores=min(length(files), ceiling(parallel::detectCores()/3)))
  
  # order each data frame by chr and pos1, add coverage column, and add tuple distance column
  methtuple_list_o <- lapply(methtuple_list, function(df) {
    df$cov <- with(df, MM + UU + UM + MU)
    df$inter_dist <- with(df, pos2-pos1)
    df[order(df[,"chr"], df[,"pos1"], df[,"pos2"]),]
  }) 
  
  # filter out low coverage and maxGap
  methtuple_list_o_f <- lapply(methtuple_list_o, function(df) {
    w <- which(df$cov >= min_coverage & df$inter_dist <= max_gap)
    return(df[w,])
  })
  
  return(methtuple_list_o_f)

}
# END OF A ############

# B # Calculate ASM score for a list of samples in the output format of the result of read_tuples
# This functions uses the following other functions: process, calcScore, calcWeight

calcASM <- function(sample_list, beta=0.5, a=0.2, return_matrix=TRUE, ...) {
  # Add 1 to every count and calculate odds ratio
  ASM_sample_list <- parallel::mclapply(sample_list, process, 
                                        mc.cores = min(length(sample_list), ceiling(parallel::detectCores()/6)))
  
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

# calculate the weight per site given beta and a
calcWeight <- function(beta, MM, UU, a) {
  
  # set quantile range
  theta = seq(0.005, 0.995, length = 1000)
  
  # calculate pbeta (note that this is the cumulative prob distr). To get the densities, use dbeta.
  posterior <- pbeta(theta, beta+MM, beta+UU, lower.tail = T)
  
  # get prob(0.3<b<0.7, ie a=0.2): our new weight
  w_1 <- which(theta<(0.5-a))
  index_1 <- w_1[length(w_1)] + 1
  
  w_2 <- which(theta<(0.5+a))
  index_2 <- w_2[length(w_2)]
  
  posterior[index_2] - posterior[index_1]
  
}


# get a vector of the posterior scores for all positions given beta and a
calcScore <- function(beta, a, df) {
  
  weights <- parallel::mcmapply(calcWeight, beta, df$MM, df$UU, a, mc.cores = min(7, ceiling(parallel::detectCores()/6)))
  df$log_odds_ratio*weights
}

process <- function(s) {
  s$strand <- "*"
  s$MM <- s$MM + 1
  s$MU <- s$MU + 1
  s$UM <- s$UM + 1
  s$UU <- s$UU + 1
  s$cov <- s$cov + 4
  # calc log odds ratio
  ratio <- (s$MM*s$UU)/(s$MU*s$UM)
  s$log_odds_ratio <-  log(ratio, base=10)
  return(s)
}


# END OF B ############

# C # Order matrix of scores by tuple median. 
# The rownames must have the following format: "chr.pos1.pos2"

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

# END OF C ############



# D # Transforms scores with square root function that preserves the signs
# Uses the following function: modulus_sqrt

transform_scores <- function(score_matrix, ...) {
  sm_t <- matrix(data=NA, nrow=nrow(score_matrix), ncol=ncol(score_matrix))
  colnames(sm_t) <- colnames(score_matrix)
  rownames(sm_t) <- rownames(score_matrix)
  for (i in 1:ncol(sm_t)) {
    sm_t[,i] <- modulus_sqrt(score_matrix[,i])
  }
  return(sm_t)
}

modulus_sqrt <- function(values) {
  t_values <- abs(values)
  sign(values)*sqrt(t_values)
}


# END OF D ############

# E # Get t-stats

get_t_stats <- function(sm_t, design, ...) {
  

  # moderated t-statistic
  coeff <- 2 # the column in the design matrix to consider
  
  fit <- limma::lmFit(sm_t, mod, method = "robust")
  
  fit2 <- limma::eBayes(fit, proportion = 0.01)
  fit2$t[, coeff]
}

# END OF E ############

# F # Smooth t-statistics
get_smoothed_t_stats <- function(betas, ...) {
  
  rows <- names(betas)
  chr <- limma::strsplit2(rows, ".", fixed=T)[,1]
  pos1 <- as.numeric(limma::strsplit2(rows, ".", fixed=T)[,2])
  pos2 <- as.numeric(limma::strsplit2(rows, ".", fixed=T)[,3])
  midpt <- floor((pos2 - pos1)/2)
  pos <- pos1 + midpt
  
  pns <- bumphunter::clusterMaker(chr, pos, maxGap = 300)
  
  verbose <- TRUE
  Q <- 0.9
  
  smooth <- bumphunter::smoother(y = betas, x = pos, cluster = pns, smoothFunction = loessByCluster, 
                     verbose = verbose)
  smooth$fitted
  
}

# END OF F ############

# G # Detect DAMEs
find_dames <- function(nm, smoothed_betas, verbose=TRUE, ...) {

  chr <- limma::strsplit2(nm, ".", fixed=T)[,1]
  pos1 <- as.numeric(limma::strsplit2(nm, ".", fixed=T)[,2])
  pos2 <- as.numeric(limma::strsplit2(nm, ".", fixed=T)[,3])
  midpt <- floor((pos2 - pos1)/2)
  pos <- pos1 + midpt
  
  pns <- bumphunter::clusterMaker(chr, pos, maxGap = 300)
  
  Q <- 0.9
  # Detect DAMEs
  K <- quantile(abs(smoothed_betas),Q, na.rm=T)
  dames <- bumphunter::regionFinder(x = smoothed_betas, chr = chr, pos = pos, cluster = pns, 
                        cutoff = K, verbose = verbose)
  
  
}

# END OF G ############







