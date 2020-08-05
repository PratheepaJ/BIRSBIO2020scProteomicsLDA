#' Aligned topic distribution across different chains
#'
#' @param theta Array. Stan output results of topic proportion.
#' @param sce SingleCellExperiment.
#' @param K Integer. Number of topics
#' @param iter Integer. Iteration used for each chain, including warm-up.
#' @param chain Integer. Number of chains used.
#' @param SampleID_name Character. The variable indicate sample ID.
#'
#' @return Matrix.
#' @export
#'

alignmentMatrix <- function(theta,
                            sce,
                            K,
                            iter = 2000,
                            chain = 4,
                            SampleID_name = "cell_id"){
  # theta is a 3 dimensional array (iterations * samples * topic)
  dimnames(theta)[[2]] = colnames(sce) %>% as.character()
  dimnames(theta)[[3]] = c(paste0("Topic_", seq(1,K)))

  # array to a dataframe
  theta_all = melt(theta)
  colnames(theta_all) = c("iteration", "Sample", "Topic", "topic.dis")
  theta_all$Chain = paste0("Chain ", rep(seq(1, chain), each = (iter/2)))
  theta_all$Sample = factor(theta_all$Sample)
  theta_all$Topic = factor(theta_all$Topic)

  # join the phyloseq sample data to theta_all
  sam = colData(sce) %>% data.frame()
  theta_all$Sample = as.character(theta_all$Sample)
  theta_all = left_join(theta_all, sam, by =c("Sample"= SampleID_name))
  theta_all$Chain = factor(theta_all$Chain)


  aligned <- matrix(nrow = K, ncol = chain)
  aligned[, 1] <- seq(1, K)
  corrTop <- numeric()
  for(j in 1:K){#Topic of first chain
    chains <- lapply(as.list(1:chain), function(x){
      for(top in 1:K){
        corrTop[top] <- cor(theta_all %>% dplyr::filter(Chain == "Chain 1") %>% dplyr::filter(Topic == paste0("Topic_", j)) %>% dplyr::select(topic.dis),
                            theta_all %>% dplyr::filter(Chain == paste0("Chain ",x)) %>% dplyr::filter(Topic == paste0("Topic_", top)) %>% dplyr::select(topic.dis))
      }
      return(which(corrTop == max(corrTop)))
    })

    aligned[j, 2:chain] <- unlist(chains)[2:chain]
  }

  return(aligned)
}
