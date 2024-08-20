# This file contains miscellaneous functions


unpack_score_list <- function(calculated_scores){
  ACC <- calculated_scores$ACC
  ACC_bal <- calculated_scores$ACC_bal
  MCC <- calculated_scores$MCC
  F1  <- calculated_scores$F1
  TPR <- calculated_scores$TPR
  TNR <- calculated_scores$TNR
  PPV <- calculated_scores$PPV
  NPV <- calculated_scores$NPV
  FPR <- calculated_scores$FPR
  FNR <- calculated_scores$FNR
  FDR <- calculated_scores$FDR
  FOR <- calculated_scores$FOR
  LRp <- calculated_scores$LRp
  LRn <- calculated_scores$LRn


  return(matrix(c(ACC,ACC_bal,MCC,F1,TPR,TNR,PPV,NPV,FPR,FNR,FDR,FOR,LRp,LRn), nrow = 1, ncol = 14))

}



#' walktrap glustering for know network
cluster_tags <- function(solution){
  g <- graph_from_adjacency_matrix(adjmatrix = solution, mode = "undirected")
  c <- cluster_walktrap(g, weights = NULL, steps = 4, merges = TRUE)
  return(as.matrix(membership(c)))
}


# from network_visualisation_and_scores.R
# author: Tuomas Hautamäki
adjacency_matrix <- function(coef_matrix, condition = "AND") {
  p <- dim(coef_matrix)[1]
  adj_matrix <- matrix(0, nrow = p, ncol = p)
  if (condition != "AND" && condition != "OR") {
    print("Wrong condition! Only AND and OR are accepted.")
    return(adj_matrix)
  }
  for (i in 1:p) {
    for (j in i:p) {
      if (i == j) {
        next
      }
      if (condition == "AND") {
        if (coef_matrix[i, j] != 0 && coef_matrix[j, i] != 0) {
          adj_matrix[i, j] <- 1
          adj_matrix[j, i] <- 1
        }
      }
      else if (condition == "OR") {
        if (coef_matrix[i, j] != 0 || coef_matrix[j, i] != 0) {
          adj_matrix[i, j] <- 1
          adj_matrix[j, i] <- 1
        }
      }
    }
  }
  return(adj_matrix)
}

# from network_visualisation_and_scores.R
# author: Tuomas Hautamäki
calculate_scores <- function(cm) {
  tp <- cm[1,1]
  tn <- cm[2,2]
  fp <- cm[2,1]
  fn <- cm[1,2]
  tpr <- tp / (tp + fn)
  tnr <- tn / (tn + fp)
  ppv <- tp / (tp + fp)
  npv <- tn / (tn + fn)
  fnr <- 1 - tpr
  fpr <- 1 - tnr
  fdr <- 1 - ppv
  FOR <- 1 - npv
  lr_plus <- tpr / fpr
  lr_neg <- fnr / tnr
  pt <- sqrt(fpr) / (sqrt(tpr) + sqrt(fpr))
  ts <- tp / (tp + fn + fp)
  acc <- (tp + tn) / (tp + tn + fn + fp)
  bal_acc <- (tpr + tnr) / 2
  F1_score <- 2 * (ppv * tpr) / (ppv + tpr)
  mcc <- (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  results <- data.frame(ACC = acc, ACC_bal = bal_acc, MCC = mcc, F1 = F1_score,
                        TPR = tpr, TNR = tnr, PPV = ppv, NPV = npv, FPR = fpr,
                        FNR = fnr, FDR = fdr, FOR = FOR, LRp = lr_plus, LRn = lr_neg)
  return(results)
}
