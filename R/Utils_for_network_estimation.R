stein_loss <- function(m1, m2){
  n <- nrow(m1)
  P <- chol2inv(chol(m2)) %*% m1
  return(sum(diag(P)) - determinant(P)$modulus[1] - n)
}

conf_matrix <- function(truth, estimation, margins = FALSE, normalize = FALSE,
                        undirected = TRUE) {
  same_edges <- truth * estimation
  diff <- truth - estimation
  summ <- truth + estimation
  p <- dim(truth)[1]
  max_edges <- (p^2 - p)
  tp <- sum(same_edges)
  tn <- sum(summ == 0) - p
  fp <- sum(diff == -1)
  fn <- sum((same_edges - truth) == -1)
  P <- sum(truth)
  N <- max_edges - P
  EP <- sum(estimation)
  EN <- max_edges - EP
  cm <- matrix(c(tp, fn, fp, tn), nrow = 2, byrow = TRUE,
               dimnames = list(c("True P", "True N"), c("Estim. P", "Estim. N")))
  if (margins) {
    cm <- matrix(c(tp, fn, P, fp, tn, N, EP, EN, max_edges), nrow = 3, byrow = TRUE,
                 dimnames = list(c("True P", "True N", "Sum"), c("Estim. P", "Estim. N", "Sum")))
  }
  if (undirected) {
    cm <- cm * 0.5
  }
  if (normalize) {
    cm <- matrix(c(tp/P, fn/P, fp/N, tn/N), nrow = 2, byrow = TRUE,
                 dimnames = list(c("True P", "True N"), c("Estim. P", "Estim. N")))
  }
  return(cm)
}

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
  fm <- sqrt(ppv * tpr)
  mk <- ppv + npv - 1
  acc <- (tp + tn) / (tp + tn + fn + fp)
  bal_acc <- (tpr + tnr) / 2
  F1_score <- 2 * (ppv * tpr) / (ppv + tpr)
  mcc <- (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  results <- data.frame(ACC = acc, ACC_bal = bal_acc, MCC = mcc, F1 = F1_score,
                        TPR = tpr, TNR = tnr, PPV = ppv, NPV = npv, FPR = fpr,
                        FNR = fnr, FDR = fdr, FOR = FOR, PT = pt, TS = ts, FM = fm, MK = mk,
                        LRp = lr_plus, LRn = lr_neg)
  return(results)
}