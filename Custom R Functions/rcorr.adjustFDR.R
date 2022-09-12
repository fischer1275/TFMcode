rcorr.adjustFDR<-function (x, type = c("pearson", "spearman"), use = c("complete.obs", 
                                                      "pairwise.complete.obs")) 
{
  opt <- options(scipen = 5)
  on.exit(options(opt))
  type <- match.arg(type)
  use <- match.arg(use)
  x <- if (use == "complete.obs") 
    as.matrix(na.omit(x))
  else as.matrix(x)
  R <- rcorr(x, type = type)
  P <- P.unadj <- R$P
  p <- P[lower.tri(P)]
  adj.p <- p.adjust(p, method = "fdr")
  P[lower.tri(P)] <- adj.p
  P[upper.tri(P)] <- 0
  P <- P + t(P)
  #P <- ifelse(P < 1e-04, 0, P)
  #P <- format(round(P, 4))
  diag(P) <- NA
  #P[c(grep("0.0000", P), grep("^ 0$", P))] <- "<.0001"
  #P[grep("0.000$", P)] <- "<.001"
  #P.unadj <- ifelse(P.unadj < 1e-04, 0, P.unadj)
  #P.unadj <- format(round(P.unadj, 4))
  diag(P.unadj) <- NA
  #P.unadj[c(grep("0.0000$", P.unadj), grep("^ 0$", P.unadj))] <- "<.0001"
  #P.unadj[grep("0.000$", P.unadj)] <- "<.001"
  result <- list(R = R, P = P, P.unadj = P.unadj, type = type)
  class(result) <- "rcorr.adjustFDR"
  result
}
