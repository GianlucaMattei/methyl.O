#' Score the annotated methylation segments. 
#' 
#' Assigns a score to annotated methylated segments resulting from getMetAnnotations() function
#'  
#' @param results results from getMetAnnotations()
#' @param active.features the most effective features considering methylation. Choosen feature must be included in results. Additional feature names can be first exons (exons1) or first intron (intron1)
#' @param score.modifier value between 0-1. It specifies how the final score is computed by assigning different weights to methylations segments affecting genes expression or to genes already involved in pathologies. By increasing this value to 1, resulting scores will be focused on discovering segments affecting gene expression. A value equal to 0 will focus the results on genes involved in pathologies not considering the effect of methylation.
#'
#' @return data.frame of annotated table with assigned scores
#' 
#' @export

scoreAnnotatedDMRs <- function(results, active.features = c("promoters", "heads"), score.modifier = 0.5) {

  scale1 <- function(x) {
    (x - min(x)) / (max(x) - min(x))
  }


  if (any(active.features == "exons1")) {
    results$exons1 <- results$exons[results$exons$rank == 1, ]
  }

  if (any(active.features == "intron1")) {
    results$introns1 <- results$introns[results$introns$intron.rank == 1, ]
  }



  # candidates
  candidates <- results$genes$tag

  # features
  score.feature <- c()
  for (cand in candidates) {
    selected.features <- results[active.features]
    vals.list <- lapply(selected.features, function(x) {
      ind <- x$tag %in% cand
      if (any(ind)) {
        cur.val <- x[,grep('perc', colnames(x))][ind]
      }
    })
    score.feature <- c(score.feature, sum(unlist(vals.list)))
  }

  # beta and length
  score.betaLength <- abs(results$genes$beta)

  # CGIs and TF
  score.cpgis <- results$genes$CGIs + results$genes$TF

  # database score
  score.database <- results$genes$database.score

  # modifier
  score.modifier.beta <- 1 - score.modifier
  score.modifier.alpha <- score.modifier

  # collect scores
  results$genes$score <- ((scale1(score.database) * score.modifier.beta) + (scale1(score.cpgis + score.betaLength + score.feature) * score.modifier.alpha)) + (scale1(results$genes$genes.perc) * score.modifier.alpha)
  
  ind.ord <- order(results$genes$score, decreasing = T)
  results$genes <- results$genes[ind.ord, ]
  results$exons1 <- NULL
  return(results)
}
