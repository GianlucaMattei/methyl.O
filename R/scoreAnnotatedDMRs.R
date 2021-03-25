#' Score the annotated methylation segments. 
#' 
#' Assigns a score to annotated methylated segments resulting from annotateDMRs() function
#'  
#' @param annotatedDMRs anotated DMRs list resultingfrom annotateDMRs()
#' @param active.features character vectors, containing features to correlate. Must be from names of resulting list from annotateDMRs. Additional feature names can be first exons (exons1) or first intron (intron1). To use more than one feature use c(). Default = c("promoters", "heads")
#' @param score.modifier numeric, value between 0-1. It specifies how the final score is computed by assigning different weights to the methylation charactersistics of enhacners or to genes already involved in pathologies. By increasing this value to 1, resulting scores will be focused on discovering segments affecting gene expression. A value equal to 0 will focus the results on enahcners involving genes associated to pathologies, not considering the effect of methylation. Default = 0.5
#'
#' @return data.frame of annotated DMRs with assigned scores
#' 
#' @export

scoreAnnotatedDMRs <- function(annotatedDMRs, active.features = c("promoters", "heads"), score.modifier = 0.5) {

  scale1 <- function(x) {
    (x - min(x)) / (max(x) - min(x))
  }


  if (any(active.features == "exons1")) {
    annotatedDMRs$exons1 <- annotatedDMRs$exons[annotatedDMRs$exons$rank == 1, ]
  }

  if (any(active.features == "intron1")) {
    annotatedDMRs$introns1 <- annotatedDMRs$introns[annotatedDMRs$introns$intron.rank == 1, ]
  }



  # candidates
  candidates <- annotatedDMRs$genes$tag

  # features
  score.feature <- c()
  for (cand in candidates) {
    selected.features <- annotatedDMRs[active.features]
    vals.list <- lapply(selected.features, function(x) {
      ind <- x$tag %in% cand
      if (any(ind)) {
        cur.val <- x[,grep('perc', colnames(x))][ind]
      }
    })
    score.feature <- c(score.feature, sum(unlist(vals.list)))
  }

  # beta and length
  score.betaLength <- abs(annotatedDMRs$genes$beta)

  # CGIs and TF
  score.cpgis <- annotatedDMRs$genes$CGIs + annotatedDMRs$genes$TF

  # database score
  score.database <- annotatedDMRs$genes$database.score

  # modifier
  score.modifier.beta <- 1 - score.modifier
  score.modifier.alpha <- score.modifier

  # collect scores
  annotatedDMRs$genes$score <- ((scale1(score.database) * score.modifier.beta) + (scale1(score.cpgis + score.betaLength + score.feature) * score.modifier.alpha)) + (scale1(annotatedDMRs$genes$genes.perc) * score.modifier.alpha)
  
  ind.ord <- order(annotatedDMRs$genes$score, decreasing = T)
  annotatedDMRs$genes <- annotatedDMRs$genes[ind.ord, ]
  annotatedDMRs$exons1 <- NULL
  return(annotatedDMRs)
}
