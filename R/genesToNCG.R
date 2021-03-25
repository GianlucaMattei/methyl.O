#' Find genes annotated in NCG.
#'
#' Assess the presence of genes in results from annotateDMRs() in NCG database.
#'
#' @param annotatedDMRs anotated DMRs list resultingfrom annotateDMRs() or scoreAnnotatedDMRs()
#' @param ncg the NCG gene vectors
#' @param return.table logical, option TRUE return a table instead a vector of presences. Default = TRUE
#'
#' @return data.frame or vector of presences
#'
#' @export

genesToNCG <- function(annotatedDMRs, ncg, return.table = FALSE) {
    # check if symbol is present
    if (any(colnames(annotatedDMRs[[1]]) == "symbol")) {
        if (return.table == FALSE) {
            NCG <- lapply(annotatedDMRs, function(x) {
                apply(x, 1, function(y) {
                    as.numeric(any(ncg[, 1] %in% y["symbol"]))
                })
            })

            for (i in 1:length(NCG)) {
                annotatedDMRs[[i]]$NCG <- NCG[[i]]
            }
        } else {
            NCG <- lapply(annotatedDMRs, function(x) {
                apply(x, 1, function(y) {
                    paste0(ncg[which(ncg[, 1] %in% y["symbol"]), 2], collapse = "/")
                })
            })
            for (i in 1:length(NCG)) {
                annotatedDMRs[[i]]$NCG_type <- NCG[[i]]
            }
        }
    }
    return(annotatedDMRs)
}
