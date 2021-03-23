#' Find genes annotated in NCG.
#'
#' Verify the presence in NCG of genes in results from getMetAnnotations()
#'
#' @param results results from getMetAnnotations()
#' @param ncg the NCG genes vector
#' @param return.table optition to return a table instead a vector of presence
#'
#' @return data.frame
#'
#' @export

genesToNCG <- function(results, ncg, return.table = FALSE) {
    # check if symbol is present
    if (any(colnames(results[[1]]) == "symbol")) {
        if (return.table == FALSE) {
            NCG <- lapply(results, function(x) {
                apply(x, 1, function(y) {
                    as.numeric(any(ncg[, 1] %in% y["symbol"]))
                })
            })

            for (i in 1:length(NCG)) {
                results[[i]]$NCG <- NCG[[i]]
            }
        } else {
            NCG <- lapply(results, function(x) {
                apply(x, 1, function(y) {
                    paste0(ncg[which(ncg[, 1] %in% y["symbol"]), 2], collapse = "/")
                })
            })
            for (i in 1:length(NCG)) {
                results[[i]]$NCG_type <- NCG[[i]]
            }
        }
    }
    return(results)
}
