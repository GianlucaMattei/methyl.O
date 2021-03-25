#' Plot EnrichR Results
#'
#' Visualize enrichment results and the contributes of hyper-metyhlated and hypo-methylated genes.
#'
#' @param enrichr.results data.frame resulting list from annotatedDMTs2Enrichr().
#' @param annotatedDMRs anotated DMRs list resultingfrom annotateDMRs() or scoreAnnotatedDMRs()
#' @param stat character, statistics to visualize. Accepted "P.value", "Adjusted.P.value" or "Overlap". Default = 'P.value'.
#' @param n numeric, value indicating the number of enrichment to plot, starting from the most enriched, to visualize. Default = 25.
#' @param plot.type  character, compute or not different linear models for upregulated and downregulated genes. Accepted: "simple" or "splitted". Default = "splitted".
#' @param col.hyper character, the color representing hyper-methylated genes. Default = 'Grey70'.
#' @param col.hypo character, the color representing hyper-methylated genes. Default = 'Grey30'.
#' @param pal.col character, the palette color for plot.type = 'lollipop'. Must be one from hcl.pals(). Default = 'Dynamic'.
#' @param thrs numeric vectors, statistical thresholds to plot. Default = c(0.01, 0.05)
#' @param thrs.cols characters vector, the colors to use for thresholds. Default = c('green', "yellow").
#'
#' @return plot of enrichR results
#'
#' @export

plotDMRs2Enrichr <- function(enrichr.results, annotatedDMRs, stat="P.value" , n = 25, plot.type = "barplot", pal.col = "Dynamic", col.hyper  = "#ff0000", col.hypo = "#00b3ff", thrs = c(0.01, 0.05), thrs.cols = c("green", "yellow")){

    if(!any((.packages()) %in% "enrichR")){
        library(enrichR)
    }

    if(nrow(enrichr.results)<n){
        n = nrow(enrichr.results)
    }

    if(plot.type == 'barplot'){
        stat.mat <- c()
        for (i in 1:n) {
            cur.genes <- unlist(strsplit(enrichr.results[i, ]$Genes, ";"))
            if(stat=="Overlap"){
                cur.stat <- as.numeric(unlist(strsplit(enrichr.results[i, stat], '/'))[1])
            } else { 
                cur.stat <- -log10(enrichr.results[i, stat])
            }
            cur.beta <- annotatedDMRs[[1]][which(annotatedDMRs[[1]]$symbol %in% cur.genes), "beta"]

            pos.tot <- sum(cur.beta > 0)
            if (pos.tot > 0) {
                perc <- length(cur.beta) / pos.tot
                cur.stat.pos <- cur.stat / perc
            } else {
                cur.stat.pos <- 0
            }

            neg.tot <- sum(cur.beta < 0)
            if (neg.tot > 0) {
                perc <- length(cur.beta) / neg.tot
                cur.stat.neg <- cur.stat / perc
            } else {
                cur.stat.neg <- 0
            }

            stat.mat <- cbind(c(cur.stat.pos, cur.stat.neg), stat.mat)
        }

        colnames(stat.mat) <- rev(enrichr.results$Term[1:n])
        stat.mat <- stat.mat[,order(colSums(stat.mat))]
        rownames(stat.mat) <- c("pos", "neg")
        par(mar = c(5, 15, 5, .1))
        barplot(tail(stat.mat, n), las = 1, horiz = T, col = c(col.hyper, col.hypo), cex.names = 0.85, main = "Enrichment Score")
        if(!is.null(thrs)){
            for(i in 1:length(thrs)){
                abline(v = -log10(thrs[i]), lty = 2, lwd = 1, col = thrs.cols[i])
            }
        }
        
    } else if (plot.type == 'lollipop'){
        enrichr.results.subset <- rev(-log10(enrichr.results[1:n, stat]))
        names(enrichr.results.subset) <- rev(enrichr.results$Term[1:n])
        enrichr.results.subset <- enrichr.results.subset[order(enrichr.results.subset)]
        par(mar = c(6, 15, 5, 0.1))
        plot(x = enrichr.results.subset, y = 1:length(enrichr.results.subset), pch = 21, bg = hcl.colors(n, pal.col), col = "black", main = "Enrichment Score", xlab = "", ylab = "", yaxt = "n")
        segments(x0 = 0, x1 = enrichr.results.subset, y0 = 1:length(enrichr.results.subset), y1 = 1:length(enrichr.results.subset), col = "black", lty = 2, lwd = .5)
        par(new=T)
        plot(x = enrichr.results.subset, y = 1:length(enrichr.results.subset), pch = 21, bg = hcl.colors(n, pal.col), col = "black", xlab = "", ylab = "", yaxt = "n", xaxt = "n", las = 2)
        ytick <- 1:length(enrichr.results.subset)
        axis(side = 2, at = ytick, names(enrichr.results.subset), las = 2, cex.axis = .85)

        if (!is.null(thrs)) {
            for (i in 1:length(thrs)) {
                abline(v = -log10(thrs[i]), lty = 2, lwd = 1, col = thrs.cols[i])
            }
        }

    }

}
