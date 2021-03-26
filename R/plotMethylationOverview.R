#' Plot a simple distribution of methylations.
#'
#' Plot distribution of beta values on chromosomes, this function has mainly an internal utility and is design for shinyApp version of package
#'
#' @param annotatedDMRs anotated DMRs list resultingfrom annotateDMRs() or scoreAnnotatedDMRs()
#' @param plot.type  character, compute or not different linear models for upregulated and downregulated genes. Accepted: "simple" or "splitted". Default = "splitted".
#' @param palette character, color palette of plot. It must be one resulting from hcl.pals()
#'
#' @return plot of enrichR results
#'
#' @export

plotMethylationOverview <- function(annotatedDMRs, plot.type, palette) {
    plot.df <- data.frame(
        chr = annotatedDMRs[[1]]$seqnames, Beta = annotatedDMRs[[1]]$beta, direction = ifelse(annotatedDMRs[[1]]$beta < 0, "down", "up"),
        stringsAsFactors = T
    )

    if (palette == "BlackWhite") {
        palette <- rep(22, "black")
    } else {
        palette <- hcl.colors(length(unique(plot.df$chr)), palette = palette)
    }

    if (plot.type == "boxplot") {
        ggplot2::ggplot(plot.df, ggplot2::aes(x = chr, y = Beta)) +
            ggplot2::geom_boxplot(position = "identity", fill=palette) +
            ggplot2::theme_bw() +
            ggplot2::theme(
                plot.background = ggplot2::element_blank(),
                panel.grid.major = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank(),
                axis.line.x = ggplot2::element_blank(),
                axis.ticks.x = ggplot2::element_blank(),
                axis.title.x = ggplot2::element_blank(),
                axis.text.x = ggplot2::element_text(angle = 90)
            ) +
            ggplot2::ylim(-1, 1) + ggplot2::ggtitle("Beta diff. values distribution by chromosome")
    } else if (plot.type == "violin") {
        ggplot2::ggplot(plot.df, ggplot2::aes(x = chr, y = Beta)) +
            ggplot2::geom_violin(ggplot2::aes(fill = chr)) +
            ggplot2::scale_fill_manual(values = palette) +
            ggplot2::theme_bw() +
            ggplot2::theme(
                legend.position = "none",
                plot.background = ggplot2::element_blank(),
                panel.grid.major = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank(),
                axis.line.x = ggplot2::element_blank(),
                axis.ticks.x = ggplot2::element_blank(),
                axis.title.x = ggplot2::element_blank(),
                axis.text.x = ggplot2::element_text(angle = 90)
            ) +
            ggplot2::ylim(-1, 1) + ggplot2::ggtitle("Beta diff. values distribution by chromosome")
    }
}
