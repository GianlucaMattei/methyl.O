#' Plot a simple distribution of methylations.
#'
#' Plot distribution of beta values on chromosomes, this function has mainly an internal utility and is design for shinyApp version of package
#'
#' @param res the input table resulting from getMetAnnotations or scoreMets.
#' @param opt type of plot. Can be 'violin' or 'boxplot'
#' @param col color palette of plot. It must be one resulting from hcl.pals()
#'
#' @return plot
#'
#' @export

plotMethylationOverview <- function(res, opt, col) {
    plot.df <- data.frame(
        chr = res[[1]]$seqnames, Beta = res[[1]]$beta, direction = ifelse(res[[1]]$beta < 0, "down", "up"),
        stringsAsFactors = T
    )

    if (col == "BlackWhite") {
        palette <- rep(22, "black")
    } else {
        palette <- hcl.colors(length(unique(plot.df$chr)), palette = col)
    }

    if (opt == "boxplot") {
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
    } else if (opt == "violin") {
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
