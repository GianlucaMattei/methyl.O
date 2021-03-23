#' Plot the beta value of TF and related target genes expression
#'
#' Barplot of the beta value of TF and the expression of the target genes
#'
#' @param associated.table data.frame. Corresponding to resulting data.frame from tfExprsCor().
#' @param symbol feature to correlate. Must be from names of resulting list from gemMetAnnotations.  Additional feature names can be first exons (exons1) or first intron (intron1). To use more than one feature use c().
#' @param col.meth color for beta value
#' @param pals.bars palette for gene expresion
#' 
#' @return plot showing methylation levels of TF and expresion of target genes
#'
#' @export


plotTFs2Exprs <- function(associated.table, symbol, col.meth = '#8e0000', pals.bars = 'Cold'){
    cur.slice = associated.table[associated.table$symbol %in% symbol,]
    cur.tf <- cur.slice[1,1]
    cur.bars <- cur.slice['target.expression']
    cur.bars <- rbind(cur.tf,0,cur.bars)
    rownames(cur.bars) <- c(cur.slice[1,2],'',cur.slice[,'target.symbol'])
    marg <- max(abs(c(min(cur.bars) , max(cur.bars))))
    ylimit <- round(marg)
    colors <- c(col.meth, hcl.colors(nrow(cur.bars)-1,pals.bars))
    # colors <- hcl.colors(nrow(cur.bars),pals.bars)
    bords <- rep('black', length(colors))
    bords[2] = 'white'
    barplot(height=cur.bars[,1],names=rownames(cur.bars), las=2, ylim = c(ylimit*-1,ylimit), col = colors, axes=FALSE, border = bords, main = paste0(cur.slice[1,2], ' beta val. and target genes expression'))
    axis(4,at=seq(ylimit*-1,ylimit))
    axis(2,at=seq(-1,1))
}
