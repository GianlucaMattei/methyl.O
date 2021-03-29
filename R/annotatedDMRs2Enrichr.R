#' Compute correlation between expression and methylation levels. 
#'
#' @param annotatedDMRs anotated DMRs list resultingfrom annotateDMRs() or scoreAnnotatedDMRs()
#' @param expressionProfile data.frame, expression profile. 
#' @param active.features character vectors, containing features to correlate. Must be from names of resulting list from annotateDMRs. Additional feature names can be first exons (exons1) or first intron (intron1). To use more than one feature use c(). Default = c("promoters", "heads")
#' @param col.genes numeric, the column of expressionProfile data.frame with gene Ids. If NULL geneIDs will be taken from rownames() of expressionProfile. Default = 0.
#' @param col.stat numeric, the column of expressionProfile data.frame with the statistics to use. Default = 6.
#' @param stat.thr numeric, threshold for statistical significance. Default = 0.05
#' @param col.logFC numeric, the column of expressionProfile data.frame with log. fold change. Default = 2
#' @param logfc.thr numeric, threshold value for log. fold change. Default = 0.
#' @param convert.genes boolean, used to indicate if gene ids have to be translated in official gene symbols. Default = FALSE
#' @param convert.from character, used annotation for gene in expressionProfile to be converted to symbols gene IDs. Accepted: c("ENTREZID" ,"EXONID" ,"GENEBIOTYPE" ,"GENEID" ,"GENENAME" ,"PROTDOMID" ,"PROTEINDOMAINID" ,"PROTEINDOMAINSOURCE" ,"PROTEINID" ,""SEQSTRAND" ,"SYMBOL" ,"TXBIOTYPE" ,"TXID" ,"TXNAME" ,"UNIPROTID")
#' @param beta.thr numeric, beta difference threshold value. Default = 0.3.
#' @param param.type character, threshold parameter to filter methylations overlapping the selected features. Accetped c("dmr.length", "overlap.length", "overlap.percentage"). Default = "overlap.length".
#' @param overlap.param.thr nuemric, threshold value for selected parameter to filter methylations overlapping the selected features. Default = 100
#' @param line.col character, color of lines at x=0, y=0. Default = "lightgray"
#' @param lmfit.col1 character, color of linear model line 1 or for simple plot. Default = "red"
#' @param lmfit.col2 character, color of linear model line 2. Default = "green"
#' @param pal character, color palette. hcl.pals() to show available. Default = "RdGy"
#' @param plot.type  character, compute or not different linear models for upregulated and downregulated genes. Accepted: "simple" or "splitted". Default = "splitted".
#' @param show.text logical, indicating if print gene names in the final plot. Accepted values: TRUE or FALSE. Default = FALSE.
#' @param filter.by.genes character vectors, gene symbols used for filtering output
#' @param cor.type character, correlation method. Available "pearson", "kendall" or "spearman" correlation, Default = "pearson"
#' @param return.table logical, TRUE return a data.frame instead a plot. Default = FALSE
#' 
#' @return plot of correlation or (if return.table = TRUE) a data.frame of genes expression associated to beta methylation values.
#'
#' @export

annotatedDMRs2Exprs <- function(annotatedDMRs, expressionProfile, active.features = c("promoters", "heads"), col.genes = 0, col.stat = 6, stat.thr = 0.05, col.logFC = 2, logfc.thr = 0, convert.genes = FALSE, convert.from, beta.thr = .3, overlap.param.thr = 100, param.type = 'overlap.length', line.col = 'lightgrey', lmfit.col1='red', lmfit.col2='green', pal = 'RdGy', plot.type = 'splitted', show.text = FALSE, cor.type = 'pearson', filter.by.genes = NULL, return.table = FALSE){

    outTable = data.frame()

    if (any(active.features == "exons1")) {
        annotatedDMRs$exons1 <- annotatedDMRs$exons[annotatedDMRs$exons$rank == 1, ]
    }

    if (any(active.features == "introns1")) {
        annotatedDMRs$introns1 <- annotatedDMRs$introns[annotatedDMRs$introns$intron.rank == 1, ]
    }

    if(param.type != 'overlap.length' ){
        warning('Consider to change value threshold if you are changing param.type')
    }

    # create objects to work on
    # df with features
    toAssociate <- c()

    if(param.type == 'overlap.percentage'){

        for(f in active.features){
            ind <- match(f, names(annotatedDMRs))
            ind.perc <- grep('perc', colnames(annotatedDMRs[[ind]]))
            curFeat <- annotatedDMRs[[f]][, c("beta", "symbol", colnames(annotatedDMRs[[f]])[ind.perc])]
            curFeat$feature <- rep(names(annotatedDMRs)[ind],nrow(curFeat))
            colnames(curFeat)=c('beta', 'symbol', 'perc.overlap', 'feature')
            toAssociate <- rbind(toAssociate, curFeat)
        }

        toAssociate <- toAssociate[toAssociate$perc.overlap >= overlap.param.thr, ]
    } else if (param.type == "overlap.length"){

        for(f in active.features){
            ind <- match(f, names(annotatedDMRs))
            ind.owidth <- grep('overlap\\.width', colnames(annotatedDMRs[[ind]]))
            curFeat <- annotatedDMRs[[f]][, c("beta", "symbol", colnames(annotatedDMRs[[f]])[ind.owidth])]
            curFeat$feature <- rep(names(annotatedDMRs)[ind],nrow(curFeat))
            colnames(curFeat)=c('beta', 'symbol', 'width.overlap', 'feature')
            toAssociate <- rbind(toAssociate, curFeat)
        }

        toAssociate <- toAssociate[toAssociate$width.overlap >= overlap.param.thr, ]

    } else if (param.type == "dmr.length"){
        
        for(f in active.features){
            ind <- match(f, names(annotatedDMRs))
            ind.width <- 4
            curFeat <- annotatedDMRs[[f]][, c("beta", "symbol", colnames(annotatedDMRs[[f]])[ind.width])]
            curFeat$feature <- rep(names(annotatedDMRs)[ind],nrow(curFeat))
            colnames(curFeat)=c('beta', 'symbol', 'width.dmr', 'feature')
            toAssociate <- rbind(toAssociate, curFeat)
        }

        toAssociate <- toAssociate[toAssociate$width.dmr >= overlap.param.thr, ]

    } else {
        stop("Impossible to find parameter for filtering")
    }
    
    if(nrow(toAssociate)==0){
        stop('None of DMRs can pass overlap parameter threshold')
    }

    toAssociate <- toAssociate[abs(toAssociate$beta)>= beta.thr ,]

    # gene expression
    if(!is.null(col.genes) & col.genes!=0){
        expression.genes <- as.character(expressionProfile[,col.genes])
    } else if (missing(col.genes) | col.genes==0) {
        expression.genes <- rownames(expressionProfile)
    }

    if(all(expression.genes[1:6] == c("1","2","3","4","5","6"))){
        stop("Please select the genes column \n retrieved genes are not symbols")
    }

    expression.stat <- expressionProfile[, col.stat]
    expression.logFC <- expressionProfile[, col.logFC]
    expression.parsed <- data.frame(stringsAsFactors = F, genes = expression.genes, logFC = expression.logFC, expression.stat)
    expression.parsed <- expression.parsed[!is.na(expression.parsed[,1]),]
    expression.parsed <- expression.parsed[expression.parsed$expression.stat < stat.thr,]

    if(convert.genes){
            to.translate <- expression.parsed$genes
            annotations <- biomaRt::select(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75, keys = as.character(to.translate), columns = "SYMBOL", keytype = convert.from)
            translated <- annotations[,2][match(expression.parsed$genes, annotations[, match(convert.from, colnames(annotations))])]
            expression.parsed$genes <- translated
            expression.parsed$sourceID <- to.translate
    }



    # start association
    outTable <- c()
    for (i in 1:nrow(toAssociate)) {
        gene <- toAssociate[i, 'symbol']
        indB <- match(gene, expression.parsed$genes)
        if(!is.na(indB)){
            cur.fc <- expression.parsed[indB, 2]
            if (abs(cur.fc) >= logfc.thr) {
                curSlice <- cbind(toAssociate[i, ], fc = cur.fc)
                outTable <- rbind(outTable, curSlice)
            }
        }
    }

    if(!is.null(filter.by.genes)){
        outTable.new <- outTable[which(outTable$symbol %in% filter.by.genes),]
        if(nrow(outTable.new)<2){
            warning('none of the gene fol filtering is in the final table. Returning non filtered table')
        } else {
            outTable <- outTable.new
        }
    }



    if(return.table==FALSE){
        if(plot.type == 'simple'){
            score <- scale((outTable[, 1] * outTable[, 5]) * -1)
            palette <- hcl.colors(length(outTable[, 2]), pal)
            ind.order <- order(score, decreasing = T)
            score <- score[ind.order]
            outTable <- outTable[ind.order, ]

            par(mar=c(4,4,4,4))
            layout(cbind(matrix(1, ncol = 5, nrow = 6), matrix(c(3,3,2,2,3,3), ncol = 1, nrow = 6)))
            plot(x = outTable[, 1], y = outTable[, 5], xlab = "beta", ylab = "logFC", xlim = c(-1, 1), ylim = c(-max(abs(outTable$fc)), max(abs(outTable$fc))), col = palette, pch=19, main = 'beta - logFC correlation', cex = 1.3)
            if(show.text){
                text(x = outTable[, 1], y = outTable[, 5], labels = outTable[, 2])
            }
            abline(v = 0, h = 0, col = line.col, lty = 2)
            abline(lm(as.numeric(outTable[, 5]) ~ as.numeric(outTable[, 1]) - 1), col = lmfit.col1, lty = 2, lwd = 2)

            correlation <- cor.test(as.numeric(outTable[, 5]), as.numeric(outTable[, 1]), method = cor.type)
            mtext(paste("R:", round(correlation$statistic, 3), " pval:", round(correlation$p.val, 3), sep = " "), side = 1, at = -0.8, line = 2.5, cex = 0.7)

            par(mar=c(2,2,2,2))
            gradient <- as.raster(matrix(palette, ncol = 1))
            plot(c(0, 2), c(min(score), max(score)), type = "n", axes = F, xlab = "", ylab = "", main = "Score")
            steps <- ((abs(max(score)) + abs(min(score))) / 4)
            text(x = 1.5, y = seq(min(score), max(score), steps), labels = round(seq(min(score), max(score), steps), 3))
            rasterImage(gradient, 0, min(score), 1, max(score))

        } else if(plot.type == 'splitted'){
            score <- scale((outTable[, 1] * outTable[, 5]) * -1)
            palette <- hcl.colors(length(outTable[, 2]), pal)
            ind.order <- order(score, decreasing = T)
            score <- score[ind.order]
            outTable <- outTable[ind.order, ]

            indPos <- outTable[, 1] > 0
            indNeg <- outTable[, 1] < 0

            par(mar = c(4, 4, 4, 4))
            layout(cbind(matrix(1, ncol = 5, nrow = 6), matrix(c(3, 3, 2, 2, 4, 4), ncol = 1, nrow = 6)))

            if (sum(indPos) > 2) {
                plot(x = outTable[indPos, 1], y = outTable[indPos, 5], xlab = "beta", ylab = "logFC", xlim = c(-1, 1), ylim = c(-max(abs(outTable$fc)), max(abs(outTable$fc))), pch = 19, cex = 1.3, col = palette[indPos])
                if(show.text){
                    text(x = outTable[indPos, 1], y = outTable[indPos, 5], labels = outTable[indPos, 2])
                }
                abline(lm(as.numeric(outTable[indPos, 5]) ~ as.numeric(outTable[indPos, 1]) - 1), col = lmfit.col1, lty = 2, lwd = 2)
                par(new = T)
                correlationPos <- cor.test(as.numeric(outTable[indPos, 5]), as.numeric(outTable[indPos, 1]), method = cor.type)
                mtext(paste("positive R:", round(correlationPos$statistic, 3), " pval:", round(correlationPos$p.val, 3), sep = " "), side = 1, at = -0.8, line = 2.5, cex = 0.7)
            } else {
                correlationPos = ''
            }
            if (sum(indNeg) > 2) {
                if(sum(indPos) > 2){
                    par(new=T)
                }
                plot(x = outTable[indNeg, 1], y = outTable[indNeg, 5], xlab = "beta", ylab = "logFC", xlim = c(-1, 1), ylim = c(-max(abs(outTable$fc)), max(abs(outTable$fc))), pch = 19, cex = 1.3, col = palette[indNeg])
                if(show.text){
                    text(x = outTable[indNeg, 1], y = outTable[indNeg, 5], labels = outTable[indNeg, 2])
                }
                abline(lm(as.numeric(outTable[indNeg, 5]) ~ as.numeric(outTable[indNeg, 1]) - 1), col = lmfit.col2, lty = 2, lwd = 2)
                correlationNeg <- cor.test(as.numeric(outTable[indNeg, 5]), as.numeric(outTable[indNeg, 1]), method = cor.type)
                mtext(paste("negative R:", round(correlationNeg$statistic, 3), " pval:", round(correlationNeg$p.val, 3), sep = " "), side = 1, at = 0.8, line = 2.5, cex = 0.7)
            } else {
                correlationNeg = ''
            }
            abline(v = 0, h = 0, col = line.col, lty = 2)
            

            
            par(mar = c(2, 2, 2, 2))
            gradient <- as.raster(matrix(palette, ncol = 1))
            plot(c(0, 2), c(min(score), max(score)), type = "n", axes = F, xlab = "", ylab = "", main = "Score")
            steps <- ((abs(max(score)) + abs(min(score))) / 4)
            text(x = 1.5, y = seq(min(score), max(score), steps), labels = round(seq(min(score), max(score), steps), 3))
            rasterImage(gradient, 0, min(score), 1, max(score))

            plot.new()
            par(mar = c(0, 0, 0, 0), xpd = T)
            legend("center", c("Beta Positive Corr.", "Beta Negative Corr."), pch = 15, col = c(lmfit.col1, lmfit.col2), box.lwd = 0 )

        }

    } else if(return.table==TRUE){
        score <- scale((outTable[, 1] * outTable[, 5]) * -1)
        ind.order <- order(score, decreasing = T)
        score <- score[ind.order]
        outTable <- outTable[ind.order, ]
        outTable$score <- score
        return(outTable)
    }
}
