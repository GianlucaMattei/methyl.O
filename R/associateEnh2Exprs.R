#' Compute expression methylation correlation 
#'
#' Compute and plot correlation between expression and methylation
#'
#' @param annotatedEnhancers Data.frame. It corresponds to resulting list from enhancerAnnotations()
#' @param hg hg19, hg38. Version of the enhancer database
#' @param expression expression data.frame
#' @param enhancer.db character, which database to use between 'FANTOM5' or '4DGenome'. Default = "FANTOM5"
#' @param col.genes number of the column of expression data.frame with gene Ids Default = 0
#' @param col.stat number of the column of expression data.frame with the selected statistics, Default = 6
#' @param stat.thr threshold for statistical significance. Default = 0.05
#' @param col.logFC number of the column of expression data.frame with lofFC. Default = 2
#' @param logfc.thr logFC threshold value. Default = 0.5
#' @param convert.genes boolean value used to indicate if gene ids have to be translated in official gene symbols
#' @param convert.from gene id annotation of expression data.frame to convert to symbols gene ids. Character: ENTREZID TXNAME TXID GENEID
#' @param beta.thr beta difference threshold value
#' @param param.type threshold parameter to filter methylations overlapping enhancers. Accetped c("dmr.length", "overlap.length", "overlap.percentage"). Default = "overlap.percentage".
#' @param show.text Logical indicating if print gene names in the final plot. Accepted values: TRUE or FALSE. Default = FALSE
#' @param overlap.param.thr threshold value for selected parameter to filter methylations overlapping enhancers. Default = 40
#' @param line.col color of lines at x=0, y=0. Defalt = lightgray
#' @param lmfit.col1 color of linear model line 1 or for simple plot. Defalt = red
#' @param lmfit.col2 color of linear model line 2. Defalt = green
#' @param pal Color palette for plot. Default = RdGy 
#' @param plot.type simple or splitted. Compute or not different linear models for hyper/hypo genes. Default = splitted
#' @param cor.type pearson, kendall or spearman correlation, Default = pearson
#' @param return.table boolean value. TRUE return a data.frame instead a plot. Default = FALSE
#'
#' @return plot of correlation or (return.table = TRUE) data.frame of expression associated to beta methylation values
#'
#' @export

associateEnh2Exprs <- function(annotatedEnhancers, expression, hg, enhancer.db = 'FANTOM5', col.genes = 0, col.stat = 6, stat.thr = 0.05, col.logFC = 2, logfc.thr = .5, convert.genes = FALSE, convert.from, beta.thr = .3, overlap.param.thr = 40, param.type = "overlap.percentage", line.col = 'lightgrey', lmfit.col1='red', lmfit.col2='green', pal = 'RdGy', plot.type = 'splitted', show.text = FALSE, cor.type = 'pearson', return.table = FALSE){

    libpath <- paste0(.libPaths()[1], "/metExplorer")

    if(param.type != 'overlap.percentage' ){
        warning('Consider to change value threshold if you are changing param.type')
    }


    # load db
    if(hg=='hg19'){
        hacer.all <- readRDS(paste0(libpath, "/data/hacer.all.enhancer.hg19.gr.RDS"))
    } else {
        hacer.all <- readRDS(paste0(libpath, "/data/hacer.all.enhancer.hg38.gr.RDS"))
    }

    # create df with enhancers
    toAssociate <- annotatedEnhancers
    toAssociate.wenhancer <- metExplorer::queryDatabase(toAssociate, hacer.all, return.table = T)
    colnames(toAssociate.wenhancer)[duplicated(colnames(toAssociate.wenhancer))] <- paste0("annot.", colnames(toAssociate.wenhancer)[duplicated(colnames(toAssociate.wenhancer))])


    if(param.type == 'overlap.percentage'){
        toAssociate.wenhancer <- toAssociate.wenhancer[toAssociate.wenhancer$enhancer.perc>=overlap.param.thr,]
    } else if (param.type == "overlap.length"){
        toAssociate.wenhancer <- toAssociate.wenhancer[toAssociate.wenhancer$overlap.width>=overlap.param.thr,]
    } else if (param.type == "dmr.length"){
        toAssociate.wenhancer <- toAssociate.wenhancer[toAssociate.wenhancer$width>=overlap.param.thr,]
    } else {
        stop("Impossible to find parameter for filtering")
    }

    if(enhancer.db == 'FANTOM5'){
        if(param.type == 'overlap.percentage'){
            toAssociate.wenhancer.working <- toAssociate.wenhancer[, c("beta", "annot.associated_gene.FANTOM5", "enhancer.perc", "seqnames","start", "end", "annot.seqnames", "annot.start", "annot.end")]
        } else if (param.type == "overlap.length"){
            toAssociate.wenhancer.working <- toAssociate.wenhancer[, c("beta", "annot.associated_gene.FANTOM5", "overlap.width", "seqnames","start", "end", "annot.seqnames", "annot.start", "annot.end")]
        } else {
            toAssociate.wenhancer.working <- toAssociate.wenhancer[, c("beta", "annot.associated_gene.FANTOM5", "width", "seqnames","start", "end", "annot.seqnames", "annot.start", "annot.end")]
        }
    } else if (enhancer.db == '4DGenome') { 
        if(param.type == 'overlap.percentage'){
            toAssociate.wenhancer.working <- toAssociate.wenhancer[, c("beta", "annot.associated_gene.4DGenome", "enhancer.perc", "seqnames","start", "end", "annot.seqnames", "annot.start", "annot.end")]
        } else if (param.type == "overlap.length"){
            toAssociate.wenhancer.working <- toAssociate.wenhancer[, c("beta", "annot.associated_gene.4DGenome", "overlap.width", "seqnames","start", "end", "annot.seqnames", "annot.start", "annot.end")]
        } else {
            toAssociate.wenhancer.working <- toAssociate.wenhancer[, c("beta", "annot.associated_gene.4DGenome", "width", "seqnames","start", "end", "annot.seqnames", "annot.start", "annot.end")]
        }
    }

    colnames(toAssociate.wenhancer.working)[2] <- "usethis"


    # gene expression
    if(!is.null(col.genes) & col.genes!=0){
        expression.genes <- as.character(expression[,col.genes])
    } else if (missing(col.genes) | col.genes==0) {
        expression.genes <- rownames(expression)
    }

    expression.stat <- expression[, col.stat]
    expression.logFC <- expression[, col.logFC]
    expression.parsed <- data.frame(stringsAsFactors = F, genes = expression.genes, logFC = expression.logFC, expression.stat)
    expression.parsed <- expression.parsed[expression.parsed$expression.stat < stat.thr, ]
    expression.parsed <- expression.parsed[!is.na(expression.parsed[,1]),]

    if(convert.genes){
        to.translate <- expression.parsed$genes
        annotations <- biomaRt::select(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75, keys = as.character(to.translate), columns = "SYMBOL", keytype = convert.from)
        translated <- annotations[,2][match(expression.parsed$genes, annotations[, match(convert.from, colnames(annotations))])]
        expression.parsed$genes <- translated
        expression.parsed$sourceID <- to.translate
    }



    # start association
    outTable = data.frame()
    toAssociate.geneRaw <- toAssociate.wenhancer.working[, "usethis"]
    toAssociate.geneVect <- unlist(strsplit(toAssociate.geneRaw, ","))
    toAssociate.geneVect.intersect <- intersect(toAssociate.geneVect, expression.parsed$genes)
    for (gene in toAssociate.geneVect.intersect) {
        indB <- grep(gene, toAssociate.geneRaw)

        if (length(indB) > 1) {
            index <- which(!duplicated(toAssociate.wenhancer.working[indB, 1]))
            if (length(index) > 1) {
                cur.beta <- cbind(toAssociate.wenhancer.working[indB[index], 1])
                cur.enh <- cbind(toAssociate.wenhancer.working[indB[index],c(7:9)])
                cur.dmrs <- cbind(toAssociate.wenhancer.working[indB[index],c(4:6)])
            } else {
                cur.beta <- toAssociate.wenhancer.working[indB[1], 1]
                cur.enh <- cbind(toAssociate.wenhancer.working[indB[1],c(7:9)])
                cur.dmrs <- cbind(toAssociate.wenhancer.working[indB[1],c(4:6)])
            }
        } else if (length(indB == 1)) {
            cur.beta <- toAssociate.wenhancer.working[indB, 1]
            cur.enh <- cbind(toAssociate.wenhancer.working[indB,c(7:9)])
            cur.dmrs <- cbind(toAssociate.wenhancer.working[indB,c(4:6)])
            
        }

        ind.expr <- which(expression.parsed$genes %in% gene)
        if (abs(expression.parsed$logFC[ind.expr]) >= logfc.thr) {
            curSlice <- cbind(target.gene = rep(gene, length(cur.beta)), beta = cur.beta, logFC = rep(expression.parsed$logFC[ind.expr], length(cur.beta)),cur.dmrs, cur.enh)
            outTable <- rbind(outTable, curSlice)
        }
    }
    # colnames(outTable)[1:3] <- c('gene', 'beta', 'logFC')
    outTable[, 2] <- as.numeric(as.character(outTable[, 2]))
    outTable[, 3] <- as.numeric(as.character(outTable[, 3]))


    if(return.table){
        score <- scale((outTable[, 2] * outTable[, 3]) * -1)
        ind.order <- order(score, decreasing = T)
        score <- score[ind.order]
        outTable <- outTable[ind.order, ]
        outTable$score <- score
        outTable <- outTable[,c("seqnames","start","end","annot.seqnames", "annot.start","annot.end", "target.gene", "beta", "logFC", "score")]
        return(outTable)


    } else {
        if(plot.type == 'simple'){
            score <- scale((outTable[, 2] * outTable[, 3]) * -1)
            palette <- hcl.colors(length(outTable[, 2]), pal)
            ind.order <- order(score, decreasing = T)
            score <- score[ind.order]
            outTable <- outTable[ind.order, ]

            par(mar=c(4,4,4,4))
            layout(cbind(matrix(1, ncol = 5, nrow = 6), matrix(c(3,3,2,2,3,3), ncol = 1, nrow = 6)))
            plot(x = outTable[, 2], y = outTable[, 3], xlab = "beta", ylab = "logFC", xlim = c(-1, 1), ylim = c(-max(abs(outTable$logFC)), max(abs(outTable$logFC))), col = palette, pch=19, main = 'beta - logFC correlation', cex = 1.3)
            if (show.text) {
                text(x = outTable[, 1], y = outTable[, 5], labels = outTable[, 2])
            }
            abline(v = 0, h = 0, col = line.col, lty = 2)
            abline(lm(as.numeric(outTable[, 3]) ~ as.numeric(outTable[, 2]) - 1), col = lmfit.col1, lty = 2, lwd = 2)
            correlation <- cor.test(as.numeric(outTable[, 3]), as.numeric(outTable[, 2]), method = cor.type)
            mtext(paste("R:", round(correlation$statistic, 3), " pval:", round(correlation$p.val, 3), sep = " "), side = 1, at = -0.8, line = 2.5, cex = 0.7)

            par(mar=c(2,2,2,2))
            gradient <- as.raster(matrix(palette, ncol = 1))
            plot(c(0, 2), c(min(score), max(score)), type = "n", axes = F, xlab = "", ylab = "", main = "Score")
            steps <- ((abs(max(score)) + abs(min(score))) / 4)
            text(x = 1.5, y = seq(min(score), max(score), steps), labels = round(seq(min(score), max(score), steps), 3))
            rasterImage(gradient, 0, min(score), 1, max(score))

        } else if(plot.type == 'splitted'){
            score <- scale((outTable[, 2] * outTable[, 3]) * -1)
            palette <- hcl.colors(length(outTable[, 2]), pal)
            ind.order <- order(score, decreasing = T)
            score <- score[ind.order]
            outTable <- outTable[ind.order, ]

            indPos <- outTable[, 2] > 0
            indNeg <- outTable[, 2] < 0

            par(mar = c(4, 4, 4, 4))
            layout(cbind(matrix(1, ncol = 5, nrow = 6), matrix(c(3, 3, 2, 2, 4, 4), ncol = 1, nrow = 6)))

            if (length(which(indPos)) > 0) {
                plot(x = outTable[indPos, 2], y = outTable[indPos, 3], xlab = "beta", ylab = "logFC", xlim = c(-1, 1), ylim = c(-max(abs(outTable$logFC)), max(abs(outTable$logFC))), pch = 19, cex = 1.3, col = palette[indPos])
                if (show.text) {
                    text(x = outTable[indPos, 1], y = outTable[indPos, 5], labels = outTable[indPos, 2])
                }
                abline(lm(as.numeric(outTable[indPos, 3]) ~ as.numeric(outTable[indPos, 2]) - 1), col = lmfit.col1, lty = 2, lwd = 2)
                par(new = T)
                correlationPos <- cor.test(as.numeric(outTable[indPos, 3]), as.numeric(outTable[indPos, 2]), method = cor.type)
                mtext(paste("positive R:", round(correlationPos$statistic, 3), " pval:", round(correlationPos$p.val, 3), sep = " "), side = 1, at = -0.8, line = 2.5, cex = 0.7)
            } else { 
                correlationPos = ''
            }
            if (length(which(indNeg)) > 0) {
                if(sum(indPos) > 2){
                    par(new=T)
                }
                plot(x = outTable[indNeg, 2], y = outTable[indNeg, 3], xlab = "beta", ylab = "logFC", xlim = c(-1, 1), ylim = c(-max(abs(outTable$logFC)), max(abs(outTable$logFC))), pch = 19, cex = 1.3, col = palette[indNeg])
                if (show.text) {
                    text(x = outTable[indNeg, 1], y = outTable[indNeg, 5], labels = outTable[indNeg, 2])
                }
                abline(lm(as.numeric(outTable[indNeg, 3]) ~ as.numeric(outTable[indNeg, 2]) - 1), col = lmfit.col2, lty = 2, lwd = 2)
                correlationNeg <- cor.test(as.numeric(outTable[indNeg, 3]), as.numeric(outTable[indNeg, 2]), method = cor.type)
                mtext(paste("negative R:", round(correlationNeg$statistic, 3), " pval:", round(correlationNeg$p.val, 3), sep = " "), side = 1, at = 0.8, line = 2.5, cex = 0.7)
            } else { correlationNeg = ''}
            abline(v = 0, h = 0, col = line.col, lty = 2)
            
            
            par(mar = c(2, 2, 2, 2))
            gradient <- as.raster(matrix(palette, ncol = 1))
            plot(c(0, 2), c(min(score), max(score)), type = "n", axes = F, xlab = "", ylab = "", main = "Score")
            steps <- ((abs(max(score)) + abs(min(score))) / 4)
            text(x = 1.5, y = seq(min(score), max(score), steps), labels = round(seq(min(score), max(score), steps), 3))
            rasterImage(gradient, 0, min(score), 1, max(score))

            plot.new()
            par(mar = c(0, 0, 0, 0), xpd = T)
            legend("center", c("Beta Positive", "Beta Negative"), pch = 15, col = c(lmfit.col1, lmfit.col2), box.lwd = 0, cex = 1.2)
        }

    }
}
