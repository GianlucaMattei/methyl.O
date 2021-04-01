#' Compute expression methylation correlation 
#'
#' Compute and plot correlation between expression and methylation
#'
#' @param annotatedEnhancers data.frame. It corresponds to resulting list from annotateEnhancers()
#' @param expressionProfile expression profile data.frame
#' @param hg character, "hg19", "hg38". Version of the enhancer database. Default = "hg19"
#' @param enhancer.db character, which database to use between 'FANTOM5' or '4DGenome'. Default = "FANTOM5"
#' @param col.genes numeric, the column of expressionProfile data.frame with gene Ids. If NULL geneIDs will be taken from rownames() of expressionProfile. Default = 0.
#' @param col.stat numeric, the column of expressionProfile data.frame with the statistics to use. Default = 6.
#' @param stat.thr threshold for statistical significance. Default = 0.05
#' @param col.logFC numeric, the column of expressionProfile data.frame with log. fold change. Default = 2
#' @param logfc.thr numeric, threshold value for log. fold change. Default = 0.
#' @param convert.genes logical, used to indicate if gene ids have to be translated in official gene symbols. Default = FALSE
#' @param convert.from character, used annotation for gene in expressionProfile to be converted to symbols gene IDs. Accepted: c("ENTREZID" ,"EXONID" ,"GENEBIOTYPE" ,"GENEID" ,"GENENAME" ,"PROTDOMID" ,"PROTEINDOMAINID" ,"PROTEINDOMAINSOURCE" ,"PROTEINID" ,""SEQSTRAND" ,"SYMBOL" ,"TXBIOTYPE" ,"TXID" ,"TXNAME" ,"UNIPROTID")
#' @param beta.thr numeric, beta difference threshold value. Default = 0.3.
#' @param param.type character, threshold parameter to filter methylations overlapping the selected features. Accetped c("dmr.length", "overlap.length", "overlap.percentage"). Default = "overlap.length".
#' @param show.text logical, indicating if print gene names in the final plot. Accepted values: TRUE or FALSE. Default = FALSE.
#' @param overlap.param.thr numeric, threshold value for selected parameter to filter methylations overlapping the selected features. Default = 100
#' @param line.col character, color of lines at x=0, y=0. Defalt = "lightgray"
#' @param lmfit.col1 character, color of linear model line 1 or for simple plot. Defalt = "red"
#' @param lmfit.col2 character, color of linear model line 2. Defalt = "green"
#' @param pal character, color palette. hcl.pals() to show available. Default = "RdGy"
#' @param plot.type  character, compute or not different linear models for upregulated and downregulated genes. Accepted: "simple" or "splitted". Default = "splitted".
#' @param cor.type character, correlation method. Available "pearson", "kendall" or "spearman" correlation, Default = "pearson"
#' @param return.table logical, TRUE return a data.frame instead a plot. Default = FALSE
#'
#' @return plot of correlation or (if return.table = TRUE) a data.frame of genes expression associated to beta methylation values of enhancers.
#'
#' @export

annotatedEnh2Exprs <- function(annotatedEnhancers, expressionProfile, hg = 'hg19', enhancer.db = 'FANTOM5', col.genes = 0, col.stat = 6, stat.thr = 0.05, col.logFC = 2, logfc.thr = .5, convert.genes = FALSE, convert.from, beta.thr = .3, overlap.param.thr = 40, param.type = "overlap.percentage", line.col = 'lightgrey', lmfit.col1='red', lmfit.col2='green', pal = 'RdGy', plot.type = 'splitted', show.text = FALSE, cor.type = 'pearson', return.table = FALSE){

    libpath <- paste0(.libPaths()[1], "/methyl.O")

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
    toAssociate.wenhancer <- methyl.O::queryDatabase(toAssociate, hacer.all, return.table = T)
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

    # gene expressionProfile
    if(!is.null(col.genes) & col.genes!=0){
        expression.genes <- as.character(expressionProfile[,col.genes])
    } else if (missing(col.genes) | col.genes==0) {
        expression.genes <- rownames(expressionProfile)
    }

    if (all(expression.genes[1:6] == c("1", "2", "3", "4", "5", "6"))) {
        stop("Please select the genes column \n retrieved genes are not symbols")
    }

    expression.stat <- expressionProfile[, col.stat]
    expression.logFC <- expressionProfile[, col.logFC]
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

                lm.line <- lm(as.numeric(outTable[indPos, 3]) ~ as.numeric(outTable[indPos, 2]) + 0)
                x1 <- 1.5
                y1 <- lm.line[1]$coefficients[[1]] * x1
                segments(0, 0, x1, y1, col = lmfit.col1, lty = 2, lwd = 2)
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

                lm.line <- lm(as.numeric(outTable[indNeg, 3]) ~ as.numeric(outTable[indNeg, 2]) + 0)
                x1 <- -1.5
                y1 <- lm.line[1]$coefficients[[1]] * x1
                segments(0, 0, x1, y1, col = lmfit.col2, lty = 2, lwd = 2)

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
            legend("center", c("Beta + Corr", "Beta - Corr"), pch = 15, col = c(lmfit.col1, lmfit.col2), box.lwd = 0, cex = 1.2)
        }

    }
}
