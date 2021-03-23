#' Associate target genes to TFs and retrieve their expression
#'
#' For each TF find and associate the target genes, within the annotation results, and retrieve expression
#'
#' @param results List. Corresponding to resulting list from getMetAnnotations() or scoreMets()
#' @param features feature to correlate. Must be from names of resulting list from gemMetAnnotations.  Additional feature names can be first exons (exons1) or first intron (intron1). To use more than one feature use c().
#' @param expression expression data.frame
#' @param col.genes number of the column of expression data.frame with gene Ids. If NULL geneIDs will be taken from rownames() of expression input file
#' @param col.stat number of the column of expression data.frame with the selected statistics
#' @param stat.thr threshold for statistical significance. Default = 0.05
#' @param col.logFC number of the column of expression data.frame with lofFC
#' @param logfc.thr logFC threshold value
#' @param convert.genes boolean value used to indicate if gene ids have to be translated in official gene symbols
#' @param convert.from gene id annotation of expression data.frame to convert to symbols gene ids. Character: ENTREZID TXNAME TXID GENEID
#' @param beta.thr beta difference threshold value
#' @param param.type threshold parameter to filter methylations overlapping the selected features. Accetped c("dmr.length", "overlap.length", "overlap.percentage"). Default = "overlap.length".
#' @param overlap.param.thr threshold value for selected parameter to filter methylations overlapping the selected features. Default = 100
#'
#' @return data.frame with TF, methylation leveles, target.genes and their expression
#'
#' @export


associateTFs2Exprs <- function(results, features = c("promoters", "heads"), expression, col.genes = 0, col.stat = 6, stat.thr = 0.05, col.logFC = 2, logfc.thr = 0, convert.genes = FALSE, convert.from,  beta.thr = .3,  overlap.param.thr = 30, param.type = 'overlap.percentage'){
    
    libpath <- paste0(.libPaths()[1], "/metExplorer") 
    tf.genes <- readRDS(paste0(libpath,"/data/ENCODE_targets.RDS"))

    results.parsed <- lapply(results, function(x){
        ind <- x$TF == 1
        x[ind,]
    })
    
    outTable = data.frame()
    scale1 <- function(x) {
        (x - min(x)) / (max(x) - min(x))
    }

    if (any(features == "exons1")) {
        results.parsed$exons1 <- results.parsed$exons[results.parsed$exons$rank == 1, ]
    }

    if (any(features == "introns1")) {
        results.parsed$introns1 <- results.parsed$introns[results.parsed$introns$intron.rank == 1, ]
    }

    # create objects to work on
    # df with features
    toAssociate <- c()

    if(param.type == 'overlap.percentage'){

        for(f in features){
            ind <- match(f, names(results.parsed))
            ind.perc <- grep('perc', colnames(results.parsed[[ind]]))
            curFeat <- results.parsed[[f]][, c("beta", "symbol", colnames(results.parsed[[f]])[ind.perc])]
            curFeat$feature <- rep(names(results.parsed)[ind],nrow(curFeat))
            colnames(curFeat)=c('beta', 'symbol', 'perc.overlap', 'feature')
            toAssociate <- rbind(toAssociate, curFeat)
        }

        toAssociate <- toAssociate[toAssociate$perc.overlap >= overlap.param.thr, ]
    } else if (param.type == "overlap.length"){

        for(f in features){
            ind <- match(f, names(results.parsed))
            ind.owidth <- grep('overlap\\.width', colnames(results.parsed[[ind]]))
            curFeat <- results.parsed[[f]][, c("beta", "symbol", colnames(results.parsed[[f]])[ind.owidth])]
            curFeat$feature <- rep(names(results.parsed)[ind],nrow(curFeat))
            colnames(curFeat)=c('beta', 'symbol', 'width.overlap', 'feature')
            toAssociate <- rbind(toAssociate, curFeat)
        }

        toAssociate <- toAssociate[toAssociate$width.overlap >= overlap.param.thr, ]

    } else if (param.type == "dmr.length"){
        
        for(f in features){
            ind <- match(f, names(results.parsed))
            ind.width <- 4
            curFeat <- results.parsed[[f]][, c("beta", "symbol", colnames(results.parsed[[f]])[ind.width])]
            curFeat$feature <- rep(names(results.parsed)[ind],nrow(curFeat))
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


    # gene expression
    expression <- expression[abs(expression[, col.logFC])>= logfc.thr,]
    expression <- expression[expression[, col.stat] < stat.thr,]

    if(col.genes!=0){
        expression.genes <- as.character(expression[,col.genes])
    } else {
        expression.genes <- rownames(expression)
    }

    if(convert.genes){
        to.translate <- expression.genes
        annotations <- biomaRt::select(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75, keys = as.character(to.translate), columns = "SYMBOL", keytype = convert.from)
        translated <- annotations[,2][match(expression.genes, annotations[, match(convert.from, colnames(annotations))])]
        ind <- which(is.na(translated))
        expression.genes <- translated[-ind]
        expression <- expression[-ind,]
    }

    toAssociate <- toAssociate[abs(toAssociate$beta)>= beta.thr ,]
    toAssociate <- toAssociate[toAssociate[,2] %in% tf.genes[,1],]

    # associate each TF to genes
    cur.TF = list()
    for(i in 1:nrow(toAssociate)){
        x <- toAssociate[i,]
        ind.hold <- tf.genes[,1] %in% x['symbol']
        associated.TF.genes <- intersect(unique(tf.genes[ind.hold,2]),expression.genes)
        associated.TF.expression <- expression[expression.genes %in% associated.TF.genes, col.logFC]
        rownames(x) = NULL
        if(length(associated.TF.expression)!=0){
            cur.TF[[i]] <- cbind(x, target.expression = associated.TF.expression, target.symbol = associated.TF.genes)
        }

    }
    associated.table <- do.call('rbind', cur.TF)    
    return(associated.table)
}
