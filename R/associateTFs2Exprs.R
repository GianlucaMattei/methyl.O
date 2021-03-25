#' Associate target genes to TFs and retrieve their expression
#'
#' For each TF find and associate the target genes, within the annotation results, and retrieve expression
#'
#' @param annotatedDMRs anotated DMRs list resultingfrom annotateDMRs() or scoreAnnotatedDMRs()
#' @param expressionProfile expression data.frame
#' @param active.features character vectors containing features to correlate. Must be from names of resulting list from annotateDMRs. Additional feature names can be first exons (exons1) or first intron (intron1). To use more than one feature use c(). Default = c("promoters", "heads")
#' @param col.genes numeric, the column of expressionProfile data.frame with gene Ids. If NULL geneIDs will be taken from rownames() of expressionProfile. Default = 0.
#' @param col.stat numeric, the column of expressionProfile data.frame with the statistics to use. Default = 6.
#' @param stat.thr numeric, threshold for statistical significance. Default = 0.05
#' @param col.logFC numeric, the column of expressionProfile data.frame with log. fold change. Default = 2
#' @param logfc.thr numeric, threshold value for log. fold change. Default = 0.
#' @param convert.genes logical, used to indicate if gene ids have to be translated in official gene symbols. Default = FALSE
#' @param convert.from character, used annotation for gene in expressionProfile to be converted to symbols gene IDs. Accepted: c("ENTREZID" ,"EXONID" ,"GENEBIOTYPE" ,"GENEID" ,"GENENAME" ,"PROTDOMID" ,"PROTEINDOMAINID" ,"PROTEINDOMAINSOURCE" ,"PROTEINID" ,""SEQSTRAND" ,"SYMBOL" ,"TXBIOTYPE" ,"TXID" ,"TXNAME" ,"UNIPROTID")
#' @param beta.thr numeric, beta difference threshold value. Default = 0.3.
#' @param param.type character, threshold parameter to filter methylations overlapping the selected features. Accetped c("dmr.length", "overlap.length", "overlap.percentage"). Default = "overlap.length".
#' @param overlap.param.thr numeric, threshold value for selected parameter to filter methylations overlapping the selected features. Default = 100
#'
#' @return data.frame with TF methylation levels, target.genes expression
#'
#' @export


associateTFs2Exprs <- function(annotatedDMRs, expressionProfile, active.features = c("promoters", "heads"), col.genes = 0, col.stat = 6, stat.thr = 0.05, col.logFC = 2, logfc.thr = 0, convert.genes = FALSE, convert.from,  beta.thr = .3,  overlap.param.thr = 30, param.type = 'overlap.percentage'){
    
    libpath <- paste0(.libPaths()[1], "/methyl.O") 
    tf.genes <- readRDS(paste0(libpath,"/data/ENCODE_targets.RDS"))

    DMRs.parsed <- lapply(annotatedDMRs, function(x){
        ind <- x$TF == 1
        x[ind,]
    })
    
    outTable = data.frame()
    scale1 <- function(x) {
        (x - min(x)) / (max(x) - min(x))
    }

    if (any(active.features == "exons1")) {
        DMRs.parsed$exons1 <- DMRs.parsed$exons[DMRs.parsed$exons$rank == 1, ]
    }

    if (any(active.features == "introns1")) {
        DMRs.parsed$introns1 <- DMRs.parsed$introns[DMRs.parsed$introns$intron.rank == 1, ]
    }

    # create objects to work on
    # df with active.features
    toAssociate <- c()

    if(param.type == 'overlap.percentage'){

        for(f in active.features){
            ind <- match(f, names(DMRs.parsed))
            ind.perc <- grep('perc', colnames(DMRs.parsed[[ind]]))
            curFeat <- DMRs.parsed[[f]][, c("beta", "symbol", colnames(DMRs.parsed[[f]])[ind.perc])]
            curFeat$feature <- rep(names(DMRs.parsed)[ind],nrow(curFeat))
            colnames(curFeat)=c('beta', 'symbol', 'perc.overlap', 'feature')
            toAssociate <- rbind(toAssociate, curFeat)
        }

        toAssociate <- toAssociate[toAssociate$perc.overlap >= overlap.param.thr, ]
    } else if (param.type == "overlap.length"){

        for(f in active.features){
            ind <- match(f, names(DMRs.parsed))
            ind.owidth <- grep('overlap\\.width', colnames(DMRs.parsed[[ind]]))
            curFeat <- DMRs.parsed[[f]][, c("beta", "symbol", colnames(DMRs.parsed[[f]])[ind.owidth])]
            curFeat$feature <- rep(names(DMRs.parsed)[ind],nrow(curFeat))
            colnames(curFeat)=c('beta', 'symbol', 'width.overlap', 'feature')
            toAssociate <- rbind(toAssociate, curFeat)
        }

        toAssociate <- toAssociate[toAssociate$width.overlap >= overlap.param.thr, ]

    } else if (param.type == "dmr.length"){
        
        for(f in active.features){
            ind <- match(f, names(DMRs.parsed))
            ind.width <- 4
            curFeat <- DMRs.parsed[[f]][, c("beta", "symbol", colnames(DMRs.parsed[[f]])[ind.width])]
            curFeat$feature <- rep(names(DMRs.parsed)[ind],nrow(curFeat))
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


    # gene expressionProfile
    expressionProfile <- expressionProfile[abs(expressionProfile[, col.logFC])>= logfc.thr,]
    expressionProfile <- expressionProfile[expressionProfile[, col.stat] < stat.thr,]

    if(col.genes!=0){
        expression.genes <- as.character(expressionProfile[,col.genes])
    } else {
        expression.genes <- rownames(expressionProfile)
    }

    if (all(expression.genes[1:6] == c("1", "2", "3", "4", "5", "6"))) {
        stop("Please select the genes column \n retrieved genes are not symbols")
    }


    if(convert.genes){
        to.translate <- expression.genes
        annotations <- biomaRt::select(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75, keys = as.character(to.translate), columns = "SYMBOL", keytype = convert.from)
        translated <- annotations[,2][match(expression.genes, annotations[, match(convert.from, colnames(annotations))])]
        ind <- which(is.na(translated))
        expression.genes <- translated[-ind]
        expressionProfile <- expressionProfile[-ind,]
    }

    toAssociate <- toAssociate[abs(toAssociate$beta)>= beta.thr ,]
    c <- toAssociate[toAssociate[,2] %in% tf.genes[,1],]

    # associate each TF to genes
    cur.TF = list()
    for(i in 1:nrow(toAssociate)){
        x <- toAssociate[i,]
        ind.hold <- tf.genes[,1] %in% x['symbol']
        associated.TF.genes <- intersect(unique(tf.genes[ind.hold,2]),expression.genes)
        associated.TF.expression <- expressionProfile[expression.genes %in% associated.TF.genes, col.logFC]
        rownames(x) = NULL
        if(length(associated.TF.expression)!=0){
            cur.TF[[i]] <- cbind(x, target.expression = associated.TF.expression, target.symbol = associated.TF.genes)
        }

    }
    associated.table <- do.call('rbind', cur.TF)    
    return(associated.table)
}
