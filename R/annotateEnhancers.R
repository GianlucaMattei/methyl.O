#' Query database to find enhancer.
#'
#' Query database of annotations results list to enhancer database in order to associate differentially methylated segments to genes.
#'
#' @param DMRsRanges the DMRs ranges, it must have the following columns: chr, start, end, beta diff. Other columns will be stored in the resulting output under the column other.
#' @param hg character, Available "hg19", "hg38". Genome assembly version. Default = "hg19"
#' @param thr.beta numeric, beta difference threshold to consider methylations. Default = 0.3
#' @param param.type character, threshold parameter to filter methylations overlapping the selected features. Accetped c("dmr.length", "overlap.length", "overlap.percentage"). Default = "overlap.length".
#' @param overlap.param.thr numeric, threshold value for selected parameter to filter methylations overlapping the selected features. Default = 100
#' @param score.modifier numeric, value between 0-1. It specifies how the final score is computed by assigning different weights to the methylation charactersistics of enhacners or to genes already involved in pathologies. By increasing this value to 1, resulting scores will be focused on discovering segments affecting gene expression. A value equal to 0 will focus the results on enahcners involving genes associated to pathologies, not considering the effect of methylation. Default = 0.5
#' @param col.betadiff numeric, column position for beta diff. in input table. Default = 4
#'
#' @return a vector of presence or a data.frame
#'
#' @export

annotateEnhancers <- function(DMRsRanges, hg='hg19', thr.beta = 0.3,  overlap.param.thr = 40, param.type = "overlap.percentage", score.modifier = 0.5, col.betadiff = 4) {

    scale1 <- function(x) {
        (x - min(x)) / (max(x) - min(x))
    }

    libpath <- paste0(.libPaths()[1], "/methyl.O")

    if(param.type != 'overlap.percentage' ){
        warning('Consider to change value threshold if you are changing param.type')
    }


    if (hg == "hg19") {
        db.gr <- readRDS(paste0(libpath, "/data/hacer.all.enhancer.hg19.gr.RDS"))
    } else if (hg == "hg38") {
        db.gr <- readRDS(paste0(libpath, "/data/hacer.all.enhancer.hg38.gr.RDS"))
    }

    # loading database
    ncg <- read.table(sep = "\t", stringsAsFactors = F, header = T, paste0(libpath, "/data/pathsToGene_NCG.tsv"))
    cosmic <- readRDS(file = paste0(libpath, "/data/DB.COSMIC.RDS"))

    DMRs.gr <- GenomicRanges::makeGRangesFromDataFrame(formatDMRsInput(DMRsRanges, thr.beta, col.betadiff = col.betadiff), keep.extra.columns = T)
    hits <- GenomicRanges::findOverlaps(DMRs.gr, db.gr)
    DMRs.gr.matched <- DMRs.gr[S4Vectors::queryHits(hits)]
    db.gr.matched <- db.gr[S4Vectors::subjectHits(hits)]
    overlap.width <- S4Vectors::width(GenomicRanges::pintersect(DMRs.gr.matched, db.gr.matched))
    db.gr.matched$enhancer.perc <- overlap.width / S4Vectors::width(db.gr.matched) * 100
    db.matched <- data.frame(db.gr.matched)
    colnames(db.matched)[1:5] <- paste0("annot.", colnames(db.matched)[1:5])
    DMRs.annotated <- cbind(data.frame(DMRs.gr.matched), db.matched[, c(1:4, 8, 9, 10, 11, 14, 17, 24:26)])
    DMRs.annotated$overlap.width <- overlap.width

    # FANTOM5
    genes.query.list <- apply(DMRs.annotated, 1, function(x) {
        unique(strsplit(x['annot.associated_gene.FANTOM5'], ",")[[1]])
    })
    # FANTOM5 - NCG
    ncg.hits <- lapply(genes.query.list, function(x) {
        x.names <- x[[1]][which(x[[1]] %in% ncg$symbol)]
        x.hits <- length(x.names[!is.na(x.names)])
        if (x.hits == 0) {
            x.names <- NA
        }

        return(cbind(symbol = x.names, hits = x.hits))
    })

    FANTOM5.NCG <- unlist(lapply(ncg.hits, function(x) {
        paste0(x[, 1], collapse = ",")
    }))
    FANTOM5.NCG.hits <- unlist(lapply(ncg.hits, function(x) {
        x[1, 2]
    }))

    # FANTOM5 - COSMIC
    cosmic.hits <- lapply(genes.query.list, function(x) {
        x.names <- x[[1]][which(x[[1]] %in% cosmic$symbol)]
        x.hits <- length(x.names[!is.na(x.names)])
        if (x.hits == 0) {
            x.names <- NA
        }

        return(cbind(symbol = x.names, hits = x.hits))
    })

    FANTOM5.COSMIC <- unlist(lapply(cosmic.hits, function(x) {
        paste0(x[, 1], collapse = ",")
    }))
    FANTOM5.COSMIC.hits <- unlist(lapply(cosmic.hits, function(x) {
        x[1, 2]
    }))


    # GENOME4D
    genes.query.list <- apply(DMRs.annotated, 1, function(x) {
        unique(strsplit(x['annot.associated_gene.4DGenome'], ",")[[1]])
    })
    # GENOME4D - NCG
    ncg.hits <- lapply(genes.query.list, function(x) {
        x.names <- x[[1]][which(x[[1]] %in% ncg$symbol)]
        x.hits <- length(x.names[!is.na(x.names)])
        if (x.hits == 0) {
            x.names <- NA
        }

        return(cbind(symbol = x.names, hits = x.hits))
    })

    GENOME4D.NCG <- unlist(lapply(ncg.hits, function(x) {
        paste0(x[, 1], collapse = ",")
    }))
    GENOME4D.NCG.hits <- unlist(lapply(ncg.hits, function(x) {
        x[1, 2]
    }))

    # FANTOM5 - COSMIC
    cosmic.hits <- lapply(genes.query.list, function(x) {
        x.names <- x[[1]][which(x[[1]] %in% cosmic$symbol)]
        x.hits <- length(x.names[!is.na(x.names)])
        if (x.hits == 0) {
            x.names <- NA
        }
        return(cbind(symbol = x.names, hits = x.hits))
    })

    GENOME4D.COSMIC <- unlist(lapply(cosmic.hits, function(x) {
        paste0(x[, 1], collapse = ",")
    }))
    
    GENOME4D.COSMIC.hits <- unlist(lapply(cosmic.hits, function(x) {
        x[1, 2]
    }))

    colhold <- c("seqnames","start","end","width","beta","annot.seqnames","annot.start","annot.end","annot.width","annot.associated_gene.FANTOM5","annot.associated_gene.4DGenome","enhancer.perc", "overlap.width", "others")
    DMRs.annotated <- cbind(DMRs.annotated[,colhold], FANTOM5.NCG, FANTOM5.NCG.hits, FANTOM5.COSMIC, FANTOM5.COSMIC.hits,  GENOME4D.NCG, GENOME4D.NCG.hits, GENOME4D.COSMIC, GENOME4D.COSMIC.hits)

    score.modifier.beta <- 1 - score.modifier
    score.modifier.alpha <- score.modifier
    DMRs.annotated$score <- (scale1(
        DMRs.annotated$beta * DMRs.annotated$enhancer.perc
    ) * score.modifier.alpha) + (scale1(
        as.numeric(DMRs.annotated[,'FANTOM5.NCG.hits']) + 
        as.numeric(DMRs.annotated[,'FANTOM5.COSMIC.hits']) + 
        as.numeric(DMRs.annotated[,'GENOME4D.NCG.hits']) +
        as.numeric(DMRs.annotated[,'GENOME4D.COSMIC.hits'])        
    ) * score.modifier.beta)
    
    DMRs.annotated <- DMRs.annotated[abs(as.numeric(as.character(DMRs.annotated$beta))) >= thr.beta,]


    if(param.type == 'overlap.percentage'){
        DMRs.annotated <- DMRs.annotated[DMRs.annotated$enhancer.perc>=overlap.param.thr,]
    } else if (param.type == "overlap.length"){
        DMRs.annotated <- DMRs.annotated[DMRs.annotated$overlap.width>=overlap.param.thr,]
    } else if (param.type == "dmr.length"){
        DMRs.annotated <- DMRs.annotated[DMRs.annotated$width>=overlap.param.thr,]
    } else {
        stop("Impossible to find parameter for filtering")
    }

    # duplicated lines are returned due different cells lines:
    DMRs.annotated <- DMRs.annotated[!duplicated(DMRs.annotated), ]

    DMRs.annotated <- DMRs.annotated[order(DMRs.annotated$score, decreasing = T), ]
    return(DMRs.annotated)
}
