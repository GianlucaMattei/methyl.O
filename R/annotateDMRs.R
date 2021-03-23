#' Annotates the methylation profiles. 
#' 
#' Annotate the differentially methylated regions 
#'  
#' @param table the input table, it must have the following columns: chr, start, end, beta diff. Other columns will be stored in the resulting output under the column other.
#' @param prom.length length of promoters
#' @param head.length length of the first part of the gene, named head, starting from the TSS. If longer than the gene, the entire txs will be considered as head. 
#' @param longest.trx option to use the longest transcript to represent the gene
#' @param annotation Ensembl or UCSC
#' @param hg genome assembly version
#' @param annotation.fast 1:1 mapping or 1:many - many:1 - many:many mapping
#' @param thr.beta the beta difference threshold to consider methylations. Default = 0.3
#' @param thr.cgis the length, in percentage, of methylated CGIs in order to be considered altered. Default = 0.4
#' @param thr.enhancer the length in percentage of methylated enhancer in order to be considered altered. Default = 0.4
#' @param col.betadiff the column position for beta diff. in input table. Default = 4
#' @param col.beta1 the column position for first sample beta values in input table
#' @param col.beta2 the column position for second sample beta values in input table
#' 
#' @return list of annotated table
#' 
#' @export

annotateDMRs = function(table , prom.length=1500, head.length=1500, longest.trx=TRUE, annotation='ensembl', hg='hg19', annotation.fast=TRUE, thr.beta=.3, thr.cgis=.4, thr.enhancer=.4, col.betadiff = 4, col.beta1 = NULL, col.beta2 = NULL) {
    libpath <- paste0(.libPaths()[1], "/metExplorer")

    ncg <- read.table(sep = "\t", stringsAsFactors = F, header = T, paste0(libpath, "/data/pathsToGene_NCG.tsv"))
    cosmic <- readRDS(file = paste0(libpath,"/data/DB.COSMIC.RDS"))

    if (hg == "hg19") {
        cosmic.CGIs <- readRDS(paste0(libpath,"/data/CGIs.19.comsic.positions.gr.RDS"))
        cgis <- readRDS(paste0(libpath,"/data/CGIs.19.gr.RDS"))
        hacer <- readRDS(paste0(libpath,"/data/hacer.all.enhancer.hg19.gr.RDS"))
        gnomad <- readRDS(file = paste0(libpath,"/data/DB.gnomad.hg19.RDS"))
        dgv <- readRDS(file = paste0(libpath,"/data/DB.DGV.hg19.RDS"))

    } else if (hg == "hg38") {
        cosmic.CGIs <- readRDS(paste0(libpath,"/data/CGIs.38.comsic.positions.gr.RDS"))
        cgis <- readRDS(paste0(libpath,"/data/CGIs.38.gr.RDS"))
        hacer <- readRDS(paste0(libpath,"/data/hacer.all.enhancer.hg19.gr.RDS"))
        gnomad <- readRDS(file = paste0(libpath,"/data/DB.gnomad.hg38.RDS"))
        dgv <- readRDS(file = paste0(libpath,"/data/DB.DGV.hg38.RDS"))
    }


    if(annotation == 'ensembl'){
        if(hg=='hg19'){
            genedb <- EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
        } else if (hg=='hg38') {
            genedb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
        }
        ktype.genes <- "GENEID"
        ktype.txs <- 'TXID'
        cosmic$ID <- cosmic$ensemblID

    } else if (annotation == 'ucsc') {
        if(hg=='hg19'){
            genedb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
        } else if (hg=='hg38') {
            stop('hg38 not supported for UCSC')
        }

        ktype.genes <- "ENTREZID"
        ktype.txs <- 'TXNAME'
        cosmic$ID <- cosmic$entrezID
    } else {
        stop('annotation not supported. Valid annotations: ensembl / ucsc')
    }

    results <- list()

    table <- formatDMRsInput(table, thr.beta = .3, col.betadiff = col.betadiff, col.beta1 = col.beta1, col.beta2 = col.beta2)

    # creating granges
    table.ranges <- GenomicRanges::makeGRangesFromDataFrame(table, keep.extra.columns = TRUE)

    # longest transcipts
    if(annotation == 'ensembl'){
        if(longest.trx==TRUE){
            txs.selected <- readRDS(file = paste0(libpath,"/data/txsSelectedLongest_ENSEMBL_",hg ,".RDS"))
        } else {    
            txs.selected <- readRDS(file = paste0(libpath,"/data/txsSelectedAll_ENSEMBL_", hg, ".RDS"))
        }
    }
    
    if(annotation == 'ucsc'){
        if(longest.trx==TRUE){
            txs.selected <- readRDS(file = paste0(libpath,"/data/txsSelectedLongest_UCSC_", hg, ".RDS"))
        } else {    
            txs.selected <- readRDS(file = paste0(libpath,"/data/txsSelectedAll_UCSC_", hg, ".RDS"))
        }
    }

    # creating database for all genes in genome
    db.genes.base <- GenomicRanges::makeGRangesFromDataFrame(txs.selected, keep.extra.columns = T)

    # mapping on genes considering promoters too
    genes.coordinates <- cbind(seqnames = GenomicRanges::seqnames(db.genes.base), data.frame(GenomicRanges::ranges(db.genes.base))[, -4], gene.id = data.frame(db.genes.base)[, "gene.id"], tx.name = data.frame(db.genes.base)[, "tx.name"], strand = GenomicRanges::strand(db.genes.base))
    strand.map <- as.logical(ifelse(genes.coordinates[, "strand"] == "+", 1, 0))
    genes.coordinates.withPromoters <- genes.coordinates
    genes.coordinates.withPromoters[which(strand.map), "start"] <- genes.coordinates[which(strand.map), "start"] - prom.length
    genes.coordinates.withPromoters[-which(strand.map), "end"] <- genes.coordinates[-which(strand.map), "end"] + prom.length
    db.genes <- GenomicRanges::makeGRangesFromDataFrame(genes.coordinates.withPromoters, keep.extra.columns = T)

    # overlap annotations and Mets
    annotated.table <- data.frame()
    default.colnames <- c("seqnames", "start", "end", "width", "beta", "others", "gene.start", "gene.end", "gene.width", "gene.strand", "gene.id", "tx.name", "genes.perc", "tag")
    hits <- GenomicRanges::findOverlaps(db.genes, table.ranges)
    table.overlapping <- data.frame(table.ranges[S4Vectors::subjectHits(hits), ], stringsAsFactors = F)[, -5]
    db.genes.overlapping <- data.frame(db.genes[S4Vectors::queryHits(hits), ], stringsAsFactors = F)
    if(length(hits)>0){
        colnames(db.genes.overlapping)[2:6] <- c("gene.start", "gene.end", "gene.width", "gene.strand", "gene.id")
        annotated.table <- cbind(table.overlapping, db.genes.overlapping[, c("gene.start","gene.end", "gene.width","gene.strand","gene.id","tx.name")])
        overlaps <- GenomicRanges::pintersect(db.genes[S4Vectors::queryHits(hits), ], table.ranges[S4Vectors::subjectHits(hits), ])
        width.overlap <- S4Vectors::width(overlaps)
        annotated.table <- cbind(annotated.table, genes.perc = width.overlap / S4Vectors::width(db.genes[S4Vectors::queryHits(hits), ]) * 100)
        annotated.table$tag <- apply(annotated.table,1,function(x){gsub(' +', '', paste0(x[1:3], collapse='_'))})
        annotated.table.ranges <- GenomicRanges::makeGRangesFromDataFrame(annotated.table, keep.extra.columns = TRUE)
    }

    results$genes <- annotated.table

    # databases overlap and database-score
    if (nrow(results$genes) > 0) {
        results$genes$dgv <- metExplorer::queryDatabase(results$genes, dgv, return.table = FALSE)
        results$genes$gnomad <- metExplorer::queryDatabase(results$genes, gnomad, return.table = FALSE)
        results$genes$NCG <- NA
        results$genes$NCG_type <- NA
        results$genes$COSMIC <- as.numeric(results$genes$gene.id %in% cosmic$ID)
        results$genes$cosmic.CGIs <- metExplorer::queryDatabase(results$genes, cosmic.CGIs, return.table = FALSE, is.genomic.ranges = TRUE,thr=thr.cgis)
        results$genes$hacer <- metExplorer::queryDatabase(results$genes, hacer, return.table = FALSE, thr = thr.enhancer)
        results$genes$CGIs <- metExplorer::queryDatabase(results$genes, cgis, return.table = FALSE, ,thr=thr.cgis)
    }


    # move "others" col to the end
    if(any(match('others', colnames(results$genes)))){
        results$genes <- cbind(results$genes[,-match('others', colnames(results$genes))],
            others = results$genes[,match('others', colnames(results$genes))]
            )
    }

    ##### GENE'S HEAD
    # extract heads 
    ind.hold <- which(txs.selected$width <= head.length)
    txs.selected.head.excluded <- txs.selected[ind.hold,]
    txs.selected.head.included <- txs.selected[-ind.hold,]

    txs.selected.head.pos <- txs.selected.head.included[txs.selected.head.included$strand=='+',]
    txs.selected.head.pos$end <- txs.selected.head.pos$start + head.length
    txs.selected.head.pos$width <- head.length
    txs.selected.head.neg <- txs.selected.head.included[txs.selected.head.included$strand=='-',]
    txs.selected.head.neg$start <- txs.selected.head.neg$end - head.length
    txs.selected.head.neg$width <-  head.length

    db.genes.head <- GenomicRanges::makeGRangesFromDataFrame(rbind(txs.selected.head.excluded, txs.selected.head.pos, txs.selected.head.neg), keep.extra.columns = T)

    # overlap annotations and Mets
    annotated.table.head <- data.frame()
    default.colnames <- c("seqnames", "start", "end", "width", "beta", "others", "head.start", "head.end", "head.width", "gene.strand", "gene.id", "tx.name", "overlap.perc.head", "tag")
    hits <- GenomicRanges::findOverlaps(db.genes.head, table.ranges)
    table.overlapping <- data.frame(table.ranges[S4Vectors::subjectHits(hits), ], stringsAsFactors = F)[, -5]
    db.genes.head.overlapping <- data.frame(db.genes.head[S4Vectors::queryHits(hits), ], stringsAsFactors = F)
    if(length(hits)>0){
        colnames(db.genes.head.overlapping)[2:6] <- c("head.start", "head.end", "overlap.width.head", "gene.strand", "gene.id")
        overlaps <- GenomicRanges::pintersect(db.genes.head[S4Vectors::queryHits(hits), ], table.ranges[S4Vectors::subjectHits(hits), ])
        width.overlap <- S4Vectors::width(overlaps)
        annotated.table.head <- cbind(table.overlapping, db.genes.head.overlapping[, c("head.start","head.end", "overlap.width.head","gene.strand","gene.id","tx.name")])

        annotated.table.head <- cbind(annotated.table.head, head.perc = width.overlap / S4Vectors::width(db.genes.head[S4Vectors::queryHits(hits), ]) * 100)
        annotated.table.head$tag <- apply(annotated.table.head,1,function(x){gsub(' +', '', paste0(x[1:3], collapse='_'))})
        annotated.table.head.ranges <- GenomicRanges::makeGRangesFromDataFrame(annotated.table.head, keep.extra.columns = TRUE)
    }

    results$heads <- annotated.table.head

    # databases overlap and database-score
    if (nrow(results$heads) > 0) {
        results$heads$dgv <- metExplorer::queryDatabase(results$heads, dgv, return.table = FALSE)
        results$heads$gnomad <- metExplorer::queryDatabase(results$heads, gnomad, return.table = FALSE)
        results$heads$NCG <- NA
        results$heads$NCG_type <- NA
        results$heads$COSMIC <- as.numeric(results$heads$gene.id %in% cosmic$ID)
        results$heads$cosmic.CGIs <- metExplorer::queryDatabase(results$heads, cosmic.CGIs, return.table = FALSE, is.genomic.ranges = TRUE,thr=thr.cgis)
        results$heads$hacer <- metExplorer::queryDatabase(results$heads, hacer, return.table = FALSE, thr = thr.enhancer)
        results$heads$CGIs <- metExplorer::queryDatabase(results$heads, cgis, return.table = FALSE, ,thr=thr.cgis)
    }


    # move "others" col to the end
    if(any(match('others', colnames(results$heads)))){
        results$heads <- cbind(results$heads[,-match('others', colnames(results$heads))],
            others = results$heads[,match('others', colnames(results$heads))]
            )
    }




    ##### EXONS
    if(annotation == 'ensembl'){
        db.exons.ranges <- readRDS(file = paste0(libpath,"/data/db.exons_ENSEMBL_",hg ,".RDS"))
    } else if (annotation=='ucsc') {
        db.exons.ranges <- readRDS(file = paste0(libpath,"/data/db.exons_UCSC_",hg ,".RDS"))
    }

    db.exons.ranges <- db.exons.ranges[txs.selected$tx.name]
    db.exons.ranges.unlist <- unlist(db.exons.ranges)

    # overlap
    annotated.table.exons <- data.frame()
    default.exons.colnames <- c("seqnames", "start", "end", "width", "beta", "gene.id", "others", "exon.start", "exon.end", "exon.width", "strand", "rank", "tx.name", "overlap.width.ex", "overlap.perc.ex", "tag")

    hits <- GenomicRanges::findOverlaps(db.exons.ranges.unlist, table.ranges)
    annotated.table.ranges.overlapping <- table.ranges[S4Vectors::subjectHits(hits), ]
    table.exons.overlapping <- data.frame(annotated.table.ranges.overlapping, stringsAsFactors = F)[, -5]
    db.exons.ranges.overlapping <- db.exons.ranges.unlist[S4Vectors::queryHits(hits), ]
    db.exons.overlapping <- data.frame(db.exons.ranges.overlapping, tx.name=names(db.exons.ranges.overlapping), stringsAsFactors = F)
    db.exons.overlapping <- db.exons.overlapping[, c("start", "end", "width", "strand", "exon_rank", "tx.name")]
    colnames(db.exons.overlapping) <- c("exon.start", "exon.end", "exon.width", "strand" ,"rank", "tx.name")

    # Overlap lengths, perc, frame
    table.exons.overlapping.tmp <- table.exons.overlapping
    overlaps <- GenomicRanges::pintersect(db.exons.ranges.overlapping, GenomicRanges::makeGRangesFromDataFrame(table.exons.overlapping.tmp, keep.extra.columns = T))
    width.overlap <- S4Vectors::width(overlaps)
    percent.overlap <- width.overlap / S4Vectors::width(db.exons.ranges.overlapping) *100
    annotated.table.exons <- cbind(table.exons.overlapping[, c('seqnames','start','end','width','beta', 'others')], db.exons.overlapping, overlap.width.ex = width.overlap, overlap.perc.ex = percent.overlap)

    # annotate
    annotations <- biomaRt::select(genedb, keys = as.character(annotated.table.exons$tx.name),columns = "GENEID", keytype = ktype.txs)
    annotated.table.exons$gene.id <- annotations$GENEID[match(annotated.table.exons$tx.name, annotations[,1])]
    annotated.table.exons$tag <- apply(annotated.table.exons,1,function(x){gsub(' +', '', paste0(x[1:3], collapse='_'))})


    # format the "other" column
    if(nrow(annotated.table.exons)>0){
        if(any(colnames(annotated.table.exons) %in% "others")){
            annotated.table.exons <- cbind(annotated.table.exons[, -match("others", colnames(annotated.table.exons))], others = annotated.table.exons[, match("others", colnames(annotated.table.exons))])
        } else {
            annotated.table.exons$others <- ""
        }
    }

    results$exons <- annotated.table.exons

    # databases overlap and database-score
    if (nrow(results$exons) > 0) {
        results$exons$dgv <- metExplorer::queryDatabase(results$exons, dgv, return.table = FALSE)
        results$exons$gnomad <- metExplorer::queryDatabase(results$exons, gnomad, return.table = FALSE)
        results$exons$NCG <- NA
        results$exons$NCG_type <- NA
        results$exons$COSMIC <- as.numeric(results$exons$gene.id %in% cosmic$ID)
        results$exons$cosmic.CGIs <- metExplorer::queryDatabase(results$exons, cosmic.CGIs, return.table = FALSE, is.genomic.ranges = TRUE, , thr = thr.cgis)
        results$exons$hacer <- metExplorer::queryDatabase(results$exons, hacer, return.table = FALSE,thr=thr.enhancer)
        results$exons$CGIs <- metExplorer::queryDatabase(results$exons, cgis, return.table = FALSE,thr=thr.cgis)
    }


    #### UTRS
    if(annotation == 'ensembl'){
        db.five.ranges <- readRDS(file = paste0(libpath,"/data/five.utrs_ENSEMBL_",hg ,".RDS"))
        db.three.ranges <- readRDS(file = paste0(libpath,"/data/three.utrs_ENSEMBL_",hg ,".RDS"))
    } else if (annotation == 'ucsc') {
        db.five.ranges <- readRDS(file = paste0(libpath,"/data/five.utrs_UCSC_",hg ,".RDS"))
        db.three.ranges <- readRDS(file = paste0(libpath,"/data/three.utrs_UCSC_",hg ,".RDS"))
    }
    db.five.ranges <- unlist(db.five.ranges[intersect(txs.selected$tx.name, names(db.five.ranges))])
    db.three.ranges <- unlist(db.three.ranges[intersect(txs.selected$tx.name, names(db.three.ranges))])

    annotated.table.fives <- data.frame()
    default.fives.colnames <- c("seqnames", "start", "end", "width", "beta", "gene.id", "others", "five.start", "five.end", "five.width", "strand", "rank", "tx.name", "overlap.width.five", "overlap.perc.five", "tag")
    
    annotated.table.threes <- data.frame()
    default.threes.colnames <- c("seqnames", "start", "end", "width", "beta", "gene.id", "others", "three.start", "three.end", "three.width", "strand", "rank", "tx.name", "overlap.width.three", "overlap.perc.three", "tag")

    # FIVE
    hits <- GenomicRanges::findOverlaps(db.five.ranges, table.ranges)
    annotated.table.ranges.overlapping <- table.ranges[S4Vectors::subjectHits(hits), ]
    if(length(S4Vectors::subjectHits(hits))>0){
        table.fives.overlapping <- data.frame(annotated.table.ranges.overlapping, stringsAsFactors = F)[, -5]
        db.fives.ranges.overlapping <- db.five.ranges[S4Vectors::queryHits(hits), ]
        db.fives.overlapping <- data.frame(db.fives.ranges.overlapping, tx.name=names(db.fives.ranges.overlapping), stringsAsFactors = F)
        db.fives.overlapping <- db.fives.overlapping[, c("start", "end", "width", "strand", "exon_rank", "tx.name")]
        colnames(db.fives.overlapping) <- c("five.start", "five.end", "five.width", "strand" ,"rank", "tx.name")
    
        # Overlap lengths, perc, frame
        table.fives.overlapping.tmp <- table.fives.overlapping

        # table.fives.overlapping.tmp[, 3] <- table.fives.overlapping.tmp[, 2] + table.fives.overlapping.tmp[, 4]
        overlaps <- GenomicRanges::pintersect(db.fives.ranges.overlapping, GenomicRanges::makeGRangesFromDataFrame(table.fives.overlapping.tmp, keep.extra.columns = T))
        width.overlap <- S4Vectors::width(overlaps)
        percent.overlap <- width.overlap / S4Vectors::width(db.fives.ranges.overlapping) *100
        frame <- '-'
        annotated.table.fives <- cbind(table.fives.overlapping[, c('seqnames','start','end','width','beta',  'others')], db.fives.overlapping, overlap.width.five = width.overlap, overlap.perc.five = percent.overlap, frame = frame)
        
        # annotate
        annotations <- biomaRt::select(genedb, keys = as.character(annotated.table.fives$tx.name),columns = "GENEID", keytype = ktype.txs)
        annotated.table.fives$gene.id <- annotations$GENEID[match(annotated.table.fives$tx.name, annotations[,1])]
        annotated.table.fives$tag = apply(annotated.table.fives,1, function(x) {gsub(' +', '', paste0(x[1:3], collapse='_'))})
    }



    # THREE
    hits <- GenomicRanges::findOverlaps(db.three.ranges, table.ranges)
    annotated.table.ranges.overlapping <- table.ranges[S4Vectors::subjectHits(hits), ]
    if(length(S4Vectors::subjectHits(hits))>0){
        table.threes.overlapping <- data.frame(annotated.table.ranges.overlapping, stringsAsFactors = F)[, -5]
        db.threes.ranges.overlapping <- db.three.ranges[S4Vectors::queryHits(hits), ]
        db.threes.overlapping <- data.frame(db.threes.ranges.overlapping, tx.name=names(db.threes.ranges.overlapping), stringsAsFactors = F)
        db.threes.overlapping <- db.threes.overlapping[, c("start", "end", "width", "strand", "exon_rank", "tx.name")]
        colnames(db.threes.overlapping) <- c("three.start", "three.end", "three.width", "strand" ,"rank", "tx.name")
    
        # Overlap lengths, perc, frame
        table.threes.overlapping.tmp <- table.threes.overlapping

        overlaps <- GenomicRanges::pintersect(db.threes.ranges.overlapping, GenomicRanges::makeGRangesFromDataFrame(table.threes.overlapping.tmp, keep.extra.columns = T))
        width.overlap <- S4Vectors::width(overlaps)
        percent.overlap <- width.overlap / S4Vectors::width(db.threes.ranges.overlapping) *100
        frame <- '-'
        annotated.table.threes <- cbind(table.threes.overlapping[, c('seqnames','start','end','width','beta',  'others')], db.threes.overlapping, overlap.width.three = width.overlap, overlap.perc.three = percent.overlap, frame = frame)
        
        # annotate
        annotations <- biomaRt::select(genedb, keys = as.character(annotated.table.threes$tx.name),columns = "GENEID", keytype = ktype.txs)
        annotated.table.threes$gene.id <- annotations$GENEID[match(annotated.table.threes$tx.name, annotations[,1])]
        annotated.table.threes$tag = apply(annotated.table.threes,1, function(x) {gsub(' +', '', paste0(x[1:3], collapse='_'))})
    }


    # distances from 'intron-five' junction
    if(nrow(annotated.table.fives)>0){

        distances <- cbind(annotated.table.fives["start"] - annotated.table.fives["five.start"], annotated.table.fives["five.end"] - annotated.table.fives["end"])
        colnames(distances) <- c("distance.left", "distance.right")
        annotated.table.fives <- cbind(annotated.table.fives, distances)

        if(any(colnames(annotated.table.fives) %in% "others")){
            annotated.table.fives <- cbind(annotated.table.fives[, -match("others", colnames(annotated.table.fives))], others = annotated.table.fives[, match("others", colnames(annotated.table.fives))])
        } else {
            annotated.table.fives$others <- ""
        }
    }

    results$fiveUTRs <- annotated.table.fives

    # distances from 'intron-three' junction
    if(nrow(annotated.table.threes)>0){

        distances <- cbind(annotated.table.threes["start"] - annotated.table.threes["three.start"], annotated.table.threes["three.end"] - annotated.table.threes["end"])
        colnames(distances) <- c("distance.left", "distance.right")
        annotated.table.threes <- cbind(annotated.table.threes, distances)

        if(any(colnames(annotated.table.threes) %in% "others")){
            annotated.table.threes <- cbind(annotated.table.threes[, -match("others", colnames(annotated.table.threes))], others = annotated.table.threes[, match("others", colnames(annotated.table.threes))])
        } else {
            annotated.table.threes$others <- ""
        }
    }

    results$threeUTRs <- annotated.table.threes


    # databases overlap and database-score
    if (nrow(results$fiveUTRs) > 0) {
        results$fiveUTRs$dgv <- metExplorer::queryDatabase(results$fiveUTRs, dgv, return.table = FALSE)
        results$fiveUTRs$gnomad <- metExplorer::queryDatabase(results$fiveUTRs, gnomad, return.table = FALSE)
        results$fiveUTRs$NCG <- NA
        results$fiveUTRs$NCG_type <- NA
        results$fiveUTRs$COSMIC <- as.numeric(results$fiveUTRs$gene.id %in% cosmic$ID)
        results$fiveUTRs$cosmic.CGIs <- metExplorer::queryDatabase(results$fiveUTRs, cosmic.CGIs, return.table = FALSE, is.genomic.ranges = TRUE, , thr = thr.cgis)
        results$fiveUTRs$hacer <- metExplorer::queryDatabase(results$fiveUTRs, hacer, return.table = FALSE, thr = thr.enhancer)
        results$fiveUTRs$CGIs <- metExplorer::queryDatabase(results$fiveUTRs, cgis, return.table = FALSE, thr = thr.cgis)
    }

    if (nrow(results$threeUTRs) > 0) {
        results$threeUTRs$dgv <- metExplorer::queryDatabase(results$threeUTRs, dgv, return.table = FALSE)
        results$threeUTRs$gnomad <- metExplorer::queryDatabase(results$threeUTRs, gnomad, return.table = FALSE)
        results$threeUTRs$NCG <- NA
        results$threeUTRs$NCG_type <- NA
        results$threeUTRs$COSMIC <- as.numeric(results$threeUTRs$gene.id %in% cosmic$ID)
        results$threeUTRs$cosmic.CGIs <- metExplorer::queryDatabase(results$threeUTRs, cosmic.CGIs, return.table = FALSE, is.genomic.ranges = TRUE, thr = thr.cgis)
        results$threeUTRs$hacer <- metExplorer::queryDatabase(results$threeUTRs, hacer, return.table = FALSE, thr = thr.enhancer)
        results$threeUTRs$CGIs <- metExplorer::queryDatabase(results$threeUTRs, cgis, return.table = FALSE, thr = thr.cgis)
    }


    ##### PROMOTERS
    # promoter coordinates
    proms.coordinates <- genes.coordinates
    proms.coordinates[which(strand.map), "end"] <- genes.coordinates[which(strand.map), "start"] + 1
    proms.coordinates[which(strand.map), "start"] <- genes.coordinates[which(strand.map), "start"] - prom.length
    proms.coordinates[-which(strand.map), "start"] <- genes.coordinates[-which(strand.map), "end"] - 1
    proms.coordinates[-which(strand.map), "end"] <- genes.coordinates[-which(strand.map), "end"] + prom.length
    db.proms.ranges <- GenomicRanges::makeGRangesFromDataFrame(proms.coordinates, keep.extra.columns = T)

    # overlap to find genes
    annotated.table.proms <- data.frame()
    annotated.proms.table.colnames <- c("seqnames", "start", "end", "width", "strand", "beta", "others", "gene.id", "prom.start", "prom.end", "prom.width", "tx.name", "overlap.width.prom", "overlap.perc.prom", "tag")

    hits <- GenomicRanges::findOverlaps(db.proms.ranges, table.ranges)
    prom.overlapping.ranges <- db.proms.ranges[S4Vectors::queryHits(hits), ]
    table.prom.overlapping <- table.ranges[S4Vectors::subjectHits(hits), ]
    table.prom.overlapping <- data.frame(table.prom.overlapping)
    table.prom.overlapping$tag <- apply(table.prom.overlapping, 1, function(x) {
        gsub(" +", "", paste0(x[1:3], collapse = "_"))
    })
    table.prom.overlapping.tmp <- table.prom.overlapping
    overlaps <- GenomicRanges::pintersect(prom.overlapping.ranges, GenomicRanges::makeGRangesFromDataFrame(table.prom.overlapping.tmp, keep.extra.columns = T))
    width.overlap <- S4Vectors::width(overlaps)
    percent.overlap <- width.overlap / S4Vectors::width(prom.overlapping.ranges) * 100
    prom.overlapping.df <- data.frame(prom.overlapping.ranges)
    colnames(prom.overlapping.df) <- c("seqnames", "prom.start", "prom.end", "prom.width", "strand", "gene.id", "tx.name")
    annotated.table.proms <- cbind(table.prom.overlapping[,c("seqnames", "start", "end", "width", "beta", "others", "tag")], prom.overlapping.df[, c("prom.start", "prom.end","prom.width","strand", "gene.id", "tx.name")],overlap.width.prom = width.overlap, overlap.perc.prom = percent.overlap)
    annotated.table.proms <- annotated.table.proms[, annotated.proms.table.colnames]

    # # distances from junction

    if(length(annotated.table.proms$others)>0){
        results$promoters <- cbind(annotated.table.proms[,-7], others = annotated.table.proms[, "others"])
    } else {
        results$promoters <- annotated.table.proms
    }

    # databases overlap and database-score
    if (nrow(results$promoters) > 0) {
        results$promoters$dgv <- metExplorer::queryDatabase(results$promoters, dgv, return.table = FALSE)
        results$promoters$gnomad <- metExplorer::queryDatabase(results$promoters, gnomad, return.table = FALSE)
        results$promoters$NCG <- NA
        results$promoters$NCG_type <- NA
        results$promoters$COSMIC <- as.numeric(results$promoters$gene.id %in% cosmic$ID)
        results$promoters$cosmic.CGIs <- metExplorer::queryDatabase(results$promoters, cosmic.CGIs, return.table = FALSE, is.genomic.ranges = TRUE, thr = thr.cgis)
        results$promoters$hacer <- metExplorer::queryDatabase(results$promoters, hacer, return.table = FALSE,thr=thr.enhancer)
        results$promoters$CGIs <- metExplorer::queryDatabase(results$promoters, cgis, return.table = FALSE, thr = thr.cgis)
    }


    ###TSS SURROUNDING
    coordinates.tss <- proms.coordinates
    cur.pos <- coordinates.tss$strand == "+"
    cur.neg <- coordinates.tss$strand == "-"
    coordinates.tss[cur.pos, "end"] = coordinates.tss[cur.pos, "end"] + head.length
    coordinates.tss[cur.neg, "start"] <- coordinates.tss[cur.neg, "start"] - head.length

    db.tss.ranges <- GenomicRanges::makeGRangesFromDataFrame(coordinates.tss, keep.extra.columns = T)

    # overlap
    annotated.table.tss <- data.frame()
    annotated.tss.table.colnames <- c("seqnames", "start", "end", "width", "strand", "beta", "others", "gene.id", "tss.sur.start", "tss.sur.end", "tss.sur.width", "tx.name", "overlap.width.tss.sur", "overlap.perc.tss.sur", "tag")

    hits <- GenomicRanges::findOverlaps(db.tss.ranges, table.ranges)
    tss.overlapping.ranges <- db.tss.ranges[S4Vectors::queryHits(hits), ]
    table.tss.overlapping <- table.ranges[S4Vectors::subjectHits(hits), ]
    table.tss.overlapping <- data.frame(table.tss.overlapping)
    table.tss.overlapping$tag <- apply(table.tss.overlapping, 1, function(x) {
        gsub(" +", "", paste0(x[1:3], collapse = "_"))
    })
    table.tss.overlapping.tmp <- table.tss.overlapping
    overlaps <- GenomicRanges::pintersect(tss.overlapping.ranges, GenomicRanges::makeGRangesFromDataFrame(table.tss.overlapping.tmp, keep.extra.columns = T))
    width.overlap <- S4Vectors::width(overlaps)
    percent.overlap <- width.overlap / S4Vectors::width(tss.overlapping.ranges) * 100
    tss.overlapping.df <- data.frame(tss.overlapping.ranges)
    colnames(tss.overlapping.df) <- c("seqnames", "tss.sur.start", "tss.sur.end", "tss.sur.width", "strand", "gene.id", "tx.name")
    annotated.table.tss <- cbind(table.tss.overlapping[, c("seqnames", "start", "end", "width", "beta", "others", "tag")], tss.overlapping.df[, c("tss.sur.start", "tss.sur.end", "tss.sur.width", "strand", "gene.id", "tx.name")], overlap.width.tss.sur = width.overlap, overlap.perc.tss.sur = percent.overlap)
    annotated.table.tss <- annotated.table.tss[, annotated.tss.table.colnames]
    
    if(any(colnames(annotated.table.tss) %in% "others")){
        annotated.table.tss <- cbind(annotated.table.tss[, -match("others", colnames(annotated.table.tss))], others = annotated.table.tss[, match("others", colnames(annotated.table.tss))])
    } else {
        annotated.table.tss$others <- ""
    }


    results$tss.surrounding <- annotated.table.tss

    # databases overlap and database-score
    if (nrow(results$tss.surrounding) > 0) {
        results$tss.surrounding$dgv <- metExplorer::queryDatabase(results$tss.surrounding, dgv, return.table = FALSE)
        results$tss.surrounding$gnomad <- metExplorer::queryDatabase(results$tss.surrounding, gnomad, return.table = FALSE)
        results$tss.surrounding$NCG <- NA
        results$tss.surrounding$NCG_type <- NA
        results$tss.surrounding$COSMIC <- as.numeric(results$tss.surrounding$gene.id %in% cosmic$ID)
        results$tss.surrounding$cosmic.CGIs <- metExplorer::queryDatabase(results$tss.surrounding, cosmic.CGIs, return.table = FALSE, is.genomic.ranges = TRUE, thr = thr.cgis)
        results$tss.surrounding$hacer <- metExplorer::queryDatabase(results$tss.surrounding, hacer, return.table = FALSE, thr = thr.enhancer)
        results$tss.surrounding$CGIs <- metExplorer::queryDatabase(results$tss.surrounding, cgis, return.table = FALSE, thr = thr.cgis)
    }
  




    ##### INTRONS
    if(annotation == 'ensembl'){
        db.introns.ranked.gr <- readRDS(paste0(libpath,"/data/intronsRankedAll_ENSEMBL_",hg,".RDS"))
    } else if (annotation == 'ucsc'){
        db.introns.ranked.gr <- readRDS(paste0(libpath,"/data/intronsRankedAll_UCSC_",hg,".RDS"))
    }

    names(db.introns.ranked.gr) <- data.frame(db.introns.ranked.gr, stringsAsFactors = F)$group_name
    db.introns.ranked.gr <- db.introns.ranked.gr[intersect(names(db.introns.ranked.gr), txs.selected$tx.name)]

    # overlap introns with SVs
    annotated.table.introns <- data.frame()
    annotated.table.introns.colnames <- c("seqnames", "start", "end", "width", "beta", "gene.id", "others", "intron.start", "intron.end", "intron.width", "intron.rank", "tx.name", "overlap.width.intron", "overlap.perc.intron", "tag")

    hits <- GenomicRanges::findOverlaps(db.introns.ranked.gr, table.ranges)
    annotated.table.ranges.overlapping <- table.ranges[S4Vectors::subjectHits(hits), ]
    table.introns.overlapping <- data.frame(annotated.table.ranges.overlapping, stringsAsFactors = F)[, -5]

    db.introns.ranges.overlapping <- db.introns.ranked.gr[S4Vectors::queryHits(hits), ]
    db.introns.overlapping <- data.frame(db.introns.ranges.overlapping, stringsAsFactors = F)

    colnames(db.introns.overlapping) <- c("seqnames", "intron.start", "intron.end", "intron.width", "strand", "group", "tx.name", "intron.rank")
    annotated.table.introns <- cbind(table.introns.overlapping, db.introns.overlapping[, -1])

    # SCORING OVERLAPS
    overlaps <- GenomicRanges::pintersect(db.introns.ranges.overlapping, GenomicRanges::makeGRangesFromDataFrame(annotated.table.introns, keep.extra.columns = T))

    width.overlap <- S4Vectors::width(overlaps)
    percent.overlap <- width.overlap / S4Vectors::width(db.introns.ranges.overlapping) * 100
    annotated.table.introns$overlap.width.intron <- width.overlap
    annotated.table.introns$overlap.perc.intron <- percent.overlap
    annotated.table.introns$tag <- apply(annotated.table.introns, 1, function(x) {
        gsub(" +", "", paste0(x[1:3], collapse = "_"))
    })

    # annotate
    annotations <- biomaRt::select(genedb, keys = as.character(annotated.table.introns$tx.name), columns = "GENEID", keytype = ktype.txs)
    annotated.table.introns$gene.id <- annotations$GENEID[match(annotated.table.introns$tx.name, annotations[,1])]
    annotated.table.introns <- annotated.table.introns[, annotated.table.introns.colnames]


    # distances from 'intron-exon' junction
    distances <- cbind(distance.left = annotated.table.introns["start"] - annotated.table.introns["intron.start"], distance.right = annotated.table.introns["intron.end"] - annotated.table.introns["end"])
    intron.annotations <- cbind(annotated.table.introns, distance.left = distances[, 1], distance.right = distances[, 2])

    if(nrow(intron.annotations)<0){
        if(any(colnames(intron.annotations) %in% "others")){
            intron.annotations <- cbind(intron.annotations[, -match("others", colnames(intron.annotations))], others = intron.annotations$others)
        }
    }

    results$introns  <- intron.annotations

    # databases overlap and database-score
    if (nrow(results$introns) > 0) {
        results$introns$dgv <- metExplorer::queryDatabase(results$introns, dgv, return.table = FALSE)
        results$introns$gnomad <- metExplorer::queryDatabase(results$introns, gnomad, return.table = FALSE)
        results$introns$NCG <- NA
        results$introns$NCG_type <- NA
        results$introns$COSMIC <- as.numeric(results$introns$gene.id %in% cosmic$ID)
        results$introns$cosmic.CGIs <- metExplorer::queryDatabase(results$introns, cosmic.CGIs, return.table = FALSE, is.genomic.ranges = TRUE, thr = thr.cgis)
        results$introns$hacer <- metExplorer::queryDatabase(results$introns, hacer, return.table = FALSE, thr = thr.enhancer)
        results$introns$CGIs <- metExplorer::queryDatabase(results$introns, cgis, return.table = FALSE, thr = thr.cgis)
    }


    # converting to gene symbols

    if(annotation.fast == TRUE){
        for(i in 1:length(results)){
            if(nrow(results[[i]])>0){
                annotations <- biomaRt::select(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75, keys = as.character(results[[i]]$gene.id), columns = "SYMBOL", keytype = ktype.genes)
                results[[i]]$symbol <- annotations$SYMBOL[match(results[[i]]$gene.id, annotations[,match(ktype.genes, colnames(annotations))])]
            }
        }
    } else {
        for(i in 1:length(results)){
            to.trans <- as.character(results[[i]]$gene.id)
            if(length(to.trans)>0){
                gene.symbols <- biomaRt::select(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75, keys = to.trans, columns = "SYMBOL", keytype = ktype.genes)
                gene.symbols.tags <- apply(gene.symbols,1,function(x) {paste0(x,collapse='')})
                gene.symbols <- gene.symbols[!duplicated(gene.symbols.tags),]
                
                # collapse duplicated annotations
                cur.ind.dup <- which(duplicated(gene.symbols[,1]))
                if(length(cur.ind.dup)>0){
                    gene.symbols.duplicated <- gene.symbols[cur.ind.dup,]
                    for(k in 1:nrow(unique(gene.symbols.duplicated))){
                        ind.gene.symbols.duplicated = which(gene.symbols[,1] %in% unique(gene.symbols.duplicated[,1])[k])
                        gene.symbols[ind.gene.symbols.duplicated,2] <- gsub(' +', '', (paste0(gene.symbols[ind.gene.symbols.duplicated,2],collapse=";")))
                    }
                    gene.symbols <- gene.symbols[-cur.ind.dup,]
                }


                # add symbols column
                results[[i]]$symbol <- ""
                for(j in 1:nrow(gene.symbols)){
                    ind.results = which(results[[i]]$gene.id %in% gene.symbols[j,1])
                    results[[i]]$symbol[ind.results] = gene.symbols[j,2]
                }
            }
        }
    }

    # annotate NCG
    results <- genesToNCG(results, ncg, return.table = FALSE)
    results <- genesToNCG(results, ncg, return.table = TRUE)

    # find TFs
    tf.genes <- readRDS(paste0(libpath, "/data/ENCODE_targets.RDS"))

    results$promoters$TF = 0
    results$introns$TF = 0
    results$threeUTRs$TF = 0
    results$fiveUTRs$TF = 0
    results$exons$TF = 0
    results$genes$TF = 0
    results$heads$TF = 0
    results$tss.surrounding$TF = 0

    results$promoters$TF[results$promoters$symbol %in% tf.genes[,1] %in% tf.genes[,1]] = 1
    results$introns$TF[results$introns$symbol %in% tf.genes[,1]] = 1
    results$threeUTRs$TF[results$threeUTRs$symbol %in% tf.genes[,1]] = 1
    results$fiveUTRs$TF[results$fiveUTRs$symbol %in% tf.genes[,1]] = 1
    results$exons$TF[results$exons$symbol %in% tf.genes[,1]] = 1
    results$genes$TF[results$genes$symbol %in% tf.genes[,1]] = 1
    results$heads$TF[results$heads$symbol %in% tf.genes[,1]] = 1
    results$tss.surrounding$TF[results$tss.surrounding$symbol %in% tf.genes[,1]] = 1



    # computing database score
    results$promoters$database.score <- apply(cbind(results$promoters$dgv, results$promoters$gnomad, results$promoters$NCG, results$promoters$COSMIC, results$promoters$cosmic.CGIs, results$promoters$TF), 1, sum)
    results$introns$database.score <- apply(cbind(results$introns$dgv, results$introns$gnomad, results$introns$NCG, results$introns$COSMIC, results$introns$cosmic.CGIs, results$introns$TF), 1, sum)
    results$threeUTRs$database.score <- apply(cbind(results$threeUTRs$dgv, results$threeUTRs$gnomad, results$threeUTRs$NCG, results$threeUTRs$COSMIC, results$threeUTRs$cosmic.CGIs, results$threeUTRs$TF), 1, sum)
    results$fiveUTRs$database.score <- apply(cbind(results$fiveUTRs$dgv, results$fiveUTRs$gnomad, results$fiveUTRs$NCG, results$fiveUTRs$COSMIC, results$fiveUTRs$cosmic.CGIs, results$fiveUTRs$TF), 1, sum)
    results$exons$database.score <- apply(cbind(results$exons$dgv, results$exons$gnomad, results$exons$NCG, results$exons$COSMIC, results$exons$cosmic.CGIs, results$exons$TF), 1, sum)
    results$genes$database.score <- apply(cbind(results$genes$dgv, results$genes$gnomad, results$genes$NCG, results$genes$COSMIC, results$genes$cosmic.CGIs, results$genes$TF), 1, sum)
    results$heads$database.score <- apply(cbind(results$heads$dgv, results$heads$gnomad, results$heads$NCG, results$heads$COSMIC, results$heads$cosmic.CGIs, results$heads$TF), 1, sum)
    results$tss.surrounding$database.score <- apply(cbind(results$tss.surrounding$dgv, results$tss.surrounding$gnomad, results$tss.surrounding$NCG, results$tss.surrounding$COSMIC, results$tss.surrounding$cosmic.CGIs, results$tss.surrounding$TF), 1, sum)

    # reordering output
    for(p in 1:length(results)){
        ind.tf <- match('TF', colnames(results[[p]]))
        ind.others <- match('others', colnames(results[[p]]))
        ind.dbs <- match('database.score', colnames(results[[p]]))
        ind.symbol <- match('symbol', colnames(results[[p]]))
        results[[p]] <- results[[p]][,c((1:ncol(results[[p]]))[c(-ind.tf, -ind.symbol, -ind.dbs, -ind.others)], ind.tf, ind.symbol, ind.dbs, ind.others)]
    }

    # sort the output
    order.output <- c("genes", "heads", "tss.surrounding", "promoters", "fiveUTRs", "exons", "introns", "threeUTRs")
    results <- results[order.output]
    return(results)

}



