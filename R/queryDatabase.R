#' Query database to find pathologic genes. 
#' 
#' Query elements of resulting annotations list to specified database in order to find pathogenic genes.
#'  
#' @param table element of results list from annotatedDMRs()
#' @param database the database to query
#' @param return.table optition to return a table instead a vector of presence
#' @param hold.columns column positions to hold when return.table=TRUE
#' @param is.genomic.ranges specifies if database is a GenomicRange objet or a data.frame
#' @param thr beta value threshold   
#' 
#' @return a vector of presence or a data.frame 
#' 
#' @export

queryDatabase <- function(table, database, return.table = TRUE, hold.columns, is.genomic.ranges = FALSE, thr=0) {

    # make Genomic Ranges
    gene.annotations.ranges <- GenomicRanges::makeGRangesFromDataFrame(table, keep.extra.columns = TRUE)

    if(is.genomic.ranges==FALSE){
        database.ranges <- GenomicRanges::makeGRangesFromDataFrame(database, keep.extra.columns = TRUE)
    } else if(is.genomic.ranges==TRUE) {
        database.ranges <- database
    }

    # overlaps features/DB
    hits <- GenomicRanges::findOverlaps(gene.annotations.ranges, database.ranges)
    gene.annotations.ranges.parsed <- gene.annotations.ranges[S4Vectors::queryHits(hits)]
    database.ranges.parsed <- database.ranges[S4Vectors::subjectHits(hits)]

    if (return.table == TRUE) {
        if(!missing(hold.columns)){
            table.overlapping.database.annotated <- cbind(data.frame(gene.annotations.ranges.parsed), data.frame(database.ranges.parsed)[,hold.columns])
            if(length(hold.columns)>1){
                col.n <- ncol(table.overlapping.database.annotated) - length(hold.columns) + 1
                colnames(table.overlapping.database.annotated)[col.n:ncol(table.overlapping.database.annotated)] <- hold.columns
            } else {
                col.n <- ncol(table.overlapping.database.annotated)
                colnames(table.overlapping.database.annotated)[col.n] <- hold.columns
            }
        } else if(missing(hold.columns)) {
            table.overlapping.database.annotated <- cbind(data.frame(gene.annotations.ranges.parsed), data.frame(database.ranges.parsed))
        }

        intrsct <- GenomicRanges::pintersect(gene.annotations.ranges.parsed, database.ranges.parsed)
        overlap.lengths <- S4Vectors::width(intrsct) / S4Vectors::width(database.ranges.parsed)
        ind.hold <- (overlap.lengths) >= thr
        if(any(ind.hold)){
            table.overlapping.database.annotated$overlap.width <- overlap.lengths[ind.hold]*100
        } else {
            stop('no region satisfy threshold value')
        }
        return(table.overlapping.database.annotated[ind.hold,])
    } else {
        intrsct <- GenomicRanges::pintersect(gene.annotations.ranges.parsed, database.ranges.parsed)
        overlap.lengths <- S4Vectors::width(intrsct) / S4Vectors::width(database.ranges.parsed)
        ind.hold <- (overlap.lengths) >= thr
        vector.presence <- rep(0, nrow(table))
        if(any(ind.hold)){
            vector.presence[S4Vectors::queryHits(hits)[ind.hold]] <- 1
        }
        return(vector.presence)
    }
}
