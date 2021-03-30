#' Query different databases to find enriched proceses. 
#' 
#' Query gene symbols from resulting annotations list in order to find enriched proceses.
#'  
#' @param annotatedDMRs anotated DMRs list resulting from annotateDMRs() or scoreAnnotatedDMRs()
#' @param active.features annotation level from the gene symbols are taken. Default = c('promoters', 'heads') 
#' @param stat.filter character indicating which type of statistics use for filtering results. Accepted values: 'P.value', 'Adjusted.P.value' or 'Overlap'. Default 'P.value'.
#' @param stat.thr numeric value indicating the threshold to use for selcted statistical test. Default 0.01.
#' @param db vector of characters indicating DBs to query. Default c("ClinVar_2019", "OMIM_Disease", "Elsevier_Pathway_Collection", "MSigDB_Hallmark_2020", "MSigDB_Oncogenic_Signatures", "GO_Biological_Process_2018", "Human_Phenotype_Ontology", "KEGG_2016", "NCI-Nature_2016", "Panther_2016", "Reactome_2016", "WikiPathways_2019_Human")
#' 
#' @return data.frame with enriched processes
#' 
#' @export

annotatedDMRs2Enrichr <- function(annotatedDMRs, active.features = c('promoters', 'heads'), stat.filter = 'P.value', stat.thr = 0.01, db = NULL){
    if (!any((.packages()) %in% "enrichR")){
        library(enrichR)
    }

    if(stat.filter=="Overlap"){
        stat.thr <- min(stat.thr,1)
    }

    if(is.null(db)){
        manifest <- c("ClinVar_2019", "OMIM_Disease", "Elsevier_Pathway_Collection", "MSigDB_Hallmark_2020", "MSigDB_Oncogenic_Signatures", "GO_Biological_Process_2018", "Human_Phenotype_Ontology", "KEGG_2016", "NCI-Nature_2016", "Panther_2016", "Reactome_2016", "WikiPathways_2019_Human")
    } else {
        manifest <- db
    }

    query.vector <- c()
    for(af in active.features){
        query.vector <- c(query.vector, annotatedDMRs[[af]]$symbol)
    } 

    enirchr.results <- enrichR::enrichr(query.vector, databases = manifest)

    # filtering
    if(stat.filter!="Overlap"){
        enirchr.results.parsed <- lapply(
                enirchr.results, function(x) {
                    x[x[stat.filter] <= stat.thr, ]
                }
        )
    } else {
        enirchr.results.parsed <- lapply(
            enirchr.results, function(x) {
                x[x[stat.filter] >= stat.thr, ]
            }
        )
    }

    # filter out not enriched db
    ind <- which(
        unlist(
            lapply(
                enirchr.results.parsed, function(x) {
                    nrow(x) != 0
                }
            )
        )
    )

    enirchr.results.table.complete <- do.call("rbind", enirchr.results.parsed[ind])
    enirchr.results.table.complete$Db <- gsub('\\.\\d+','',rownames(enirchr.results.table.complete))
    enirchr.results.table.complete <- enirchr.results.table.complete[order(enirchr.results.table.complete[,stat.filter]),]
    return(enirchr.results.table.complete[,c("Db","Term", "Overlap", "P.value", "Adjusted.P.value","Genes")])
}
