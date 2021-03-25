#' Query different databases to find enriched proceses. 
#' 
#' Query TF's targeted genes in order to find enriched proceses.
#'  
#' @param associatedTFs2Expr data.frame. Corresponding to resulting data.frame from associateTFs2Exprs().
#' @param logfc.thr numeric value indicating logFC threshold. Default = 1.
#' @param stat.filter character indicating which type of statistics use for filtering results. Accepted values: 'P.value' or 'Adjusted.P.value'. Default 'P.value'.
#' @param stat.thr numeric value indicating the threshold to use for selcted statistical test. Default 0.01.
#' @param db vector of characters indicating DBs to query. Default c("ClinVar_2019", "OMIM_Disease", "Elsevier_Pathway_Collection", "MSigDB_Hallmark_2020", "MSigDB_Oncogenic_Signatures", "GO_Biological_Process_2018", "Human_Phenotype_Ontology", "KEGG_2016", "NCI-Nature_2016", "Panther_2016", "Reactome_2016", "WikiPathways_2019_Human")
#' @return data.frame with enriched processes
#' 
#' @export

tfs2Enrichr <- function(associatedTFs2Expr, logfc.thr = 1, stat.filter = 'P.value', stat.thr = 0.01, db = NULL){
    if(!any((.packages()) %in% "enrichR")){
        library(enrichR)
    }

    associatedTFs2Expr <- associatedTFs2Expr[abs(associatedTFs2Expr$target.expression)>=logfc.thr,]

    if(is.null(db)){
        manifest <- c("ClinVar_2019", "OMIM_Disease", "Elsevier_Pathway_Collection", "MSigDB_Hallmark_2020", "MSigDB_Oncogenic_Signatures", "GO_Biological_Process_2018", "Human_Phenotype_Ontology", "KEGG_2016", "NCI-Nature_2016", "Panther_2016", "Reactome_2016", "WikiPathways_2019_Human")
    } else {
        manifest <- db
    }
   
    enirchr.results <- enrichR::enrichr(as.character(associatedTFs2Expr$target.symbol), databases = manifest)

    # stat filtering
    enirchr.results.parsed <- lapply(
        enirchr.results, function(x) {
            if(nrow(x)>0){
                x[x[stat.filter] < stat.thr, ]
            }
        }
    )

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
