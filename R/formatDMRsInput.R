#' Convert input table to proper format. 
#' 
#' Convert input table to proper format. The firts three column of the input table must have chr, start, end coordinates. 
#'  
#' @param DMRsRanges the DMRs ranges, it must have the following columns: chr, start, end, beta diff. Other columns will be stored in the resulting output under the column other.
#' @param thr.beta numeric, beta difference threshold to consider methylations. Default = 0.3
#' @param col.betadiff numeric, column position for beta diff. in input table. Default = 4
#' @param col.beta1 numeric, column position for first sample beta values in input table
#' @param col.beta2 numeric, column position for second sample beta values in input table
#' 
#' @return data.frame of DMR ranges
#' 
#' @export


formatDMRsInput <- function(tableIn, thr.beta, col.betadiff = 4, col.beta1 = NULL, col.beta2 = NULL){

    tableIn[, 2] <- as.numeric(as.character(tableIn[, 2]))
    tableIn[, 3] <- as.numeric(as.character(tableIn[, 3]))
    tableIn[, col.betadiff] <- as.numeric(as.character(tableIn[, col.betadiff]))
    tableIn <- tableIn[abs(tableIn[, col.betadiff]) >= thr.beta, ]

    if(!is.null(col.beta1)){
        colnames(tableIn)[col.beta1] <- "beta1"
    }

    if (!is.null(col.beta2)) {
        colnames(tableIn)[col.beta2] <- "beta2"
    }


    indChr <- grep("chr", tableIn[, 1])
    if(length(indChr) == 0){
        tableIn[, 1]= paste("chr", tableIn[, 1], sep = "")
    }
    tableIn <- tableIn[grep("^chr\\d{1,2}$", tableIn[, 1]), ]
    
    # format the tableIn in order to get beta values after the 3rd column
    indBeta = grep('^[b,B]eta$',colnames(tableIn))
    if(length(indBeta)==0){
        indBeta = 4
    }
    tableInCore = cbind(tableIn[,1:3], beta=tableIn[,indBeta])
    if(ncol(tableInCore)!=ncol(tableIn)){
            tableInExtra = tableIn[,-c(indBeta, 1:3)]
            tableIn = cbind(tableInCore, tableInExtra)
    } else {
        tableIn = tableInCore
    }

    # formatting extra annotations
    if (ncol(tableIn[-c(1:3, indBeta)]) > 1) {
        tableIn[, 5] <- apply(tableIn, 1, function(x) {
            gsub(' *','', paste(names(x)[5:length(x)],' = ',x[5:length(x)], collapse = ";", sep = ""))
        })

        # adding standard colnames to the tableIn
        tableIn <- tableIn[,c(1:5)]
        colnames(tableIn)[5] <- "others"
    } else if (ncol(tableIn[-c(1:3, indBeta)]) == 1) {
        tableIn[,5] <- paste0(colnames(tableIn)[5], ' = ', tableIn[,5])
        colnames(tableIn)[5] <- "others"
    } else {
        tableIn$others = NA
    }
    colnames(tableIn)[1:4] <- c("seqname", "start", "end", "beta")
   return(tableIn)
}
