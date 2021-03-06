% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/associateTFs2Exprs.R
\name{associateTFs2Exprs}
\alias{associateTFs2Exprs}
\title{Associate target genes to TFs and retrieve their expression}
\usage{
associateTFs2Exprs(
  annotatedDMRs,
  expressionProfile,
  active.features = c("promoters", "heads"),
  col.genes = 0,
  col.stat = 6,
  stat.thr = 0.05,
  col.logFC = 2,
  logfc.thr = 0,
  convert.genes = FALSE,
  convert.from,
  beta.thr = 0.3,
  overlap.param.thr = 30,
  param.type = "overlap.percentage"
)
}
\arguments{
\item{annotatedDMRs}{anotated DMRs list resultingfrom annotateDMRs() or scoreAnnotatedDMRs()}

\item{expressionProfile}{expression data.frame}

\item{active.features}{character vectors containing features to correlate. Must be from names of resulting list from annotateDMRs. Additional feature names can be first exons (exons1) or first intron (intron1). To use more than one feature use c(). Default = c("promoters", "heads")}

\item{col.genes}{numeric, the column of expressionProfile data.frame with gene Ids. If NULL geneIDs will be taken from rownames() of expressionProfile. Default = 0.}

\item{col.stat}{numeric, the column of expressionProfile data.frame with the statistics to use. Default = 6.}

\item{stat.thr}{numeric, threshold for statistical significance. Default = 0.05}

\item{col.logFC}{numeric, the column of expressionProfile data.frame with log. fold change. Default = 2}

\item{logfc.thr}{numeric, threshold value for log. fold change. Default = 0.}

\item{convert.genes}{logical, used to indicate if gene ids have to be translated in official gene symbols. Default = FALSE}

\item{convert.from}{character, used annotation for gene in expressionProfile to be converted to symbols gene IDs. Accepted: c("ENTREZID" ,"EXONID" ,"GENEBIOTYPE" ,"GENEID" ,"GENENAME" ,"PROTDOMID" ,"PROTEINDOMAINID" ,"PROTEINDOMAINSOURCE" ,"PROTEINID" ,""SEQSTRAND" ,"SYMBOL" ,"TXBIOTYPE" ,"TXID" ,"TXNAME" ,"UNIPROTID")}

\item{beta.thr}{numeric, beta difference threshold value. Default = 0.3.}

\item{overlap.param.thr}{numeric, threshold value for selected parameter to filter methylations overlapping the selected features. Default = 100}

\item{param.type}{character, threshold parameter to filter methylations overlapping the selected features. Accetped c("dmr.length", "overlap.length", "overlap.percentage"). Default = "overlap.length".}
}
\value{
data.frame with TF methylation levels, target.genes expression
}
\description{
For each TF find and associate the target genes, within the annotation results, and retrieve expression
}
