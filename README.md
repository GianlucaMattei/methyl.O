# Methyl.O
## Annotate, Score and Visualize Differentially Methylated Regions

## Introduction
During the last years the adaptation of next generation sequencing technologies to epigenetic studies has revolutionized the ability to study the methylation status of DNA. Despite the development of new techniques and the increasing interest in epigenomic, the landscape of available tools does not offer a wide range of possibilities for annotation. Here we present methyl.O, a R package and an online tool focused on interpretations of methylation analysis. While the practical user interface benefits non- bioinformatic users, the R package allows those more experienced to perform methylation analyses with more in-depth control over parameters. The user interface can also be run in local by R, to facilitate preliminary analyses and to easily visualize data distribution and methylated gene regions. A peculiarity of methyl.O results is the ranked list of the genes most affected by methylation and, according to available databases, most involved in pathologies. The ranked list suggests the genes to focus on to understand methylation effects and to consider for further analyses. Furthermore methyl.O allows, when expression profiles are provided, to compare transcriptomics to methylome to have a deeper insight into the resulting phenotype. methyl.O also offers the possibility to perform enrichment analysis, based on the package EnrichR, for a better interpretation of results. Finally, among the peculiarities that we implemented in our tool, the possibility to annotate enhancers and TFs to genes in order to predict the effects on expression. Moreover, both enhancers and TFs results can be integrated with expression profiles for a better understanding of the impact that methylation of these elements may lead to.

## Getting started
#### Installation:
The R release version of methyl.O is available via Bioconductor and GitHub. 
It can be installed as follow:

Bioconductor:

```Bioconductor
if (!requireNamespace("BiocManager", quietly=TRUE)) 
install.packages("BiocManager") 
BiocManager::install(methyl.O)
```

GitHub:

```GitHub
devtools::install_github(GianlucaMattei/methyl.O)
```

#### Run interface:
Once the package is installed, using the function `runOnDesktop()` it is possible to initialize the user interface. 

#### Online interface:
methyl.O is also available online at: “INSERIRE DOMINIO”


## Required inputs:
#### **R:**	 
To run methyl.O you need a data.frame object with the genomic coordinates (chr, start and end) on the first three columns and the beta difference values on the fourth column. Specific beta values for the samples, used to score the beta difference values, are optional but recommended. The Beta difference values must be a fraction of 1. Additional columns will be stored in output under the "others" field

seqnames|	start    	|	end  		|   betadiff	|   betaTum 	|   betaHlt	|   pval
--------|---------------|---------------|---------------|---------------|-----------|--------
chr10	|	38300277	|   38300763	|   -0.355	    |   0.474		|   0.031	|   0.030
chr10	|	43725576	|   43725884	|   0.301	    |   0.333		|   0.8125	|   0.010
chr10	|	50886863	|   50887851	|   0.396	    |   0.222		|   0.785	|   0.001
chr10	|	105343883	|   105344827	|   0.315	    |   0.058		|   0.35	|   0.003
chr10	|	105728786	|   105730463	|   -0.386	    |   0.75		|   0.181	|   0.032
chr10	|	115999356	|   115999894	|   0.384	    |   0		    |   0.636	|   0.020




#### **GUI:**	
Otherwise in the GUI’s Homepage is possible to upload custom files specifying the delimiter type and if it already includes the names of the columns (header) and, if expressed as a percentage , to convert the beta difference values to a fraction of 1.

## Basic Concepts
#### Gene Model and Gene’s Features:
Each is based on the classical gene model, which includes the promoter, the 5’/3’ UTRs and the exons/introns, and additional features specific to the methylation: the head of the gene and the TSS surrounding region. While the head corresponds to the first part of the gene, the TSS surrounding region results from considering the promoter and the head as a single range where methylation may have a greater effect on gene expression. As shown in the figure 1, by considering TSS surrounding or head regions for the annotation it is possible to avoid overlapping features. In fact, in some cases, the 5’ and the 3’ UTRs may overlap part of the first and last exon respectively (A), they can coincide with the first and the last exon or they can be absent (B) and in other cases they can overlap more than one exons (C). Thus considering only the 5’UTRs will exclude some genes while considering both 5’UTR and the first exons will lead to redundant results.

<img src="https://github.com/GianlucaMattei/methyl.O/blob/main/vignettes/genemodels.png" alt=" **Figure 1**: Schematic gene model " width="1280" height="720">

![**Figure 1**: Schematic gene model](https://github.com/GianlucaMattei/methyl.O/blob/main/vignettes/genemodels.png)

## Annotate DMRs:
#### **R:**	
Genes and gene’s features can be annotated by methyl.O using the function annotateDMRs. An example dataset is included in the package and can be loaded by

```DMRsSubset
data("DMRsSubset", package="methyl.O")
```

The function uses as default options Ensembl database to annotate and the hg19 version of human assembly. It is possible also to use hg38 for Ensembl or UCSC database and hg19 to annotate. At the moment UCSC can not be used to annotate the hg38 version of the human assembly. The settings for annotation permit to customize results. Running annotateDMRs allows: - to modify the length of promoters (prom.length, default 1500bp) and the length of the heads (head.length, default 1500bp), - to decide whether to use the longest transcript for each gene or to perform the analysis on all transcripts for each gene (longest.trx, default TRUE), - to define the database to be used for annotation (annotation, default = Ensembl), - to set the assembly version (hg, default hg19), the beta threshold for each DMRs to be considered in the annotation (thr.beta, default = 0.3), the percentage threshold for CpG Islands (CGIs) to be considered differentially methylated(thr.cgis, default = 0.4) and the percentage threshold for enhancers to be considered differentially methylated (thr.enhancer, default =0.4). Annotation.fast option concerns the method to return the gene Symbols. In fact, the ranges used as input are first annotated as Ensembl gene IDs, if Ensembl is used, or as Entrez gene IDs if UCSC is used, which in turn can be converted 1:1 to Symbols, speeding up the process, or can be 1:many (or many:1) for more accurate IDs conversion (annotation.fast, default = TRUE). The last three options permit, if the input data are not correctly formatted, to select the column number where to find the beta difference (col.betadiff, default 4), needed to perform the analysis and the columns of beta values of the two compared samples (col.beta1 and col.beta2, default = NULL), recommended but not necessary. 

```annotatedDMRs
annotatedDMRs <- annotateDMRs(DMRsSubset , prom.length=1500, head.length=1500, longest.trx=TRUE, 
annotation='ensembl', hg='hg19', annotation.fast=TRUE, thr.beta=.3, thr.cgis=.4, thr.enhancer=.4, 
col.betadiff = 4, col.beta1 = NULL, col.beta2 = NULL)
```

#### **GUI:**
In the “Annotate Methylated Regions” tab is possible to annotate the DMRs.The same page also displays the ranked list of the genes with scores, but this function will be explained in the next paragraph. On the left side are placed the main settings, while e the results are displayed in the main part of the pag. 

Command| Description |
--|--
Select Annotations to Use | specifies the database to use for annotation  
Select Assembly version | Set the assembly version to use 
Use Longest transcript | Yes to use the longest transcript for each gene or use all the transcripts for each gene 
Compute Fast Annotation | Select the type of conversion to gene symbol from ensembl gene IDs, if ensemble is used, or from entrez gene IDs if UCSC is used, It can be setted to 1:1, speeding up the process, or can be 1:many (or many:1) for more accurate IDs conversion 
Select Promoter Lengths | Set the promoter lengths 
Select Head Lengths | Set the head lengths 
Length Percentage of Altered Methylated CGIs | Set the length, in percentage, of overlaps between the the DMRs and the CGIs to be considered during the analysis 
Length Percentage of Altered Methylated Enhancers | Set the length, in percentage, of overlaps between the the DMRs and the enhancers to be considered during the analysis 
Beta Diff. Threshold | The beta threshold for each DMRs to be considered in the annotation 
Column Position of Beta Diff. | Select the column number where to find the beta difference 
Column's Position of Beta Values of Sample 1 | Select the column number where to find the beta values of sample 1
Column's Position of Beta Values of Sample 2 | Select the column number where to find the beta values of sample 1

**Table 2**: Parameters for Annotate Methylated Regions tab


In the GUI we also implemented some plots to have a quick overview on results. These plots are in order: the beta difference values distribution by chromosome, the number distribution of DMRs by chromosome, the distribution of widths, ranks, beta values and database scores, for each annotation (Genes, Heads, TSS surrounding, Promoters, Exons, 5’UTRs, 3’UTRs, Introns). Two additional annotations can be plotted (First exons, First introns). The last plot compares the number of DMRs for each annotated region. Each of these plots has some specific customization settings by clicking on the red gear button. Other common settings are placed in the left column of the page:

Command | Description
-|-
Met Width Min | Minimum methylation width threshold in bp 
Met Width Max | Maximum methylation width threshold in bp
Selected Feature Percentage Min | Minimum methylation width threshold in %
Selected Feature Percentage Max | Maximum methylation width threshold in %
Feature Rank Min | Minimum rank threshold 
Feature Rank Max | Maximum rank threshold

: (\##tab:table3)

**Table 3:** Additional parameters for plots and for results table in Annotate Methylated Regions tab

The last two options, Feature Rank Min and Max refer to the ranking position of the features in the gene model. For example, considering the exons, the first is rank 1, the second rank 2 and so on. The same for the introns. Of course these parameters are not useful for other features as promoters, heads and other non rankable features, thus they will not affect the plots related to these features. 



## Score Methylated Regions:
methyl.O implements a customizable score system aiming to suggest the most impacting DMRs. This score system is based on both databases and overlapped regions. In fact, the databases return the implication of the overlapped genes in certain pathologies as well as the presence of regulatory elements as CGIs and TF, while the different overlapped regions may affect the expression in different ways. The returned score aims to integrate these two informations. Overlapped regions affecting the gene expression can be set and must be one or more elements from annotation results, first exons and first introns. The two pieces of information can be weighted by the option score modifier. The score modifier ranges from 0 to 1, where 0 returns a score based on the database only and 1 a score based on the overlapped features, focusing the results on the effects of DMRs on gene expression.

#### **R:**	
The function scoreAnnotatedDMRs accepts as input object the list resulting from annotateDMRs, as input options active.features which specifies the features considered to affect most the gene expression and the score.modifier,

```annotatedDMRs
annotatedDMRs <- scoreAnnotatedDMRs(annotatedDMRs, active.features = c("promoters", "heads"), score.modifier = 0.5)
```

#### **GUI:**	
The score is computed automatically during the annotation process and is displayed in the Annotate Methylated Regions tab. Active features can be selected by the red dashboard at the top of the page while the score weights can be modified by the slider on the left panel.

## Results of annotateDMRs and scoreAnnotatedDMRs:
#### **R / GUI:**	
In R the results from annotations are stored in a list object where each element corresponds to a feature. The same object is shown in the GUI where the elements, therefore each feature, can be displayed by clicking on the red gear button. The results shown are affected by parameters explained in the table 3. An example of results retrieved by R or by the GUI, is shown below. The table shows the annotations for genes, the first element of the list of results, but the same columns can be found in the other elements of the list.  


seqnames|start|end|width|beta|gene.start|gene.end|gene.width|gene.strand|gene.id|tx.namegenes.perc|tag|dgv|gnomad|NCG|NCG_type|COSMIC|cosmic.CGIs|hacer|CGIs|TF|symbol|database.scoreothers
-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|
chr7|1894240|1895255|1016|-0.4587633|1855430|2274378|418949|-|ENSG00000002822|ENST00000406869|0.2425116|chr7_1894240_1895255|1|1|0|0|0|1|0|1|MAD1L1|3|betaHlt=0.831983805668016;betaTum=0.118055555555556;pval=0.000343062585523649
chr17|26712105|26712449|345|0.4203297|26689878|26725151|35274|+|ENSG00000004139|ENST00000379061|0.9780575|chr17_26712105_26712449|1|1|0|0|0|0|0|1|SARM1|3|betaHlt=0.150375939849624;betaTum=0.541958041958042;pval=0.00132628844389729
chr17|42462135|42462494|360|0.4608637|42449550|42468373|18824|-|ENSG00000005961|ENST00000353281|1.9124522|chr17_42462135_42462494|1|1|0|0|0|0|0|1|ITGA2B|3|betaHlt=0.14957264957265;betaTum=0.822649572649573;pval=0.0449512427024494
chr16|3067875|3068113|239|0.3508772|3066946|3072087|5142|+|ENSG00000006327|ENST00000573001|4.6479969|chr16_3067875_3068113|1|1|0|0|0|0|0|1|TNFRSF12A|3|betaHlt=0.178571428571429;betaTum=0.601503759398496;pval=0.000310800310800311
chr1|55266381|55266722|342|0.3076923|55245385|55268440|23056|-|ENSG00000006555|ENST00000371276|1.4833449|chr1_55266381_55266722|1|1|0|0|0|0|1|1|TTC22|3|betaHlt=0.0588235294117647;betaTum=0.533333333333333;pval=0.00369137539839573
chr19|42056859|42057411|553|-0.4423077|42054386|42093196|38811|+|ENSG00000007129|ENST00000407170|1.4248538|chr19_42056859_42057411|1|1|0|0|0|0|0|1|CEACAM21|3|betaHlt=0.777777777777778;betaTum=0.25;pval=0.0080063295572896

**Table 4:** Resulting annotation for genes

The first five columns (seqnames, start, end, width and beta) refer to the annotated DMRs, symbol refers to the overlapping gene, score is the ranking score assigned to the DMRs and is specific for genes only, others contains all the additional information found in the input table, gene start, gene end, gene width and gene strand are coordinates and characteristics regarding the overlapped gene. For other elements these characteristics will be specific for the selected feature For example for introns we will have intron start, intron end, intron width and intron strand. The second part of the table contains the gene id and transcript name according to the database used for annotations, genes perc is the percentage of the gene (or the current feature selected) overlapped by the DMR, tag is an ID referred to the DMR, dgv, gnomad, NCG, NCG type, COSMIC, hacer, CGIs and TF are the information retrieved from database querying, where 1 is used when the DMR’s range is present in the database, and finally database score is the computed score for database as described in the paper. According to the selected element/feature, additional columns are shown: these are one returning the overlap in bp between the DMR and the feature, and rank which returns the position of the annotated element in the gene model. Moreover the GUI permits to sort the displaying table by clicking the desidered column and to search words, as the gene ids or transcript ids, within the table.



## Annotate Methylated Enhancers:
In order to assess the effect of methylation on gene expression, enhancers cannot be excluded from the analysis. Enhancer elements are regions along the DNA able to mediate the recruitment of TFs enhancing the expression of distal genes; for this reason we implemented the possibility to annotate this type of regions. The annotation is based on Hacer database which in turn is based on FANTOM5 and 4DGenome for enhancer-gene association. The resulting table is similar to the one in figure 2 and also in this case the returned results are ranked by a score that can be weighted by the option score modifier. The score modifier ranges from 0 to 1, where 0 returns a score based only on the database and 1 on the overlapped enhancers shifting the results on the effect of DMRs on gene expression. Other information returned are the presence (1 or 0) of the enhancer associated genes, with the corresponding gene symbols, in the NCG and COSMIC databases. 

#### **R:**	
The DMRs can be annotated by the function annotateEnhancers which accepts as input the data.frame object with the genomic coordinates (chr, start and end) on the first three columns and a beta difference value on the fourth column. Other options permit to specify the genome assembly version (hg), the beta difference value threshold for filtering DMRs (thr.beta), the type of metric to filter the DMRs, the threshold value for the selected metric, the score modifier (score.modifier) which works as previously explained. The metrics to filter the DMRs can be 1) the percentage of the enhancer overlapped with the DMR, 2) the length, in bp, of the overlap between the enhancers and 3) the DMR and the length, in bp of the DMR. Finally, the column where to find the beta difference values can be set in col.betadiff option. 

```annotatedEnhancers
annotatedEnhancers <- annotateEnhancers(DMRsSubset, hg='hg19', thr.beta = 0.3, overlap.param.thr = 40, 
param.type = "overlap.percentage", score.modifier = 0.5, col.betadiff = 4)
```

#### **GUI:**	
The tab Annotated Methylated Enhancers returns a table similar to that in Annotate Methylated Regions. In the left panel the customizable options are in order

Command | Description
-|-
Type of Parameter to filter Overlapping Methylation | The metrics to filter the DMRs can be 1) the percentage of the enhancer overlapped with the DMR, 2) the length, in bp, of the overlap between the enhancer and 3) the DMR and the length, in bp of the DMR.
Overlapping Methylation Value Threshold for Filtering | Threshold value for the selected parameter.
Beta Diff. Threshold | The beta threshold for each DMRs to be considered in the annotation
Score Modifier | Ranging from 0 to 1, where 0 returns a score based only on the databases and 1 on the characteristics of overlaps.
Column Position of Beta Diff. | Select the column number where to find the beta difference

**Table 5:** Parameters of Annotate Methylated Enhancers tab

Additional plots in the GUI are returned: one for the selected parameter distribution and one for the beta difference values distribution. Both are customizable by the red gear button above the plots.




## Visualize Methylations:
The package includes the option to visualize one or more DMRs occurring on one gene, displaying the features, the corresponding transcripts, the CGIs and the beta differences. The graphic representation helps to better evaluate the effect of the methylated segments on gene expression. 


![**Figure2**: Starting from the top of the plot are shown: the chromosome location, the CGIs ranges, the promoter range and its direction, one or more transcripts associated with the gene and the information about the DMR including the beta difference value and length.  ](https://github.com/GianlucaMattei/methyl.O/blob/main/vignettes/visualizeMethylation.png)




#### **R:**	
To visualize the desired gene, the function plotDMRs needs the resulting annotated list from annotateDMRs and the symbol ID of the gene. Other optional parameters permit to select the database to use between Ensembl and UCSC (database), the genome assembly version (hg), the width of the promoter (prom.width), the possibility to show all transcripts or the longest (show.all.transcripts) and to select the genomic range to show (coord.zoom). By this function it is possible to plot the beta values of samples in order to evaluate the methylation levels between two conditions. In this case this information must be included in the annotated list and in the initial input data.frame of ranges. To plot the single beta values for each sample it is necessary to indicate the name of the columns in the initial input data.frame (beta1.name, beta2.name) which can be found also under the column “other” in the annotated list object. For example by using the provided dataset you can use beta1.name=”betaTum”, beta2.name=”betaHlt”. Other graphical parameters permit: - to set the colors of the plotted beta values (beta.colors) if the beta difference value is plotted beta.colors then one color must be specified for beta.colors otherwise two must be chosen, - to set if the output plot should be in greyscale (blackandwhite), - to decide if the plotted range should be zoomed out for a better visualization (smartzoom) and other settings to save the plot including the path (path) the height and width of the output pdf (height.pdf, width.pdf).

```plotDMRs
plotDMRs(annotatedDMRs,”CTF1”)
```

#### **GUI:**	
The tab Visualize Methylation permits to plot the DMR. The main part of the page displays the plot as shown in figure 5. In the left panel are noted the customizable settings 

Command | Description
-|-
Gene Symbol | The gene symbol to be displayed
Beta 1 Name / Beta 2 Name | In order to plot the singles beta values for each sample it is necessary to indicate the name of the columns in the initial input data.frame (beta1.name, beta2.name) which can be found also under the column “other” in the annotated list object
Select Beta Color | Sets the color of the beta difference value or the color for beta value for sample 1 
Select Beta 2 Color | Sets the color for beta value for sample 2
Show All Transcripts| If TRUE all transcripts are shown for the selected gene
Plot in Grayscale|If TRUE the plot will be displayed in grayscale 
Smart Zoom|Zoom out the range for a better visualization
Promoters Length|Sets the promoter length
Zoom Coordinate - From / to |Sets custom coordinates for the range to display

**Table 6:** Parameters for Visualize Methylation tab 


## Expression Data Integration:
In order to evaluate how methylation affects gene activity it is important to integrate expression profiles, when available. methyl.O permits the integration of expression data and computes both correlation and a ranking score. The integration process assigns a score to each gene based on the beta value and the log. fold change. The higher the values the higher is the score, positive in case the beta and the expression values are discordant, negative if concordant. It is known that hyper-methylated sequences repress gene expression, however hypo-methylated sequences can not be directly associated with upregulation of genes. In fact, hyper-methylation usually is a long term regulation involved in cells' fate commitment, while hypo-methylated genes still have the possibility to be expressed but undergo several types of other regulation mechanisms. For this reason methyl.O can also compute the correlation between the expression and the beta differences of hyper-methylated and hypo-methylated genes separately as shown in figure 6. methyl.O allows to plot the correlations and to assign to each methylated gene its expression. 

![**Figure 3**:  Correlation plot between methylation and expression. Linear regression of gene expression and hyper-methylated genes is represented by the orange dashed line. The grey dashed line represents the linear regression of gene expression with hypo-methylated genes. In this case we can observe an inverse correlation between methylation levels and expression in both cases](https://github.com/GianlucaMattei/methyl.O/blob/main/vignettes/MethylationVsExpression.png)

#### **R:**	
The function associateFeat2Exprs, needs the annotated DMRs and the normalized expression profile to perform integration of data. The function parameters permit to select the methylated features considered most affecting the expression (features), to indicate the column of the expression file where to find the gene IDs (col.genes), the statistics to use (col.stat), the log. fold change (col.logFC), the statistical (stat.thr) and the log. fold change (logfc.thr) thresholds. You can also set a threshold for the beta value (beta.thr), the metrics to filter DMRs (param.type), as explained in Annotate Methylated Enhancers paragraph, and the threshold value (overlap.param.thr). The correlation is customizable by setting the type of correlation to compute and plot (plot.type), and the method to compute the correlation coefficient (cor.type). The function accepts by default the symbol IDs but can also convert automatically other types of IDs to symbols (convert.genes). In this case it is necessary to indicate from which type of IDs they must convert the genes to symbols. The function can return a data.frame (return.table) where the beta difference and the log. fold change are associated with each gene, otherwise the function can return a plot displaying the two correlations. The customizable parameters permit to set the correlation line colors (lmfit.col1, lmfit.col2), the axis lines colors (line.col), the palette reflecting the score values for the plot (pal) and whether show the gene names next to each dot (show.text). For both data.frame and plots it is possible to filter the genes (filter.by.genes) to study how the methylation can affect for instance the gene expression of a specific pathway. 

```expressionSubset
data(expressionSubset)
```

```associateFeat2Exprs
associateFeat2Exprs(annotatedDMRs, features = c("promoters", "heads"), expressionSubset, col.genes =1 , 
col.stat = 6, stat.thr = 0.05, col.logFC = 2, logfc.thr = 0, beta.thr = .3, plot.type = 'splitted', 
cor.type = 'pearson', return.table = FALSE)
```



#### **GUI:**	
The tab Methylation vs Expression permits to integrate expression profiles to methylation data. In the main page are shown both the plot of correlations and the data frame where each annotated gene is associated with the beta value and the log. fold change. The upper part of the page permits to select the methylated gene regions to associate to the expression while in the left panel it is possible to upload the expression file. The following parameters are available:

Command | Description
-|-
Type of Parameter to filter Overlapping Methylation | The metrics to filter the DMRs can be 1) the percentage of the enhancer overlapped with the DMR, 2) the length, in bp, of the overlap between the enhancer and 3) the DMR and the length, in bp of the DMR.
Overlapping Methylation Value Threshold for Filtering | Threshold value for the selected parameter.
Statistic Threshold | Threshold value for the selected statistics.
LogFC Threshold | Threshold value for expression log. fold change.
Beta Diff Threshold | The beta threshold for each DMRs to be considered in the annotation
Filter by DB| Filters the genes from a specific database, resulting from enrichment analysis, to study how the methylation can affect the gene expression.
Select Correlation Type | Set the method to compute the correlation
Select Path for filtering Genes | Filters the genes to study how the methylation can affect the gene expression of a specific pathway resulting from enrichment analysis.| Column Position of Gene ID
Set the column position in the expression file where to find gene IDs | Column Position of Used Statistics
Set the column position in the expression file where to find the statistics (p.value or adj p.value) | Column position of logFC
Set the column position in the expression file where to get the log. fold changes values. | Select TRUE if Gene IDs are not Symbols
If selected, gene IDs are not official symbols| Select Gene IDs Annotation Type to Translate | If the above option is TRUE, then select the type of IDs provided.

**Table 7:** Parameters for Methylation vs Expression tab

Other settings to customize the plot can be found in the red gear button above the plot. Here can be also found the option to compute splitted correlation for hyper-methylated and hypo-methylated genes.


## Integration of Expression Data with Annotated Enhancers:
We implemented in methyl.O the possibility to investigate by data integration the effects on gene expression of methylation levels of enhancer elements. The results of this analysis can be plotted or retrieved in a data frame. The data frame contains informations about the genomic coordinates of DMRs (seqnames, start, end), about the coordinates of enhancer (annot.seqnames, annot.start, annot.end), the target genes (target.gene), the beta difference occurring on DMRs (beta), the expression log. fold change (logFC) and the score (score). The score is computed as described in Expression Data Integration. As introduced in Annotate Methylated Enhancers, also in this case we introduced the possibility to correlate the hyper and hypo-methylated enhancers separately. In fact. hyper-methylation of these regions prevents the binding with TFs and consequently the expression enhancement of distal genes but expression of distal genes can still be mediated by TFs binding to their promoters. 

#### **R:**	
The function associateEnh2Expr is used to correlate the annotated enhancers resulting from the function annotateEnhancers to expression data. The genome assembly needs to be specified in the option (hg) as well as the enhancer database to be used (enhancer.db). As for the associateFeat2Exprs function, it is possible to specify the column in the expression file under which to find the gene IDs (col.genes), the log. fold change (col.logFC), the column where to find the statistics (col.stat), choose whether gene IDs need to be converted (convert.genes), define the ID type of genes to convert (convert.from) and set the threshold for statistical significance and for the log. fold change. Other parameters permit to customize the methylation characteristics to filter the DMRs: the beta difference threshold (beta.thr), the method to compute the overlap between the DMRs and the enhancers (param.type) and the value threshold for overlaps (overlap.param.thr). The remaining options, that allow to define the methods to compute the correlation, to set the graphical parameters and if return a table or a plot, work as in the associateFeat2Exprs function where they have been described in detail.

#### **GUI:**	
The tab Methylated Enhancers vs Expression permits to correlate the methylation levels of enhancers to expression profile of target genes. In the main part of the page are shown the plot and the table of results. In the left panel is possible to load the expression profile file and are shown the following options:



Command | Description
-|-
Type of Parameter to filter Overlapping Methylation | The metrics to filter the DMRs can be 1) the percentage of the enhancer overlapped with the DMR, 2) the length, in bp, of the overlap between the enhancer and 3) the DMR and the length, in bp of the DMR. | Overlapping Methylation Value Threshold for Filtering | Threshold value for the selected parameter.
Statistic Threshold |Threshold value for the selected statistics.
LogFC Threshold|Threshold value for expression log. fold change.
Enhancer Methylation Beta Threshold|The beta threshold for each DMRs to be considered in the annotation
Select enhancer DB|Select the DB to use
Select Correlation Type|Set the method to compute the correlation
Column Position of Gene ID|Set the column position in the expression file where to find gene IDs
Column Position of Used Statistics | Set the column position in the expression file where to find the statistics (p.value or adj p.value)
Column position of logFC|Set the column position in the expression file where to get the log. fold changes values.
Select TRUE if Gene IDs are not Symbols |If selected, gene IDs are not official symbols
Select Gene IDs Annotation Type to Translate|If the above option is TRUE, then select the type of IDs provided.

**Table 7**: Parameters for Methylated Enhancers vs Expression tab


## Enrichment of Methylated Regions:
Since the deregulation of one gene expression unlikely affects the final phenotype, the methylation profile effect should be evaluated as a whole in order to retrieve as small changes can cooperate to emerge in the final results. For this reason we implemented in methyl.O the possibility to perform enrichment analyses using 12 databases. Since the effects on gene expression depend on the overlapped regions, genes to query the databases can be chosen based on the annotated features. Therefore it is possible to query only the genes with DMRs overlapping their promoter, or their head, and so on, or a combination of the available features as well as all the genes returned by the annotation. Results from this analysis permit to detect the pathways where the differential methylation brings the major contributes to their perturbation. Moreover the genes from the resulting pathways can be extracted to further evaluate the correlation with methylation. 

#### **R:**	
The function annotatedDMRs2Enrichr is designed to easly perform enrichment analysis, based on the R package enrichR, using as input the annotated list object resulting from annotateDMRs. The second option permits to specify the features considered to select the genes to query the databases. Since the statistic is computed for each enriched pathway, other options define the parameters for filtering the results. Once set which statistic to use (stat.filter) between the p.value, the adj. p.value or the overlap between the queried genes and the genes of the pathway, the threshold can be defined (stat.thr). By default 12 different databases are used by default for the enrichment otherwise it is possible to select specific databases (db). The function plotDMRs2Enrichr permits to plot the results of annotatedDMRs2Enrichr. Two different types of plots can be returned: the barplot and the lollipop plot, both displaying the -log10 of p.value or adj. p.value of the enriched pathways. Therefore the smaller is the p.value the greater is the value shown in the plot. The barplot also returns additional information regarding the fraction of hyper and hypo-methylated genes among the genes that contributed to the enrichment of a specific pathway. This information can be used to better understand if a pathway is repressed, upregulated or just altered in the studied conditions. The function needs as inputs the enrichment results from annotatedDMRs2Enrichr and the annotation results from annotateDMRs. Other parameters permit to set the statistics, between the p.value, the adj. p.value or the overlap (stat), to show and sort the results, to set the number of pathways to display (n) and to set the plot type between barplot and lollipop (plot.type). The graphical parameters define the colors of the hyper and hypo-methylated genes (col.hyper; col.hypo) and the colors palette for lollipop plot. Finally it is possible to plot vertical lines indicating the statistical threshold (thrs) and set their colors (thrs.col)


#### **GUI:**	
The Methylation Enrichment tab performs the enrichment for annotated genes resulting from the Annotate Methylated Region tab returning both a table and a plot of results. Two different types of plots can be returned: the barplot and the lollipop plot, both displaying the number of genes found in the pathway or the -log10 of the p.value or the adj. p.value of the enriched pathways. Therefore the smaller is the p.value the greater er is the value shown in the plot. This value is also used to sort the results. The barplot also returns one additional information regarding the fraction of hyper and hypo-methylated genes among the genes which contributed to the enrichment of a specific pathway. This information can be used to better understand if a pathway is repressed, upregulated or just altered in the studied conditions. Moreover, in order to study how methylation can alter the expression of an enriched pathway, genes from each pathway can be extracted to study the correlation between methylation and expression. This can be done under the tab Methylation vs Expression: on the left panel the “Select Path for filtering Genes” shows a dropdown menu with the resulting pathways from the current analysis. The parameters for the enrichment analysis can be found on the left panel: 


Command | Description
-|-
Statistics for Results Filtering | For each enriched pathway is computed a p. value and a adj. p.value, this parameter permits to set which statistic to use.
Select Annotation from where gene symbols are taken | Select the differentially methylated features from where to pick the genes.
DBs to Query | Select the database to query among the 12 available. If left blank, all 12 databases will be used.
LogFC Threshold | The log. fold change threshold for genes to be considered.
Value Threshold for Filtering by Statistics | Select the statistical threshold for the selected statistics. 

**Table 8:** Parameters for Methylation Enrichment tab

Finally, since the enrichment is performed on results by the annotation process, changing the parameters used for annotation itself will affect the resulting enriched pathways.



## TFs Analysis:
While alteration of methylation on genes alone can lead to little changes, a DMR occuring on TF can cause a major perturbation in overall gene expression profile. For this reason methyl.O offers the possibility to perform specific analyses on TFs. First, the annotated DMRs are used to query the ENCODE database in order to select all TFs overlapped by DMRs, then the same database is used to associate target genes to each TFs. The tool offers also the possibility to compare methylation levels of TFs and expression of target genes to better evaluate their impact (**Figure 4**). Finally the differentially expressed targets can be used to perform enrichment analysis by enrichR.



![**Figure 4**: The plot shows on the left the beta Difference of Methylation between two compared conditions of specific TFs, STAT1 in this case, and the expression of target genes. ](https://github.com/GianlucaMattei/methyl.O/blob/main/vignettes/EnrichrTFsPlot.png)

#### **R:**	
First it is necessary to annotate the TFs to the target genes. This can be done by the function associateTFs2Exprs which accepts as input the annotated DMRs and the expression profile. The specific parameters for input methylation and for expression profiles are the same as used in other functions which integrate expression and methylation data. 

```TFs2Targets
TFs2Targets <- associateTFs2Exprs(annotatedDMRs, features = c("promoters", "heads"), expressionSubset, col.genes = 1,  
col.stat = 6, stat.thr = 0.05, col.logFC = 2, logfc.thr = 0, beta.thr = .3, 
overlap.param.thr = 30, param.type = 'overlap.percentage')
```

Once the association between TFs and targets is retrieved, it is possible to plot the beta difference of a methylation occurring on a specific TF and the expression levels of target genes, as shown in figure 7, using the function plotTFs2Exprs. It needs as input the output of the function associateTFs2Exprs and the TF symbol to plot. The graphical parameters allow to choose the color of the methylation bar and the color palette of the expression values. 

```plotTFs2Exprs
plotTFs2Exprs(TFs2Targets, symbol, col.meth = '##8e0000', pals.bars = 'Cold')
```

The results from the function associateTFs2Exprs can also be used to assess if the methylation of few TFs can affect any pathway by querying enrichR with the target genes. This can be easily done using the function tfs2Enrichr

```tfs2Enrichr
tfs2Enrichr(TFs2Targets, logfc.thr = 1, stat.filter = 'P.value', stat.thr = 0.01)
```


#### **GUI:**	
The TFs Analysis tab permits to study TFs. In the main part of the page it is shown the table of association between TFs and correspective target genes with relative beta differences and expression levels, the overlapped features and the characteristics of the overlap. This table can be customized by the parameters on the upper side of the page which permit to filter TFs according to the overlapped region by the DMRs. Other parameters will be discussed below. In the main page it is possible to plot the methylation of the desidered TFs and the expression of target genes. This plot can be customized by clicking the red gear button on its top left side. The enrichment analysis results are shown in the second part of the main page by a table and a plot that reflect the table results. As previously described, the enrichment plots show on the x-axis the -log10 of p. value, therefore, the smaller the p. value the greater the number shown on the x-axis. Also in this case the plot can be customized by clicking on the red gear button on its top left side. The left side panel is divided in three parts: the first part, named Expression File Parameters, is related to the parameters for using the expression data: 

Command | Description
-|-
Column Position of Gene ID | Set the column position in the expression file where to find gene IDs
Column Position of Used Statistics | Set the column position in the expression file where to find the statistics (p.value or adj p.value)
Column position of logFC | Set the column position in the expression file where to get the log. fold changes values.
Select TRUE if Gene IDs are not Symbols | If selected, gene IDs are not official symbols
Select Gene IDs Annotation Type to Translate | If the above option is TRUE, then select the type of IDs provided.


the second part, named TFs Targets/Expression Correlation Parameters, permits to set parameters for filtering DMRs and expressed genes during the association process:

Command | Description
-|-
LogFC Threshold | Threshold value for expression log. fold change.
Beta Diff Threshold|The beta threshold for each DMRs to be considered in the annotation
Statistic Threshold|Threshold value for the selected statistics
Type of Parameter to filter Overlapping Methylation | The metrics to filter the DMRs can be 1) the percentage of the enhancer overlapped with the DMR, 2) the length, in bp, of the overlap between the enhancer and 3) the DMR and the length, in bp of the DMR.
Overlapping Methylation Value Threshold for Filtering | Threshold value for the selected parameter.

The third part of the left panel, named TFs Targets Enrichment Parameters, is related to the enrichment analysis:

Command | Description
-|-
DBs to Query | Select the database to query among the 12 available. If left blank, all 12 databases will be used.
LogFC Threshold | The log. fold change threshold for target genes to be considered for the association.
Statistics Type for Results Filtering | For each enriched pathway is computed a p. value and an adj. p.value, this parameter permits to set which statistic to use.
Statistics Threshold for Results Filtering | Select the statistical threshold for the selected statistics. 
