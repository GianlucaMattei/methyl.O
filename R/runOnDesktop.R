#' Start Graphical User Interface
#'
#' @return GUI
#'
#' @export

runOnDesktop <- function(){
    options(spinner.color = "#0275D8", spinner.color.background = "#ffffff", spinner.size = 2)
    ui <- shiny::fluidPage(
        shiny::tags$head(shiny::tags$meta(charset="UTF-8")), 
        shiny::tags$meta(name="description", content="Methyl.O is a R package including several utilities for smart approaches, including the integration of expression data, to study the impact of differentially methylated segments of DNA between two conditions. Link to methyl.O repo: www.github.com/GianlucaMattei/methyl.O Link to the browser version of methyl.O: www.genomica.pro"), 
        shiny::tags$meta(name="keywords", content="methyl.o, methylo, DMRs, differentially methylated, expression integration, methylation analysis, methylation software, methylation tool"),
        theme = shinythemes::shinytheme("simplex"),
        shiny::navbarPage(
            "methyl.O",
            shiny::tabPanel("Homepage",
                icon = shiny::icon("home"), 
                shiny::mainPanel(
                    align = "center", width = 8,
                    shiny::imageOutput("logo", height = 200),
                    shiny::h5("Start to navigate to use an example dataset or upload your data", style = "color:black"),
                    shiny::hr(style = "border-top: 0px solid white;"),
                    shiny::hr(style = "border-top: 1px solid #E9EAEA"),
                    shiny::hr(style = "border-top: 0px solid white;"),
                    shiny::splitLayout(
                        shiny::fluidPage(
                            shiny::fileInput(inputId = "bedfile", label = "Upload your file here", placeholder = "select a file", multiple = FALSE),
                            shinyWidgets::prettyRadioButtons("sep", "Select Delimiter Type", c("Comma" = ",", "Space" = " ", "Tab" = "\t"), selected = "\t", shape = "round", inline = TRUE),
                            shiny::checkboxInput("header", "Header", TRUE),
                            shiny::checkboxInput("fraction1", "Are beta Values a Fraction of 1?", TRUE),
                            shiny::tags$style(".btn-file {background-color: red; border-color: red;}"),
                            shiny::hr(style = "border-top: 0px solid #000000;")
                        ),
                        shiny::imageOutput("map", height = 600)
                    ),
                    shiny::tags$footer(
                        shiny::HTML(
                            "<!-- Footer -->
                            <footer class='page-footer font-large indigo'>
                            <div class='footer-copyright text-center py-3'>
                            Methyl.O is a R package including several utilities for smart approaches, including the integration of expression data, to study the impact of differentially methylated segments of DNA between two conditions. Link to methyl.O repo: <a href='https://github.com/GianlucaMattei/methyl.O'> GianlucaMattei/methyl.O</a> Link to the browser version of methyl.O: <a href='https://genomica.pro'> www.genomica.pro</a>
                            Support: gianluca.mattei@unifi.it
                            </div>
                            </footer>
                            <!-- Footer -->"
                        )
                    ),
                    shiny::hr(style = "border-top: 0px solid white;")
                )
            ),


            shiny::tabPanel(
                "Annotate Methylated Regions",
                icon = shiny::icon("map-pin"),
                shiny::sidebarPanel(
                    # PARAMETERS
                    shinyWidgets::actionBttn("annotate", "Annotate", style = 'material-flat', size = 'md', color = 'primary'),
                    shiny::hr(style = "border-top: 0px solid white;"),
                    shinyWidgets::prettyRadioButtons("annotation", "Select Annotations to Use", c("Ensembl" = "ensembl", "UCSC" = "ucsc"), selected = "ensembl", shape = "curve"),
                    shinyWidgets::prettyRadioButtons("hg", "Select Assembly version", c("hg19" = "hg19", "hg38" = "hg38"), selected = "hg19", shape = "curve"),
                    shinyWidgets::prettyRadioButtons("longestTrxs", "Use Longest transcript?", c("YES" = "TRUE", "NO" = "FALSE"), selected = "TRUE", shape = "curve"),
                    shinyWidgets::prettyRadioButtons("annotationFast", "Compute Fast Annotation?", c("YES" = "TRUE", "NO" = "FALSE"), selected = "TRUE", shape = "curve"),
                    shiny::sliderInput("scoreMinMax", label = "Score Modifier", min = 0, max = 1, value = 0.5),
                    shiny::sliderInput("promLength", label = "Select Promoter Lengths", min = 0, max = 5000, value = 1500),
                    shiny::sliderInput("headLength", label = "Select Head Lengths", min = 0, max = 5000, value = 1500),
                    shiny::sliderInput(inputId = "thrCGIs", label = "Length Percentage of Altered Methylated CGIs", min = 0, max = 1, step = 0.01, value = 0.3),
                    shiny::sliderInput(inputId = "thrBeta", label = "Beta Diff. Threshold", min = 0, max = 1, step = 0.01, value = 0.3),
                    shiny::numericInput(inputId = "colBetaDiff", "Column Position of Beta Diff.", value = 4),
                    shiny::numericInput(inputId = "betacol1", "Column's Position of Beta Values of Sample 1", value = NULL),
                    shiny::numericInput(inputId = "betacol2", "Column's Position of Beta Values of Sample 2", value = NULL),

                    # WIDTH
                    shiny::textInput("filterTabWm", "Met Width Min", value = "0"),
                    shiny::textInput("filterTabWM", "Met Width Max", value = "Inf"),

                    # GENE-FEATURE BODY PERCENTAGE
                    shiny::textInput("filterTabPm", "Selected Feature Percentage Min", value = "0"),
                    shiny::textInput("filterTabPM", "Selected Feature  Percentage Max", value = "100"),

                    # FEATURE RANK
                    shiny::textInput("filterTabRm", "Feature Rank Min", value = "0"),
                    shiny::textInput("filterTabRM", "Feature Rank Max", value = "100"),

                    shiny::downloadButton("downloadData", "Download table"),
                    shiny::hr(),
                    width = 2
                ),


                shiny::mainPanel(
                    shinyWidgets::checkboxGroupButtons(
                        inputId = "scoreFeatureSelected",
                        label = "Select Features to Compute the Score",
                        choices = list("Promoters" = "promoters", "Heads" = "heads","TSS Surnd.ing" = "tss.surrounding" ,"Five UTRs" = "fiveUTRs", "Exons" = "exons", "First Exons" = "exons1", "Introns" = "introns","First Introns" = "introns1", "Three UTRs" = "threeUTRs"),
                        status = "primary",
                        selected = c("promoters", "heads"),
                        size = "sm",
                        individual = FALSE,
                        justified = TRUE
                    ),
                    shiny::hr(style = "border-top: 1px solid #000000;"),
                    # main table
                    shinyWidgets::dropdownButton(
                        shiny::selectInput(inputId = "choosenFeature", label = "Table", choices = c("genes", "heads", "tss.surrounding" ,"exons", "fiveUTRs", "threeUTRs", "promoters", "introns"), selected = "genes"),
                        circle = TRUE, status = "primary", icon = shiny::icon("cogs"), size = "sm",
                        tooltip = shinyWidgets::tooltipOptions(title = "Click to see inputs !")
                    ),
                    shiny::hr(style = "border-top: 0px solid white;"),

                    shinycssloaders::withSpinner(DT::dataTableOutput(outputId = "resultsTab"), type = 1),

                    # overview plot 
                    shinyWidgets::dropdownButton(
                        shinyWidgets::prettyRadioButtons("methylationOverview", "Overview Plot", c("Violin Plot" = "violin", "Boxplot" = "boxplot"), selected = "violin", shape = "curve"),
                        shiny::selectInput(inputId = "pals", label = "Overview Plot Palette", choices = c(hcl.pals(), "BlackWhite"), selected = "Dark 2"),
                        circle = TRUE, status = "primary", icon = shiny::icon("cogs"), size = "sm",
                        tooltip = shinyWidgets::tooltipOptions(title = "Click to see inputs !")
                    ),

                    shiny::hr(style = "border-top: 0px solid white;"),
                    shinycssloaders::withSpinner(shiny::plotOutput(outputId = "overviewPlot"), type = 1),
                    shiny::hr(style = "border-top: 0px solid white;"),
                    shiny::downloadButton("downloadOverviewPlot", "Download Plot"),
                    shiny::hr(style = "border-top: 0px solid white;"),
                    shiny::hr(),
                    
                    # chr distribution
                    shinyWidgets::dropdownButton(
                        circle = TRUE, status = "primary", icon = shiny::icon("cogs"), size = "sm", tooltip = shinyWidgets::tooltipOptions(title = "Click to see inputs !"),
                        shiny::selectInput(inputId = "col1", label = "Select Bars Color", choices = c(hcl.pals(), "BlackWhite"), selected = "Mint")
                    ),
                    shiny::hr(style = "border-top: 0px solid white;"),

                    shinycssloaders::withSpinner(shiny::plotOutput(outputId = "resultsPlotChr"), type = 1, color = getOption("spinner.color", default = "red")),
                    shiny::hr(style = "border-top: 0px solid white;"),
                    shiny::downloadButton("downloadResultsPlotChr", "Download Plot"),
                    shiny::hr(style = "border-top: 0px solid white;"),
                    shiny::hr(),

                    # feature overview
                    shinyWidgets::dropdownButton(
                        shiny::selectInput(inputId = "col2", label = "Select Bars Color", choices = c(hcl.pals(), "BlackWhite"), selected = "Mint"),
                        shiny::selectInput(inputId = "choosenFeatureHist", label = "Table", choices = c( "Genes"="genes","Heads"="heads", "TSS Surnd.ing" = "tss.surrounding" , "Exons"="exons", "Five UTRs"="fiveUTRs", "Three UTRs"="threeUTRs", "Promoters"="promoters", "Introns"="introns", "First Introns"="introns1", "First Exons"="exons1"), selected = "genes"),
                        shiny::selectInput(inputId = "choosenValueSWRHist", label = "Barplot Feature", choices = c("Width"="width", "Rank"="rank","Beta"="beta", "Database Score"="database.score"), selected = "width"),
                        shiny::textInput(inputId = "histoBreaks", label = "Histogram Bins", value = 250),
                        circle = TRUE, status = "primary", icon = shiny::icon("cogs"), size = "sm",
                        tooltip = shinyWidgets::tooltipOptions(title = "Click to see inputs !")
                    ),

                    shinycssloaders::withSpinner(shiny::plotOutput(outputId = "resultsPlotSWR"), type = 2),
                    shiny::hr(style = "border-top: 0px solid white;"),
                    shiny::downloadButton("downloadResultsPlotSWR", "Download Plot"),
                    shiny::hr(style = "border-top: 0px solid white;"),
                    shiny::hr(),


                    # pie
                    shinyWidgets::dropdownButton(
                        shinyWidgets::checkboxGroupButtons(inputId = "choosenFeaturePie", label = "Add Additions Features?", choices = c("No" = "No", "First Introns"="introns1", "First Exons"="exons1"), selected = "No"),
                        shiny::selectInput(inputId = "piePalette", label = "Select Colors Palette", choices = hcl.pals(), selected = "ag_GrnYl"),
                        shiny::selectInput(inputId = "featureDistType", label = "Select Chart Type", choices = c("Pie", "Barplot"), selected = "Pie"),
                        circle = TRUE, status = "primary", icon = shiny::icon("cogs"), size = "sm"
                    ),

                    shiny::plotOutput(outputId = "resultsPie"),
                    shiny::hr(style = "border-top: 0px solid white;"),
                    shiny::downloadButton("downloadResultsPie", "Download Plot"),
                    shiny::hr(style = "border-top: 0px solid white;"),
                    shiny::hr()

                )
            ),

            shiny::tabPanel(
                icon = shiny::icon("dna"),
                "Annotate Methylated Enhancers",
                shiny::sidebarPanel(
                    # PARAMETERS
                    shinyWidgets::actionBttn("annotateEnh", "Annotate", style = 'material-flat', size = 'md', color = 'primary'),                    
                    shiny::hr(style = "border-top: 0px solid white;"),
                    shiny::selectInput(inputId = "paramTypeEnhAnnot", label = "Type of Parameter to filter Overlapping Methylation", choices = c("DMR Length" = "dmr.length", "Overlap Length" = "overlap.length", "Overlap Perc." = "overlap.percentage"), selected = "overlap.percentage"),
                    shiny::numericInput("overlapParamThrEnhAnnot", "Overlapping Methylation Value Threshold for Filtering", value = 40),
                    shiny::sliderInput(inputId = "thr.beta.enhancer", label = "Beta Diff. Threshold", min = 0, max = 1, step = 0.05, value = 0.3),
                    shiny::sliderInput(inputId = "scoreModifierenh", label = "Score Modifier", min = 0, max = 1, step = 0.05, value = 0.5),
                    shiny::numericInput(inputId = "colBetaDiff", "Column Position of Beta Diff.", value = 4),
                    shiny::downloadButton("downloadDataEnhancer", "Download table"),
                    shiny::hr(),
                    width = 2
                ),

                shiny::mainPanel(
                    shinycssloaders::withSpinner(DT::dataTableOutput(outputId = "resultsTabEnhancer"), type = 1),
                    shiny::hr(style = "border-top: 0px solid white;"),
                    shiny::hr(),

                    shinyWidgets::dropdownButton(
                        shinyWidgets::prettyRadioButtons("plotDistEnhancer", "Enhancer Plot Type", c("Density" = "density", "Histogram" = "histogram"), selected = "density", shape = "curve"),
                        colourpicker::colourInput("colDistEnhancer", "Select Color", "#20f099"),
                        shiny::sliderInput(inputId = 'binDist', "Histogram Bins", min = 2, max = 1000, step = 10, value = 100),
                        circle = TRUE, status = "primary", icon = shiny::icon("cogs"), size = "sm"
                    ),
                    shinycssloaders::withSpinner(shiny::plotOutput(outputId = "distEnhancer"), type = 1),
                    shiny::downloadButton("downloadDistEnhancer", "Download Plot"),
                    shiny::hr(style = "border-top: 0px solid white;"),
                    shiny::hr(),

                    shinyWidgets::dropdownButton(
                        shinyWidgets::prettyRadioButtons("plotDistBeta", "Beta Plot Type", c("Density" = "density", "Histogram" = "histogram"), selected = "density", shape = "curve"),
                        colourpicker::colourInput("colDistBeta", "Select Color", "#3a77a0"),
                        circle = TRUE, status = "primary", icon = shiny::icon("cogs"), size = "sm"
                    ),
                    shinycssloaders::withSpinner(shiny::plotOutput(outputId = "distBeta"), type = 1),
                    shiny::downloadButton("downloadDistBeta", "Download Plot"),
                    shiny::hr(style = "border-top: 0px solid white;"),
                    shiny::hr()

                )
            ),

            shiny::tabPanel(
                "Visualize Methylation",
                icon = shiny::icon("search"),
                shiny::sidebarPanel(
                    shinyWidgets::actionBttn("visualizeMeth", "Visualize", style = 'material-flat', size = 'md', color = 'primary'),                    
                    shiny::hr(style = "border-top: 0px solid white;"),
                    # PARAMETERS
                    shiny::textInput(inputId = "symbolVisualize", label = "Gene Symbol", value = ""),
                    shiny::textInput(inputId = "betaName1", label = "Beta 1 Name", value = NULL),
                    shiny::textInput(inputId = "betaName2", label = "Beta 2 Name", value = NULL),
                    colourpicker::colourInput("colbeta1", "Select Beta Color ", "red"),
                    colourpicker::colourInput("colbeta2", "Select Beta 2 Color ", "#3a77a0"),
                    shinyWidgets::prettyRadioButtons("showAllTrxs", "Show All Transcripts", c("TRUE" = "True", "FALSE" = "False"), selected = "False", shape = "curve"),
                    shinyWidgets::prettyRadioButtons("blackwhite", "Plot in Gray Scale", c("TRUE" = "True", "FALSE" = "False"), selected = "False", shape = "curve"),
                    shinyWidgets::prettyRadioButtons("smartZoomTrxs", "Smart Zoom", c("TRUE" = "True", "FALSE" = "False"), selected = "True", shape = "curve"),
                    shiny::sliderInput(inputId = "promWidthTrxs", label = "Promoters Length", min = 0, max = 10000, step = 100, value = 1500),
                    shiny::numericInput(inputId = "zoomCoordinatesL", label = "Zoom Coordinate - From", NULL),
                    shiny::numericInput(inputId = "zoomCoordinatesR", label = "Zoom Coordinate - To", NULL),
                    width = 2
                ),

                shiny::mainPanel(
                    shinycssloaders::withSpinner(shiny::plotOutput(outputId = "methylationVisualizationPlot"), type = 2)
                )
            ),

            shiny::tabPanel(
                "Methylation vs Expression",
                icon = shiny::icon("chart-bar"),
                shiny::sidebarPanel(
                    shinyWidgets::actionBttn("correlate", "Correlate", style = 'material-flat', size = 'md', color = 'primary'),                    
                    shiny::hr(style = "border-top: 0px solid white;"),
                    shiny::fileInput(inputId = "expressionFile", label = "Upload Expression Profile", placeholder = "select a file", multiple = FALSE),
                    shinyWidgets::prettyRadioButtons("sepExpr", "Select Delimiter Type", c("Comma" = ",", "Space" = " ", "Tab" = "\t"), selected = ",", shape = "round", inline = TRUE),
                    shiny::hr(),
                    # PARAMETERS 
                    shiny::selectInput(inputId = "paramType", label = "Type of Parameter to filter Overlapping Methylation", choices = c("DMR Length" = "dmr.length", "Overlap Length" = "overlap.length", "Overlap Perc." = "overlap.percentage"), selected = "overlap.length"),
                    shiny::numericInput("overlapParamThr", "Overlapping Methylation Value Threshold for Filtering", value = 100),
                    shiny::numericInput("expressionStatThr", "Statistic Threshold", value = 0.05, step = 0.01),
                    shiny::numericInput("expressionLFCthr", "LogFC Threshold", value = 0.5, step = 0.05),
                    shiny::numericInput("metExprBetaThr", "Beta Diff Threshold", value = .3, step = 0.05),
                    shiny::selectInput(inputId = "metExpressionCorDB", label = "Filter by DB", choices = c("dgv", "gnomad", "NCG", "COSMIC", "cosmic.CGIs", "CGIs"), multiple = TRUE),
                    shiny::selectInput(inputId = "metExpressionCorType", label = "Select Correlation Type", choices = c("pearson", "kendall", "spearman"), selected = "pearson"),
                    shiny::selectInput(inputId = "featExprFilt", label = "Select Path for filtering Genes", choices = c('select a path', '-'), multiple = TRUE),
                    shiny::hr(),
                    shiny::numericInput("expressionColSymbol", "Column Position of Gene ID", value = 1),
                    shiny::numericInput("expressionColStat", "Column Position of Used Statistics", value = 6),
                    shiny::numericInput("expressionColLogFC", "Column position of logFC", value = 2),
                    shinyWidgets::prettyRadioButtons("convertGeneExpression", "Select TRUE if Gene IDs are not Symbols", c("TRUE" = "True", "FALSE" = "False"), selected = "False", shape = "curve"),
                    shiny::selectInput(inputId = "convertGeneExpressionFrom", label = "Select Gene IDs Annotation Type to Translate", choices = c("NONE","ENTREZID", "EXONID", "GENEBIOTYPE", "GENEID", "GENENAME", "PROTDOMID", "PROTEINDOMAINID", "PROTEINDOMAINSOURCE", "PROTEINID", "SEQNAME", "SEQSTRAND", "SYMBOL", "TXBIOTYPE", "TXID", "TXNAME", "UNIPROTID"), selected = "NONE"),
                    width = 2
                ),

                shiny::mainPanel(
                    shinyWidgets::checkboxGroupButtons(
                        inputId = "featureMetExpr",
                        label = "Select Methylated Features to Consider",
                        choices = list("Promoters" = "promoters","Heads" = 'heads', "TSS Surnd.ing" = "tss.surrounding" ,"Five UTRs" = "fiveUTRs", "Exons" = "exons", "First Exons" = "exons1", "Introns" = "introns","First Introns" = "introns1", "Three UTRs" = "threeUTRs"),
                        status = "primary",
                        selected = c("promoters", "heads"),
                        size = "sm",
                        individual = FALSE,
                        justified = TRUE
                    ),
                    shiny::hr(style = "border-top: 1px solid #000000;"),
                    shinyWidgets::dropdownButton(
                        colourpicker::colourInput("colrsMetExpression", "Select Color 1", "lightgrey"),
                        colourpicker::colourInput("lmfitCol1", "Select LM Fit 1 Color", "#ff8800"),
                        colourpicker::colourInput("lmfitCol2", "Select LM Fit 2 Color", "#4b4d4b"),
                        shiny::selectInput(inputId = "metExpressionPlotType", label = "Plot Type", choices = c("simple", "splitted"), selected = "splitted"),
                        shiny::selectInput(inputId = "showText", label = "Show Gene Symbols", choices = c('TRUE'=TRUE, 'FALSE'=FALSE), selected = FALSE),
                        shiny::selectInput(inputId = "metExpressionPal", label = "Plot Palette", choices = hcl.pals(), selected = "RdGy"),
                        circle = TRUE, status = "primary", icon = shiny::icon("cogs"), size = "sm"
                    ),
                    shinycssloaders::withSpinner(shiny::plotOutput(outputId = "metExpressionPlot", width = "auto", height = 600), type = 1),
                    shiny::hr(style = "border-top: 0px solid white;"),
                    shiny::downloadButton("downloadMetExpressionPlot", "Download Plot"),
                    shiny::hr(),
                    shiny::hr(style = "border-top: 0px solid white;"),
                    shinycssloaders::withSpinner(DT::dataTableOutput(outputId = "metExpressionTable"), type = 1),
                    shiny::hr(style = "border-top: 0px solid white;"),
                    shiny::downloadButton("downloadMetExpressionTable", "Download Table"),
                    shiny::hr(style = "border-top: 0px solid white;")

                )
            ),

            shiny::tabPanel(
                "Methylated Enhancers vs Expression",
                icon = shiny::icon("project-diagram"),
                shiny::sidebarPanel(
                    shinyWidgets::actionBttn("correlateEnh", "Correrlate", style = 'material-flat', size = 'md', color = 'primary'),
                    shiny::hr(style = "border-top: 0px solid white;"),
                    shiny::fileInput(inputId = "expressionFile", label = "Upload Expression Profile", placeholder = "select a file", multiple = FALSE),
                    shinyWidgets::prettyRadioButtons("sepExpr", "Select Delimiter Type", c("Comma" = ",", "Space" = " ", "Tab" = "\t"), selected = ",", shape = "round", inline = TRUE),
                    shiny::hr(),
                    # PARAMETERS
                    shiny::selectInput(inputId = "paramTypeEnh", label = "Type of Parameter to filter Overlapping Methylation", choices = c("DMR Length" = "dmr.length", "Overlap Length" = "overlap.length", "Overlap Perc." = "overlap.percentage"), selected = "overlap.percentage"),
                    shiny::numericInput("overlapParamThrEnh", "Overlapping Methylation Value Threshold for Filtering", value = 40),
                    shiny::numericInput("expressionStatThr", "Statistic Threshold", value = 0.05, step = 0.01),
                    shiny::numericInput("expressionLFCthr", "LogFC Threshold", value = 0.5, step = 0.05),
                    shiny::sliderInput(inputId = "enhExprBetaThr", label = "Enhancer Methylation Beta Threshold", min = 0, max = 1, step = 0.05, value = 0.3),
                    shinyWidgets::prettyRadioButtons("enhancerDB", "Select Enhancer DB", c("4DGenome" = "4DGenome", "FANTOM5" = "FANTOM5"), selected = "FANTOM5", shape = "curve"),
                    shiny::selectInput(inputId = "enhExpressionCorType", label = "Select Correlation Type", choices = c("pearson", "kendall", "spearman"), selected = "pearson"),
                    shiny::hr(),
                    shiny::numericInput("expressionColSymbol", "Column Position of Gene ID", value = 1),
                    shiny::numericInput("expressionColStat", "Column Position of Used Statistics", value = 6),
                    shiny::numericInput("expressionColLogFC", "Column position of logFC", value = 2),
                    shinyWidgets::prettyRadioButtons("convertGeneExpression", "Select TRUE if Gene IDs are not Symbols", c("TRUE" = "True", "FALSE" = "False"), selected = "False", shape = "curve"),
                    shiny::selectInput(inputId = "convertGeneExpressionFrom", label = "Select Gene IDs Annotation Type to Translate", choices = c("NONE","ENTREZID", "EXONID", "GENEBIOTYPE", "GENEID", "GENENAME", "PROTDOMID", "PROTEINDOMAINID", "PROTEINDOMAINSOURCE", "PROTEINID", "SEQNAME", "SEQSTRAND", "SYMBOL", "TXBIOTYPE", "TXID", "TXNAME", "UNIPROTID"), selected = "NONE"),
                    width = 2
                ),

                shiny::mainPanel(
                    shinyWidgets::dropdownButton(
                        colourpicker::colourInput("colrsEnhExpression", "Select Color 1", "lightgrey"),
                        colourpicker::colourInput("lmfitCol1Enh", "Select LM Fit 1 Color", "#ff8800"),
                        colourpicker::colourInput("lmfitCol2Enh", "Select LM Fit 2 Color", "#4b4d4b"),
                        shiny::selectInput(inputId = "enhExpressionPlotType", label = "Plot Type", choices = c("simple", "splitted"), selected = "splitted"),
                        shiny::selectInput(inputId = "showTextEnh", label = "Show Gene Symbols", choices = c('TRUE'=TRUE, 'FALSE'=FALSE), selected = FALSE),
                        shiny::selectInput(inputId = "enhExpressionPal", label = "Plot Palette", choices = hcl.pals(), selected = "RdGy"),
                        circle = TRUE, status = "primary", icon = shiny::icon("cogs"), size = "sm"
                    ),
                    shinycssloaders::withSpinner(shiny::plotOutput(outputId = "enhExpressionPlot"), type = 1),

                    shiny::hr(style = "border-top: 0px solid white;"),
                    shiny::downloadButton("downloadEnhExpressionPlot", "Download Plot"),
                    shiny::hr(),
                    shiny::hr(style = "border-top: 0px solid white;"),
                    shinycssloaders::withSpinner(DT::dataTableOutput(outputId = "enhExpressionTable"), type = 1),
                    shiny::downloadButton("downloadEnhExpressionTable", "Download"),

                )
            ),


            shiny::tabPanel(
                icon = shiny::icon("puzzle-piece"),
                "Methylation Enrichment",
                shiny::sidebarPanel(
                    shinyWidgets::actionBttn("startEnr", "Enrichment", style = 'material-flat', size = 'md', color = 'primary'),
                    shiny::hr(style = "border-top: 0px solid white;"),
                    # PARAMETERS
                    shiny::selectInput(inputId = "enrStat", label = "Statistics for Results Filtering", choices = c("P.Value" = "P.value", "Adj.P.Value" = "Adjusted.P.value", "Overlap" = "Overlap"), selected = "P.value"),
                    shiny::selectInput(inputId = "activeFeaturesEnr", label = "Select Annotation from where gene symbols are taken", choices = c("genes","heads","tss.surrounding","promoters","fiveUTRs","exons","introns","threeUTRs"), multiple = T, selected = c("genes")),
                    shiny::selectInput(inputId = "enrMetAnnDB", label = "DBs to Query", choices = c("ClinVar_2019", "OMIM_Disease", "Elsevier_Pathway_Collection", "MSigDB_Hallmark_2020", "MSigDB_Oncogenic_Signatures", "GO_Biological_Process_2018", "Human_Phenotype_Ontology", "KEGG_2016", "NCI-Nature_2016", "Panther_2016", "Reactome_2016", "WikiPathways_2019_Human"), multiple = TRUE),
                    shiny::numericInput("enrThr", "Value Threshold for Filtering by Statistics", value = 0.05, step = 0.001),                    
                    width = 2
                ),

                shiny::mainPanel(
                    shinycssloaders::withSpinner(DT::dataTableOutput(outputId = "enrichrTable"), type = 1),
                    shiny::hr(style = "border-top: 0px solid white;"),
                    shiny::downloadButton("downloadEnrichrTable", "Download Table"),
                    shiny::hr(),
                    shinyWidgets::dropdownButton(
                        shiny::numericInput("numPlotEnr", "Num of Paths to Plot", value = 20, step = 1),
                        shiny::selectInput("enrhPlotType", "Plot Type", choices = c("Barplot" = "barplot", "Lollipop" = "lollipop"), selected = "barplot"),
                        shiny::selectInput("enrhPalCol", "Color Palette for Lollipop", choices = hcl.pals(), selected = "Viridis"),
                        colourpicker::colourInput("colEnrPlotHyper", "Hyper Color Barplot", "#ff0000"),
                        colourpicker::colourInput("colEnrPlotHypo", "Hypo Color Barplot", "#00b3ff"),
                        shiny::textInput("enrThrPlot", "Values for Statistics Thresholds to Plot, use ',' for multiple values", value = "0.01, 0.05"),
                        shiny::selectInput("enrThrPlotCol", "Colors for Statistics Thresholds to Plot", choices = c("red", "yellow", "green", "orange", "magenta", "brown", "black", "grey", "grey20", "grey70"), multiple =TRUE, selected = c("green","orange")),
                        circle = TRUE, status = "primary", icon = shiny::icon("cogs"), size = "sm"
                    ),
                    shinycssloaders::withSpinner(shiny::plotOutput(outputId = "enrPlot", height = "700px"), type = 1),
                    shiny::hr(style = "border-top: 0px solid white;"),
                    shiny::downloadButton("downloadEnrichrPlot", "Download Plot"),
                    shiny::hr(style = "border-top: 0px solid white;")
                )
            ),

            shiny::tabPanel(
                icon = shiny::icon("caret-square-down"),
                "TFs Analysis",

                shiny::sidebarPanel(
                    shinyWidgets::actionBttn("startTF", "Compute", style = 'material-flat', size = 'md', color = 'primary'),
                    shiny::hr(style = "border-top: 0px solid white;"),
                    # PARAMETERS 
                    shiny::fileInput(inputId = "expressionFile", label = "Upload Expression Profile", placeholder = "select a file", multiple = FALSE),
                    shinyWidgets::prettyRadioButtons("sepExpr", "Select Delimiter Type", c("Comma" = ",", "Space" = " ", "Tab" = "\t"), selected = ",", shape = "round", inline = TRUE),
                    shiny::hr(),
                    shiny::h5("Expression File Parameters", style = "color:black"),                    
                    shiny::numericInput("expressionColSymbol", "Column Position of Gene ID", value = 1),
                    shiny::numericInput("expressionColStat", "Column Position of Used Statistics", value = 6),
                    shiny::numericInput("expressionColLogFC", "Column position of logFC", value = 2),
                    shinyWidgets::prettyRadioButtons("convertGeneExpression", "Select TRUE if Gene IDs are not Symbols", c("TRUE" = "True", "FALSE" = "False"), selected = "False", shape = "curve"),
                    shiny::selectInput(inputId = "convertGeneExpressionFrom", label = "Select Gene IDs Annotation Type to Translate", choices = c("NONE","ENTREZID", "EXONID", "GENEBIOTYPE", "GENEID", "GENENAME", "PROTDOMID", "PROTEINDOMAINID", "PROTEINDOMAINSOURCE", "PROTEINID", "SEQNAME", "SEQSTRAND", "SYMBOL", "TXBIOTYPE", "TXID", "TXNAME", "UNIPROTID"), selected = "NONE"),
                    shiny::hr(),
                    shiny::h5("TFs Targets/Expression Correlation Parameters", style = "color:black"),                    
                    shiny::numericInput("expressionLFCthrTF", "LogFC Threshold", value = 0.584, step = 0.05),
                    shiny::numericInput("metExprBetaThrTF", "Beta Diff Threshold", value = .3, step = 0.05),
                    shiny::numericInput("expressionStatThrTF", "Statistic Threshold", value = 0.05, step = 0.01),
                    shiny::selectInput(inputId = "paramTypeTF", label = "Type of Parameter to filter Overlapping Methylation", choices = c("DMR Length" = "dmr.length", "Overlap Length" = "overlap.length", "Overlap Perc." = "overlap.percentage"), selected = "overlap.percentage"),
                    shiny::numericInput("overlapParamThrTF", "Overlapping Methylation Value Threshold for Filtering", value = 30),
                    shiny::hr(),
                    shiny::h5("TFs Targets Enrichment Parameters", style = "color:black"),
                    shiny::selectInput(inputId = "enrMetAnnDBTF", label = "DBs to Query", choices = c("ClinVar_2019", "OMIM_Disease", "Elsevier_Pathway_Collection", "MSigDB_Hallmark_2020", "MSigDB_Oncogenic_Signatures", "GO_Biological_Process_2018", "Human_Phenotype_Ontology", "KEGG_2016", "NCI-Nature_2016", "Panther_2016", "Reactome_2016", "WikiPathways_2019_Human"), multiple = TRUE),
                    shiny::numericInput("enrLFCthrTF", "LogFC Threshold", value = .584, step = 0.01),
                    shiny::selectInput(inputId = "enrStatTFDf", label = "Statistics Type for Results Filtering", choices = c("P.Value" = "P.value", "Adj.P.Value" = "Adjusted.P.value"), selected = "P.value"),
                    shiny::numericInput("enrStatThrTF", "Statistics Threshold for Results Filtering", value = 0.01, step = 0.01),
                    width = 2
                ),                    

                shiny::mainPanel(
                    # association DF
                    shinyWidgets::checkboxGroupButtons(
                        inputId = "featureMetTFExpr",
                        label = "Select Methylated Features to Consider",
                        choices = list("Promoters" = "promoters","Heads" = 'heads', "TSS Surnd.ing" = "tss.surrounding" ,"Five UTRs" = "fiveUTRs", "Exons" = "exons", "First Exons" = "exons1", "Introns" = "introns","First Introns" = "introns1", "Three UTRs" = "threeUTRs"),
                        status = "primary",
                        selected = c("promoters", "heads"),
                        size = "sm",
                        individual = FALSE,
                        justified = TRUE
                    ),

                    shinycssloaders::withSpinner(DT::dataTableOutput(outputId = "tfExprCorOut"), type = 1),
                    shiny::hr(style = "border-top: 0px solid white;"),
                    shiny::downloadButton("downloadTfExprCorOut", "Download Table"),
                    shiny::hr(),
                    
                    # plot TF meth targets expression
                    shinyWidgets::dropdownButton(
                        shiny::textInput("TfSymbol", "Select TF to visualize"),
                        colourpicker::colourInput("colMetTf", "Methylation Bar Color", "#8e0000"),
                        shiny::selectInput(inputId = "palTFexprPlot", label = "Expression Color Palette", choices = hcl.pals(), selected = "Cold"),
                        circle = TRUE, status = "primary", icon = shiny::icon("cogs"), size = "sm"
                    ),
                    shinycssloaders::withSpinner(shiny::plotOutput(outputId = "tfExprsPlotOut", height = 300, width = 350), type = 1),
                    shiny::hr(style = "border-top: 0px solid white;"),
                    shiny::downloadButton("downloadTfExprsPlotOut", "Download Table"),
                    shiny::hr(style = "border-top: 0px solid white;"),
                    shiny::hr(),

                    # TF targets enrichment data.frame
                    shinycssloaders::withSpinner(DT::dataTableOutput(outputId = "enrichrTFsOut"), type = 1),
                    shiny::hr(style = "border-top: 0px solid white;"),
                    shiny::downloadButton("downloadEnrichrTFsOut", "Download Table"),
                    shiny::hr(),



                    # TF targets enrichment plot 
                    shinyWidgets::dropdownButton(
                        shiny::selectInput(inputId = "enrStatTF", label = "Statistics for Results Filtering", choices = c("P.Value" = "P.value", "Adj.P.Value" = "Adjusted.P.value"), selected = "P.value"),
                        shiny::numericInput("numPlotEnrTF", "Num of Paths to Plot", value = 20, step = 1),
                        shiny::selectInput("enrhPalColTF", "Color Palette for Lollipop", choices = hcl.pals(), selected = "Viridis"),
                        shiny::textInput("enrThrPlotTF", "Values for Statistics Thresholds to Plot, use ',' for multiple values", value = "0.01, 0.05"),
                        shiny::selectInput("enrThrPlotColTF", "Colors for Statistics Thresholds to Plot", choices = c("red", "yellow", "green", "orange", "magenta", "brown", "black", "grey", "grey20", "grey70"), multiple =TRUE, selected = c("green","orange")),
                        circle = TRUE, status = "primary", icon = shiny::icon("cogs"), size = "sm"
                    ),
                    
                    shinycssloaders::withSpinner(shiny::plotOutput(outputId = "enrPlotTF", height = "500px"), type = 1),
                    shiny::hr(style = "border-top: 0px solid white;"),
                    shiny::downloadButton("downloadEnrichrPlotTf", "Download Table"),
                    shiny::hr(),
                )
            )
        )
    )


    #####################################################################
    #######################         SERVER          #####################
    #####################################################################

    server <- function(input = input, output = output, server = server) {
        libpath <- paste0(.libPaths()[1], "/methyl.O")

        output$logo <- shiny::renderImage(
            {
                logo <- paste0(.libPaths()[1], "/methyl.O/data/logo6.png")
                list(
                    src = logo,
                    width = "40%",
                    contentType = "image/png"
                )
            },
            deleteFile = F
        )

        output$map <- shiny::renderImage(
            {
                map <- paste0(.libPaths()[1], "/methyl.O/data/methyloMap.png")
                list(
                    src = map,
                    width = "100%",
                    contentType = "image/png"
                )
            },
            deleteFile = F
        )




        tableIn <- shiny::reactive({
            if (!is.null(input$bedfile)) {
                inTable <- read.csv(input$bedfile$datapath, header = input$header, sep = input$sep, stringsAsFactors = F, quote = "")
                if(!input$fraction1){
                    inTable[,4] <- as.numeric(as.character(inTable[,4]))/100
                }
                return(inTable)
            } else if (is.null(input$bedfile)) {
                data("DMRsSubset")
                DMRsSubset$seqnames <- as.character(DMRsSubset$seqnames)
                return(DMRsSubset)
            }
        })



        checkDefault <- paste0(1500, 1500, "TRUE", "ensembl", "hg19", "TRUE", 0.3, 0.3, 4)
        checkCurent <- shiny::reactive({
            paste0(input$promLength,input$headLength, input$longestTrxs, input$annotation, input$hg, input$annotationFast, input$thrBeta, input$thrCGIs, input$colBetaDiff)
        })

        resultsAnnot <- shiny::eventReactive(input$annotate, {
            if (!is.null(input$bedfile)) {
                methyl.O::annotateDMRs(tableIn(), prom.length = input$promLength, head.length = input$headLength, longest.trx = input$longestTrxs, annotation = input$annotation, hg = input$hg, annotation.fast = input$annotationFast, thr.beta = input$thrBeta, thr.cgis = input$thrCGIs, col.betadiff = input$colBetaDiff, col.beta1 = input$betacol1, col.beta2 = input$betacol2)
            } else {
                if (checkDefault != checkCurent()) {
                    methyl.O::annotateDMRs(tableIn(), prom.length = input$promLength, head.length = input$headLength, longest.trx = input$longestTrxs, annotation = input$annotation, hg = input$hg, annotation.fast = input$annotationFast, thr.beta = input$thrBeta, thr.cgis = input$thrCGIs, col.betadiff = input$colBetaDiff, col.beta1 = 5, col.beta2 = 6)
                } else {
                    readRDS(paste0(libpath, "/data/exampleAnnoTable.RDS"))
                }
            }
        })

        checkAnnot <- shiny::reactive ({
            if(!is.null(input$bedfile)){
                date()
            } else {
                if (checkDefault != checkCurent()) {
                    date()
                } else {
                    'default'
                }
            }
        })


        resultsScore <- shiny::eventReactive(input$annotate,{
            methyl.O::scoreAnnotatedDMRs(resultsAnnot(), score.modifier = input$scoreMinMax, active.features = input$scoreFeatureSelected)
        })


        resultsTab <- shiny::reactive({
            indTab <- match(input$choosenFeature, names(resultsScore()))
            finalSV <- resultsScore()[[indTab]]

            finalWmInd <- which(finalSV[, "width"] >= as.numeric(input$filterTabWm))
            finalWMInd <- which(finalSV[, "width"] <= as.numeric(input$filterTabWM))

            finalPmInd <- which(finalSV[, grep("perc", colnames(finalSV))] >= as.numeric(input$filterTabPm))
            finalPMInd <- which(finalSV[, grep("perc", colnames(finalSV))] <= as.numeric(input$filterTabPM))

            hold <- intersect(
                finalWmInd,
                intersect(
                    finalWMInd,
                    intersect(finalPmInd, finalPMInd)
                )
            )

            if (any(input$choosenFeature %in% c("introns", "exons"))) {
                finalRmInd <- which(finalSV[, grep("rank", colnames(finalSV))] >= as.numeric(input$filterTabRm))
                finalRMInd <- which(finalSV[, grep("rank", colnames(finalSV))] <= as.numeric(input$filterTabRM))
                hold <- intersect(hold, intersect(finalRmInd, finalRMInd))
            }

            return(finalSV[hold, ])
        })

        output$resultsTab <- DT::renderDataTable({
            if (input$choosenFeature == "genes") {
                DT::datatable(resultsTab()[order(resultsTab()$score, decreasing = T), c(1:5, 21, 25, 6:20,22:24)], rownames = FALSE, options = list(
                    autoWidth = TRUE, scrollX = TRUE
                ))
            } else {
                DT::datatable(resultsTab(), rownames = FALSE, options = list(
                    autoWidth = TRUE, scrollX = TRUE
                ))
            }
        })

        overviewPlot <- shiny::reactive({
            methyl.O::plotMethylationOverview(resultsScore(), plot.type = input$methylationOverview, palette = input$pals)
        })

        output$overviewPlot <- shiny::renderPlot({
            overviewPlot()
        })

        output$downloadOverviewPlot <- shiny::downloadHandler(
            filename = function() {
                paste0("BetaOverviewPlot ", Sys.Date(), ".pdf")
            },
            content = function(file) {
                device <- function() grDevices::pdf(res = 300)
                ggplot2::ggsave(file, plot = overviewPlot())
            }
        )



        resultsPlotChr <- shiny::reactive({
            resTab <- c()
            for (i in 1:length(resultsScore())) {
                indTab <- match("seqnames", colnames(resultsScore()[[i]]))
                resTab <- c(resTab, as.character(resultsScore()[[i]][[indTab]]))
            }
            resTab <- as.factor(resTab)
            resVec <- sort(summary(resTab), decreasing = T)
            color <- hcl.colors(22, palette = input$col1)
            barplot(resVec, las = 2, col = color, main = "Chromosomes Distribution of Methylated Regions")
        })

        output$resultsPlotChr <- shiny::renderPlot({
            resultsPlotChr()
        })

        output$downloadResultsPlotChr <- shiny::downloadHandler(
            filename = function() {
                paste0("LenghtsOverviewPlot ", Sys.Date(), ".pdf")
            },
            content = function(file) {
                pdf(file)
                    resTab <- c()
                    for (i in 1:length(resultsScore())) {
                        indTab <- match("seqnames", colnames(resultsScore()[[i]]))
                        resTab <- c(resTab, as.character(resultsScore()[[i]][[indTab]]))
                    }
                    resTab <- as.factor(resTab)
                    resVec <- sort(summary(resTab), decreasing = T)
                    color <- hcl.colors(22, palette = input$col1)
                    barplot(resVec, las = 2, col = color, main = "Chromosomes Distribution of Methylated Regions")
                dev.off()
            }
        )

        observeEvent(input$choosenFeatureHist == "genes" && input$choosenValueSWRHist == "rank", {
            shinyalert::shinyalert("Oppss!, There are not ranks in genes! Please select a feature", imageUrl = "https://www.google.com/url?sa=i&url=https%3A%2F%2Fmassivelolz.com%2Fcat-poker-face%2F&psig=AOvVaw0HNLoZE02AK3CfuEFyYNgC&ust=1612609770265000&source=images&cd=vfe&ved=0CAIQjRxqFwoTCIDR7bvO0u4CFQAAAAAdAAAAABAD", type = "error")
        })


        resultsPlotSWR <- shiny::reactive({
            tmpTab <- resultsScore()

            if (input$choosenFeatureHist == "exons1") {
                tmpTab$exons1 <- tmpTab$exons[tmpTab$exons$rank == 1, ]
            }

            if (input$choosenFeatureHist == "introns1") {
                tmpTab$introns1 <- tmpTab$introns[tmpTab$introns$intron.rank == 1, ]
            }

            indTab <- match(input$choosenFeatureHist, names(tmpTab))
            resTab <- tmpTab[[indTab]]
            indPlot <- match(input$choosenValueSWRHist, colnames(resTab))
            resPlot <- resTab[, indPlot]
            color <- hcl.colors(length(resPlot), palette = input$col2)
            hist.param <- hist(
                main = paste0("Distribution of ", input$choosenValueSWRHist),
                resPlot, col = color, breaks = as.numeric(input$histoBreaks), xlim = c(0, max(resPlot) + (max(resPlot) * .1)),
                xlab = input$choosenValueSWRHist
            )
        })

        output$resultsPlotSWR <- shiny::renderPlot ({
            resultsPlotSWR()
        }) 
        


        output$downloadResultsPlotSWR <- shiny::downloadHandler(
            filename = function() {
                paste0("FeatureOverview", Sys.Date(), ".pdf")
            },
            content = function(file) {
                pdf(file)
                    indTab <- match(input$choosenFeatureHist, names(resultsScore()))
                    resTab <- resultsScore()[[indTab]]
                    indPlot <- match(input$choosenValueSWRHist, colnames(resTab))
                    resPlot <- resTab[, indPlot]
                    color <- hcl.colors(length(resPlot), palette = input$col2)
                    hist.param <- hist(
                        main = paste0("Distribution of ", input$choosenValueSWRHist),
                        resPlot, col = color, breaks = as.numeric(input$histoBreaks), xlim = c(0, max(resPlot) + (max(resPlot) * .1)),
                        xlab = input$choosenValueSWRHist
                    )
                dev.off()
            }
        )


        resultsPie <- shiny::reactive({
            tmpTab <- resultsScore()
            selected <- c("promoters",  "fiveUTRs", "exons", "introns", "Three UTRs" = "threeUTRs")

            if (any(input$choosenFeaturePie == "exons1")) {
                tmpTab$exons1 <- tmpTab$exons[tmpTab$exons$rank == 1, ]
                selected <- c(selected, "exons1")
            }

            if (any(input$choosenFeaturePie == "introns1")) {
                tmpTab$introns1 <- tmpTab$introns[tmpTab$introns$intron.rank == 1, ]
                selected <- c(selected, "introns1")
            }

            if (any(input$choosenFeaturePie == "heads")) {
                selected <- c(selected, "heads")
            }

            if (any(input$choosenFeaturePie == "tss.surrounding")) {
                selected <- c(selected, "tss.surrounding")
            }


            pieValues <- unlist(lapply(tmpTab[selected], nrow))
            pieNames <- paste(names(pieValues), "\n", pieValues)
            if(input$featureDistType == 'Pie'){
                pie(pieValues, main = "Methylated Features", labels = pieNames, col = hcl.colors(length(pieNames), palette = input$piePalette))
            } else {
                barplot(pieValues, main = "Methylated Features", col = hcl.colors(length(pieNames), palette = input$piePalette), ylim = c(0, max(pieValues) + (max(pieValues)*.1)))
            }
        })

        output$resultsPie <- shiny::renderPlot({
            resultsPie()
        })


        output$downloadResultsPie <- shiny::downloadHandler(
            filename = function() {
                paste0("MethylatedFaturesSummary", Sys.Date(), ".pdf")
            },
            content = function(file) {
                pdf(file)
                    tmpTab <- resultsScore()
                     selected <- c("promoters", "fiveUTRs", "exons", "introns", "Three UTRs" = "threeUTRs")

                    if (any(input$choosenFeaturePie == "exons1")) {
                        tmpTab$exons1 <- tmpTab$exons[tmpTab$exons$rank == 1, ]
                        selected <- c(selected, "exons1")
                    }

                    if (any(input$choosenFeaturePie == "introns1")) {
                        tmpTab$introns1 <- tmpTab$introns[tmpTab$introns$intron.rank == 1, ]
                        selected <- c(selected, "introns1")
                    }

                    if (any(input$choosenFeaturePie == "heads")) {
                       selected <- c(selected, "heads")
                    }

                    pieValues <- unlist(lapply(tmpTab[selected], nrow))
                    pieNames <- paste(names(pieValues), "\n", pieValues)
                    if (input$featureDistType == "Pie") {
                        pie(pieValues, main = "Methylated Features", labels = pieNames, col = hcl.colors(length(pieNames), palette = input$piePalette))
                    } else {
                        barplot(pieValues, main = "Methylated Features", col = hcl.colors(length(pieNames), palette = input$piePalette), width = 0.5, ylim = c(0, max(pieValues) + (max(pieValues)*.1)))
                    }
                dev.off()
            }
        )





        output$downloadData <- shiny::downloadHandler(
            filename = function() {
                paste0("methyl.O ", Sys.Date(), ".csv")
            },
            content = function(file) {
                write.csv(resultsTab(), file)
            }
        )


        shinyBS::addTooltip(session = getDefaultReactiveDomain(), id = "scoreMinMax", title = "0 to focus the score on Db annotations. 1 to focus the score on methylation affecting genes expression", placement = "right", trigger = "hover")

        ### ENHANCER ANNOTATION
        checkDefaultEnh <- paste0("hg19",0.3,40, "overlap.percentage")
        checkCurentEnh <- shiny::reactive({
            paste0(input$hg,input$thr.beta.enhancer, input$overlapParamThrEnhAnnot, input$paramTypeEnhAnnot)
        })

        enhancerAnn <- shiny::eventReactive(input$annotateEnh, {
            if (!is.null(input$bedfile)) {
                methyl.O::annotateEnhancers(tableIn(), hg = input$hg, thr.beta = input$thr.beta.enhancer, overlap.param.thr = input$overlapParamThrEnhAnnot, param.type = input$paramTypeEnhAnnot, score.modifier = input$scoreModifierenh, col.betadiff = input$colBetaDiff)
            } else {
                defaultTableEnh <- readRDS(paste0(libpath, "/data/enhancerExample.RDS"))
                if (checkDefaultEnh != checkCurentEnh()) {
                    methyl.O::annotateEnhancers(tableIn(), hg = input$hg, thr.beta = input$thr.beta.enhancer, overlap.param.thr = input$overlapParamThrEnhAnnot, param.type = input$paramTypeEnhAnnot, score.modifier = input$scoreModifierenh, col.betadiff = input$colBetaDiff)
                } else {
                    defaultTableEnh
                }
            }
        })

		selected.param <- shiny::reactive ({
				if (input$paramTypeEnhAnnot == "overlap.percentage") {
					"enhancer.perc"
				} else if (input$paramTypeEnhAnnot == "overlap.length") {
					"overlap.width"
				} else if (input$paramTypeEnhAnnot == "dmr.length") {
					"width"
				}
		})


        output$resultsTabEnhancer <- DT::renderDataTable({
            DT::datatable(enhancerAnn()[order(enhancerAnn()[,selected.param()], decreasing = T), ], rownames = FALSE, options = list(autoWidth = TRUE, scrollX = TRUE))
        })

        distEnhancer <- shiny::reactive({
            if (input$plotDistEnhancer == "density") {
                densEnh <- density(enhancerAnn()[,selected.param()])
                plot(densEnh, main = "Methylation Parameter Distribution", xlab = "Distribution", xlim=c(0,(max(enhancerAnn()[,selected.param()])*.05) + max(enhancerAnn()[,selected.param()])))
                polygon(densEnh, col = input$colDistEnhancer)
            } else {
                hist(enhancerAnn()[,selected.param()], xlim=c(0,(max(enhancerAnn()[,selected.param()])*.05) + max(enhancerAnn()[,selected.param()])), col = input$colDistEnhancer, xlab = "Distribution", main = "Methylation Parameter Distribution", input$binDist)
            }
        })

        output$distEnhancer <- shiny::renderPlot({
            distEnhancer()
        })
        

        output$downloadDistEnhancer <- shiny::downloadHandler(
            filename = function() {
                paste0("EnhancerMethylationOverlapPercentageOverview", Sys.Date(), ".pdf")
            },
            content = function(file) {
                pdf(file)
                    if (input$plotDistEnhancer == "density") {
                        densEnh <- density(enhancerAnn()[,selected.param()])
                        plot(densEnh, main = "Methylation Overlap Percentage Distribution", xlab = "Distribution")
                        polygon(densEnh, col = input$colDistEnhancer)
                    } else {
                        hist(enhancerAnn()[,selected.param()], col = input$colDistEnhancer, xlab = "Distribution")
                    }
                dev.off()
            }
        )


        distBeta <- shiny::reactive({
            if (input$plotDistBeta == "density") {
                densBeta <- density(enhancerAnn()$beta)
                plot(densBeta, main = "Beta Difference Distribution")
                polygon(densBeta, col = input$colDistBeta, xlim = c(-1,1))
            } else {
                hist(enhancerAnn()$beta, col = input$colDistBeta, main = "Beta Difference Distribution",  xlim = c(-1,1))
            }
        })

        output$distBeta <- shiny::renderPlot({
            distBeta()
        })

        output$downloadDistBeta <- shiny::downloadHandler(
            filename = function() {
                paste0("EnhancerMethylationBetaOverview", Sys.Date(), ".pdf")
            },
            content = function(file) {
                pdf(file)
                if (input$plotDistBeta == "density") {
                    densBeta <- density(enhancerAnn()$beta)
                    plot(densBeta, main = "Beta Difference Distribution")
                    polygon(densBeta, col = input$colDistBeta, xlim = "Distribution")
                } else {
                    hist(enhancerAnn()$beta, col = input$colDistBeta, xlim = "Distribution")
                }
                dev.off()
            }
        )



        output$downloadDataEnhancer <- shiny::downloadHandler(
            filename = function() {
                paste0("enhancerAnnotation ", Sys.Date(), ".csv")
            },
            content = function(file) {
                write.csv(enhancerAnn(), file)
            }
        )


        # METHYLATION VISUALIZATION

        # adapt inputs to function
        betaName1 <- shiny::reactive({
            if (input$betaName1 == "") {
                out <- NULL
                out
            } else {
                input$betaName1
            }
        })

        betaName2 <- shiny::reactive({
            if (input$betaName2 == "") {
                out <- NULL
                out
            } else {
                input$betaName2
            }
        })

        observeEvent(input$symbolVisualize, {
            updateTextInput(session = getDefaultReactiveDomain(), inputId = "zoomCoordinatesL", value = as.numeric(resultsScore()[[1]][match(toupper(input$symbolVisualize), resultsScore()[[1]]$symbol), "gene.start"]))
            updateTextInput(session = getDefaultReactiveDomain(), inputId = "zoomCoordinatesR", value = as.numeric(resultsScore()[[1]][match(toupper(input$symbolVisualize), resultsScore()[[1]]$symbol), "gene.end"]))
        })

        methylationVisualizationPlot <- shiny::eventReactive(input$visualizeMeth, {
            methyl.O::plotDMRs(resultsScore(),
                input$symbolVisualize,
                annotation = input$annotation,
                hg = input$hg,
                beta1.name = betaName1(),
                beta2.name = betaName2(),
                beta.colors = c(input$colbeta1, input$colbeta2),
                blackandwhite = input$blackwhite,
                show.all.transcripts = input$showAllTrxs,
                prom.width = input$promWidthTrxs,
                coord.zoom = c(input$zoomCoordinatesL, input$zoomCoordinatesR),
                smartzoom = input$smartZoomTrxs
            )
        })

        output$methylationVisualizationPlot <- shiny::renderPlot({ 
                   methylationVisualizationPlot()
        })



        ### EXPRESSION METHYLATION
        expressionProfile <- shiny::eventReactive(input$correlate, {
            if (!is.null(input$expressionFile)) {
                read.table(input$expressionFile$datapath, header = T, stringsAsFactors = F, sep = input$sepExpr)
            } else if (is.null(input$bedfile) & is.null(input$expressionFile)) {
                data("expressionSubset")
                expressionSubset
            } else if(!is.null(input$bedfile) & is.null(input$expressionFile)) { 
                NULL
            }
        })


        if(!is.null(expressionProfile)){

            metExpressionPlot <- shiny::reactive ({
                annotatedDMRs2Exprs(
                    annotatedDMRs = resultsScore(), active.features = input$featureMetExpr,
                    expressionProfile = expressionProfile(), col.genes = input$expressionColSymbol,
                    col.stat = input$expressionColStat, stat.thr = input$expressionStatThr, col.logFC = input$expressionColLogFC, logfc.thr = input$expressionLFCthr, convert.genes = input$convertGeneExpression, convert.from = input$convertGeneExpressionFrom, beta.thr = input$metExprBetaThr,
                    overlap.param.thr = input$overlapParamThr, param.type = input$paramType, line.col = input$colrsMetExpression, lmfit.col1 = input$lmfitCol1, lmfit.col2 = input$lmfitCol2, pal = input$metExpressionPal,
                    plot.type = input$metExpressionPlotType, show.text = input$showText, filter.by.genes = featExprFilt.genes(), cor.type = input$metExpressionCorType, return.table = FALSE
                )
            })

            output$metExpressionPlot <- shiny::renderPlot({
                metExpressionPlot()
            })


            output$metExpressionTable <- DT::renderDataTable({
                methyl.O::annotatedDMRs2Exprs(
                    annotatedDMRs = resultsScore(), active.features = input$featureMetExpr, expressionProfile = expressionProfile(),
                    col.genes = input$expressionColSymbol,
                    input$expressionColStat, input$expressionStatThr, input$expressionColLogFC, input$expressionLFCthr, input$convertGeneExpression, input$convertGeneExpressionFrom, input$metExprBetaThr,
                    input$overlapParamThr, param.type = input$paramType, input$colrsMetExpression, input$lmfitCol1, input$lmfitCol2, input$metExpressionPal,
                    input$metExpressionPlotType, show.text = input$showText, input$metExpressionCorType, return.table = TRUE
                )
            })



            output$downloadMetExpressionPlot <- shiny::downloadHandler(
                filename = function() {
                    paste0("MethylationVsExpression", Sys.Date(), ".pdf")
                },
                content = function(file) {
                    pdf(file, width = 13)
                    methyl.O::annotatedDMRs2Exprs(
                        resultsScore(), active.features = input$featureMetExpr,
                        expressionProfile = expressionProfile(), col.genes = input$expressionColSymbol,
                        col.stat = input$expressionColStat, stat.thr = input$expressionStatThr, col.logFC = input$expressionColLogFC, logfc.thr = input$expressionLFCthr, convert.genes = input$convertGeneExpression, convert.from = input$convertGeneExpressionFrom, beta.thr = input$metExprBetaThr,
                        overlap.param.thr = input$overlapParamThr, param.type = input$paramType, line.col = input$colrsMetExpression, lmfit.col1 = input$lmfitCol1, lmfit.col2 = input$lmfitCol2, pal = input$metExpressionPal,
                        plot.type = input$metExpressionPlotType,show.text = input$showText, filter.by.genes = featExprFilt.genes(), cor.type = input$metExpressionCorType, return.table = FALSE
                    )
                    dev.off()
                }
            )




            output$downloadMetExpressionTable <- shiny::downloadHandler(
                filename = function() {
                    paste0("MethylationVsExpression ", Sys.Date(), ".csv")
                },
                content = function(file) {
                    write.csv(metExpressionTable(), file)
                }
            )

        } else {
            shiny::h1("Please upload expression file")
        }



        # ENHANCER METHYLATION EXPRESSION

        enhExpressionPlot <- shiny::eventReactive(input$correlateEnh, { 
            methyl.O::annotatedEnh2Exprs(
                annotatedEnhancers = enhancerAnn(), expressionProfile = expressionProfile(), hg = input$hg,
                enhancer.db = input$enhancerDB, col.genes = input$expressionColSymbol, col.stat = input$expressionColStat,
                stat.thr = input$expressionStatThr, col.logFC = input$expressionColLogFC, logfc.thr = input$expressionLFCthr, convert.genes = input$convertGeneExpression, convert.from = input$convertGeneExpressionFrom, beta.thr = input$enhExprBetaThr,
                overlap.param.thr = input$overlapParamThrEnh, param.type = input$paramTypeEnh, line.col = input$colrsEnhExpression, lmfit.col1 = input$lmfitCol1Enh, lmfit.col2 = input$lmfitCol2Enh, 
                pal = input$enhExpressionPal, plot.type = input$enhExpressionPlotType, show.text = input$showTextEnh, cor.type = input$enhExpressionCorType, return.table = FALSE
            )
        })


        output$enhExpressionPlot  <- shiny::renderPlot({
            enhExpressionPlot()
        })


        output$downloadEnhExpressionPlot <- shiny::downloadHandler(
            filename = function() {
                paste0("EnahncersMethylationVsExpression ", Sys.Date(), ".pdf")
            },
            content <- function(file) {
                pdf(file, width = 13)
                methyl.O::annotatedEnh2Exprs(
                    annotatedEnhancers = enhancerAnn(), expressionProfile = expressionProfile(), hg = input$hg,
                    enhancer.db = input$enhancerDB, col.genes = input$expressionColSymbol, col.stat = input$expressionColStat,
                    stat.thr = input$expressionStatThr, col.logFC = input$expressionColLogFC, logfc.thr = input$expressionLFCthr, convert.genes = input$convertGeneExpression, input$convertGeneExpressionFrom, beta.thr = input$enhExprBetaThr,
                    overlap.param.thr = input$overlapParamThrEnh, param.type = input$paramTypeEnh,  line.col = input$colrsEnhExpression, lmfit.col1 = input$lmfitCol1Enh, lmfit.col2 = input$lmfitCol2Enh, 
                    pal = input$enhExpressionPal, plot.type = input$enhExpressionPlotType, show.text = input$showTextEnh, cor.type = input$enhExpressionCorType, return.table = FALSE
                )
                dev.off()
            }
        )



        enhExpressionTable <- shiny::eventReactive(input$correlateEnh, {
            methyl.O::annotatedEnh2Exprs(
                annotatedEnhancers = enhancerAnn(), expressionProfile = expressionProfile(), hg = input$hg,
                enhancer.db = input$enhancerDB, col.genes = input$expressionColSymbol, col.stat = input$expressionColStat,
                stat.thr = input$expressionStatThr, col.logFC = input$expressionColLogFC, logfc.thr = input$expressionLFCthr, convert.genes = input$convertGeneExpression, input$convertGeneExpressionFrom, beta.thr = input$enhExprBetaThr,
                overlap.param.thr = input$overlapParamThrEnh, param.type = input$paramTypeEnh, line.col = input$colrsEnhExpression, lmfit.col1 = input$lmfitCol1Enh, lmfit.col2 = input$lmfitCol2Enh, 
                pal = input$enhExpressionPal, plot.type = input$enhExpressionPlotType, show.text = input$showTextEnh, cor.type = input$enhExpressionCorType, return.table = TRUE
            )
        })


        output$enhExpressionTable <- DT::renderDataTable({
            DT::datatable(enhExpressionTable(), rownames = FALSE)
        })



        output$downloadEnhExpressionTable <- shiny::downloadHandler(
            filename = function() {
                paste0("EnhancerVsExpression", Sys.Date(), ".csv")
            },
            content = function(file) {
                write.csv(enhExpressionTable(), file)
            }
        )


        ####### ENRICHR

        enrMetAnnDB <- c("ClinVar_2019", "OMIM_Disease", "Elsevier_Pathway_Collection", "MSigDB_Hallmark_2020", "MSigDB_Oncogenic_Signatures", "GO_Biological_Process_2018", "Human_Phenotype_Ontology", "KEGG_2016", "NCI-Nature_2016", "Panther_2016", "Reactome_2016", "WikiPathways_2019_Human")
        enrMetAnnDBUpdated <- shiny::reactive({
            if(length(input$enrMetAnnDB)==0){
                enrMetAnnDB
            } else {
                input$enrMetAnnDB
            }
        })

        checkDefaultEnr <- paste('P.value', 0.01, enrMetAnnDB, 'default', "promoters", "heads", sep ='_')
        checkCurentEnr <- shiny::reactive({
            gsub(' +', '', paste0(input$enrStat,input$enrThr, enrMetAnnDBUpdated(), checkAnnot(), input$activeFeaturesEnr, collapse = '_'))
        })

        # analysis and parsing
        enirchr.results.table.complete <- shiny::eventReactive(input$startEnr, {
            if (!is.null(input$bedfile)) {
                methyl.O::annotatedDMRs2Enrichr(resultsAnnot(), active.features = input$activeFeaturesEnr,  stat.filter = input$enrStat, stat.thr = input$enrThr, db = enrMetAnnDBUpdated())
            } else if (is.null(input$bedfile)) {
                if(checkDefaultEnr==checkCurentEnr()){
                    readRDS(paste0(libpath, "/data/enrichrExample.RDS"))
                } else {
                    methyl.O::annotatedDMRs2Enrichr(resultsAnnot(), active.features = input$activeFeaturesEnr, stat.filter = input$enrStat, stat.thr = input$enrThr, db = enrMetAnnDBUpdated())
                }
            }
        })


        output$enrichrTable <- DT::renderDataTable({
            DT::datatable(enirchr.results.table.complete(), rownames = FALSE, options = list(autoWidth = TRUE, scrollX = TRUE))
        })

        # extract genes from enrichment for filtering expression methylation correlation
        observeEvent(is.null(input$featExprFilt), {
            updateSelectInput(session = getDefaultReactiveDomain(), inputId = "featExprFilt", choices = unique(enirchr.results.table.complete()$Term), selected = input$featExprFilt )
        })

        genes.path <- shiny::reactive({
            if(!is.null(input$featExprFilt)){
                genes.path = c()
                for(cur.path in input$featExprFilt){
                    genes.path <- c(genes.path,unlist(strsplit(enirchr.results.table.complete()$Genes[enirchr.results.table.complete()$Term %in% cur.path],';')))
                }
                genes.path <- genes.path[!is.na(genes.path)]
                genes.path
            } else {
                genes.path = NULL
                genes.path
            }
        })

        # extract genes from results$DB
        genes.db <- shiny::reactive({
            if (!is.null(input$metExpressionCorDB)) {
                genes.db = c()
                for(cur.db in input$metExpressionCorDB){
                    genes.db <- c(genes.db, results[[1]][results[[1]][, cur.db] == 1, "symbol"])
                }
            genes.db

            } else {
                genes.db = NULL
                genes.db
            }
        })

        featExprFilt.genes <- shiny::reactive({
            c(genes.db(), genes.path())
        })

        output$downloadEnrichrTable <- shiny::downloadHandler(
            filename = function() {
                paste0("MethylationEnrichR", Sys.Date(), ".csv")
            },
            content = function(file) {
                write.csv(enirchr.results.table.complete(), file)
            }
        )

        enrPlot <- shiny::reactive({
            methyl.O::plotDMRs2Enrichr(enirchr.results.table.complete(), resultsAnnot(), stat = input$enrStat, n = input$numPlotEnr, plot.type = input$enrhPlotType, pal.col = input$enrhPalCol, col.hyper = input$colEnrPlotHyper, col.hypo = input$colEnrPlotHypo, thrs = as.numeric(unlist(strsplit(input$enrThrPlot,','))), thrs.cols = input$enrThrPlotCol)
        })

        output$enrPlot <- shiny::renderPlot({
            enrPlot()
        })


        output$downloadEnrichrPlot <- downloadHandler(
            filename = function() {
                paste0("EnrichrPlot ", Sys.Date(), ".pdf")
            },
            content <- function(file) {
                pdf(file, width = 13)
                    methyl.O::plotDMRs2Enrichr(enirchr.results.table.complete(), resultsAnnot(), stat = input$enrStat, n = input$numPlotEnr, plot.type = input$enrhPlotType, pal.col = input$enrhPalCol, col.hyper = input$colEnrPlotHyper, col.hypo = input$colEnrPlotHypo, thrs = as.numeric(unlist(strsplit(input$enrThrPlot, ","))), thrs.cols = input$enrThrPlotCol)
                dev.off()
            }
        )


        ########## TF
        # TFs targets expression correlation
        tfExprCorOut <- shiny::eventReactive(input$startTF,{
            associateTFs2Exprs(resultsAnnot(), active.features = input$featureMetTFExpr, expressionProfile = expressionProfile(), col.genes = input$expressionColSymbol, col.stat = input$expressionColStat, stat.thr = input$expressionStatThrTF, 
            col.logFC = input$expressionColLogFC, logfc.thr = input$expressionLFCthrTF, 
            convert.genes = input$convertGeneExpression, convert.frominput$convertGeneExpressionFrom,  beta.thr = input$metExprBetaThrTF,
            overlap.param.thr = input$overlapParamThrTF, param.type = input$paramTypeTF)
        })

        output$tfExprCorOut <- DT::renderDataTable({
            DT::datatable(tfExprCorOut(), rownames = FALSE)
        })

        output$downloadTfExprCorOut <- shiny::downloadHandler(
            filename = function() {
                paste0("TFsExpressionCorrelation", Sys.Date(), ".csv")
            },
            content = function(file) {
                write.csv(tfExprCorOut(), file)
            }
        )


        # TFs meth expression plot
        tfExprsPlotOut <- shiny::eventReactive(input$startTF,{
            methyl.O::plotTFs2Exprs(tfExprCorOut(), input$TfSymbol, col.meth = input$colMetTf , pals.bars = input$palTFexprPlot)
        })

        output$tfExprsPlotOut <- shiny::renderPlot({
            if (input$TfSymbol != "") {
                tfExprsPlotOut()
            } else {
                plot.new()
                text(x = .5, y = .5, cex = 1, col = "Blue", "Waiting For Gene Symbol")
            }
        })


        output$downloadTfExprsPlotOut <- downloadHandler(
            filename = function() {
                paste0("EnrichrTFsPlot ", Sys.Date(), ".pdf")
            },
            content <- function(file) {
                pdf(file, width = 13)
                    methyl.O::plotTFs2Exprs(tfExprCorOut(), input$TfSymbol, col.meth = input$colMetTf , pals.bars = input$palTFexprPlot)
                dev.off()
            }
        )

        # enrichment df
        enrichrTFsOut <- shiny::eventReactive(input$startTF,{
            methyl.O::tfs2Enrichr(tfExprCorOut(), logfc.thr = input$enrLFCthrTF, stat.filter = input$enrStatTFDf, stat.thr = input$enrStatThrTF, db = input$enrMetAnnDBTF)
        })

        output$enrichrTFsOut <- DT::renderDataTable({
            DT::datatable(enrichrTFsOut(), rownames = FALSE, )
        })

        output$downloadEnrichrTFsOut <- shiny::downloadHandler(
            filename = function() {
                paste0("EnrichmentTFsTargets", Sys.Date(), ".csv")
            },
            content = function(file) {
                write.csv(enrichrTFsOut(), file)
            }
        )

        #enrichment plot
        enrPlotTF <- shiny::eventReactive(input$startTF,{
            methyl.O::plotDMRs2Enrichr(enrichrTFsOut(), resultsAnnot(), stat = input$enrStatTF, n = input$numPlotEnrTF, plot.type = 'lollipop', pal.col = input$enrhPalColTF, col.hyper = NULL, col.hypo = NULL, thrs = as.numeric(unlist(strsplit(input$enrThrPlotTF, ","))), thrs.cols = input$enrThrPlotColTF)
        })

        output$enrPlotTF <- shiny::renderPlot({
            enrPlotTF()
        })


        output$downloadEnrichrPlotTf <- downloadHandler(
            filename = function() {
                paste0("EnrichrTfTargets", Sys.Date(), ".pdf")
            },
            content <- function(file) {
                pdf(file, width = 13)
                    methyl.O::plotDMRs2Enrichr(enrichrTFsOut(), resultsAnnot(), stat = input$enrStatTF, n = input$numPlotEnrTF, plot.type = 'lollipop', pal.col = input$enrhPalColTF, col.hyper = NULL, col.hypo = NULL, thrs = as.numeric(unlist(strsplit(input$enrThrPlotTF, ","))), thrs.cols = input$enrThrPlotColTF)
                dev.off()
            }
        )


    }

    shiny::shinyApp(server=server, ui=ui)
}
