require(shiny)
require(shinyjs)
require(shinythemes)
require(DT)

panel_width = 3
result_with = 10

appCSS <- ".mandatory_star { color: red; }"

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(theme = shinytheme("united"),
    navbarPage("MEGA-RVs",

             tabPanel("MEGA-RVs",

                      titlePanel("MEGA-RVs"),
                      h4("Mutational Enrichment Gene set Analysis of Rare Variants"),
                      shinyjs::useShinyjs(),
                      shinyjs::inlineCSS(appCSS),

                      sidebarLayout(
                        sidebarPanel(
                                 fileInput('cohort_A', 'Cohort A (tab-separated file)*', accept = c('text/tab-separated-values', c('.tsv','.gz') )),
                                 fileInput('cohort_B', 'Cohort B (tab-separated file)*', accept = c('text/tab-separated-values', c('.tsv','.gz') )),
                                 selectInput("gs_dataset", "Select MsigDB GeneSet:",choices = list("KEGG"=1,"GO biological process"=2,"GO cellular component"=3,"GO molecular function"=4,"Reactome"=5,"BioCarta"=6,"Positional gene sets"=7,"Hallmark gene sets"=8,"Transcription factor targets"=9,"Cancer gene neighborhoods"=10,"Cancer modules"=11,"Oncogenic signatures"=12,"Disease associated"=13),selected = 1),
                                 selectInput("stat.test", 'Statistical Test', choices=c("Wilcoxon"=2,"GLM (Negative Binomial)"=1,"Kolmogorov–Smirnov"=3), selected=2),
                                 checkboxInput("customGS", label = tags$b("Use Custom Gene Set"), value = FALSE),
                                 conditionalPanel(condition = "input.customGS == true",fileInput('gene_set', label = NULL,accept = c('Gene Matrix Transposed','.gmt'))),
                                 selectInput("motecarlo", 'Monte Carlo Simulations', choices=c("True","False"), selected="False"),
                                 selectInput("bootstrapping", 'Bootstrapping', choices=c("True","False"), selected="False"),
                                 numericInput("nsim", "Number of random sampling (Monte Carlo and/or Bootstrapping)", value=1000),
                                 p("Mandatory fields are marked with *"),
                                 actionButton("submit", "Run MEGA-RVs", class = "btn-primary"),
                                 downloadButton('downloadData', 'Download Results')
                               ),

                        # MAIN PANEL
                        mainPanel(
                          dataTableOutput('mega_results')
                         )
                       )
             ),

             tabPanel("Help",
                      titlePanel("Help page"),
                      withTags({
                        div(class="header", checked=NA,
                            p("MEGA-RVs was developed to identify predefined gene
                              sets (e.g. genes involved in the same pathway, or predisposing
                              to specific diseases) that show a significantly higher number
                              of mutations in a group of samples as compared to another
                              group of samples.")
                        )}),
                      withTags({
                        div(class="body", checked=NA,
                            h4("Arguments"),
                            p(strong("A and B")," are data frame objects contaning the mutations counts. Coloums are samples, while rows are mutations. The first coloumn must always contain the name of the gene in which the mutation fall."),
                            p("Example:"),
                            p("Symbol\tS1\tS2\tS3"),
                            p("GeneA\t1\t0\t0"),
                            p("GeneB\t1\t0\t1"),
                            p("GeneC\t0\t0\t1"),
                            p("GeneD\t0\t1\t0"),
                            p("GeneE\t1\t1\t1"),
                            p(strong("gene.sets")," = List of gene sets. Each element of the list is set of genes and
                              the name of each element of the list must be the name of the gene set."),
                            p("Example:"),
                            p(code("gene.sets = list(g1=c('MAST2','ABCA10','ASPM'),g2=c('TP53','ZNF572','MYC'))")),

                            p(strong("fdr_th")," = false discovery rate threshold (default is 0.1)"),

                            p(strong("bootstrapping"), " = If equal to true the bootstrapping strategy is performed
                              to assess the effect of the sample size"),

                            p(strong("nsim")," = number of iterations used in the bootsrapping strategy (default is 1000)")
                        )})




                      )

    )
  )
)
