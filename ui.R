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
    navbarPage("MEGA-V",

             tabPanel("MEGA-V",

                      titlePanel("MEGA-V"),
                      h4("Mutational Enrichment Gene set Analysis of Variants"),
                      shinyjs::useShinyjs(),
                      shinyjs::inlineCSS(appCSS),

                      sidebarLayout(
                        sidebarPanel(
                                 fileInput('cohort_A', 'Cohort A (tab-separated file)*', accept = c('text/tab-separated-values', c('.tsv','.gz') )),
                                 fileInput('cohort_B', 'Cohort B (tab-separated file)*', accept = c('text/tab-separated-values', c('.tsv','.gz') )),
                                 selectInput("gs_dataset", "Select MsigDB GeneSet:",choices = list("KEGG"=1,"GO biological process"=2,"GO cellular component"=3,"GO molecular function"=4,"Reactome"=5,"BioCarta"=6,"Positional gene sets"=7,"Hallmark gene sets"=8,"Transcription factor targets"=9,"Cancer gene neighborhoods"=10,"Cancer modules"=11,"Oncogenic signatures"=12,"Disease associated"=13),selected = 1),
                                 # selectInput("stat.test", 'Statistical Test', choices=c("Wilcoxon"=2,"GLM (Negative Binomial)"=1,"Kolmogorov???Smirnov"=3), selected=2),
                                 selectInput("stat.test", 'Statistical Test', choices=c("Wilcoxon"=2,"Kolmogorov-Smirnov"=3,"Monte Carlo Permutations"=4), selected=2),
                                 checkboxInput("customGS", label = tags$b("Use Custom Gene Set"), value = FALSE),
                                 conditionalPanel(condition = "input.customGS == true",fileInput('gene_set', label = NULL,accept = c('Gene Matrix Transposed','.gmt'))),
                                 #selectInput("motecarlo", 'Monte Carlo Simulations', choices=c("True","False"), selected="False"),
                                 radioButtons("genome","Genome for Monte Carlo",choices=c("HG19","HG38"),selected="HG19",inline = TRUE),
                                 selectInput("bootstrapping", 'Bootstrapping', choices=c("True","False"), selected="False"),
                                 numericInput("nsim", "Number of random sampling (Monte Carlo or Bootstrapping)", value=1000),
                                 p("Mandatory fields are marked with *"),
                                 actionButton("submit", "Run MEGA-RVs", class = "btn-primary"),
                                 downloadButton('downloadData', 'Download Results')
                               ),

                        # MAIN PANEL
                        mainPanel(
                          dataTableOutput('mega_results')
                         )
                       )
             )

    )
  )
)
