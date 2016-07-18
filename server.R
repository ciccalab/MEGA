source("./functions/shiny_libs.R")
library(shiny)
library(DT)

options(shiny.maxRequestSize = 30*1024^2)
options(scipen=3)

make_pathway_list = function(pathMsigDbFile) {
  inputFile <- pathMsigDbFile
  con  <- file(inputFile, open = "r")

  c = 1
  pathway.list <- vector(mode="list",length=0)
  while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    myVector <- do.call("rbind",strsplit(oneLine, "\t"))
    t = vector(mode="list",length=1)
    t[[1]] = myVector[3:length(myVector)]
    names(t) = myVector[1]
    pathway.list = c(pathway.list,t)
    c = c+1
  }

  close(con)
  return(pathway.list)
}

shinyServer(function(input, output) {
      current = reactiveValues(res =NULL)


      RUN = eventReactive(input$submit, {
              A             = read.delim(input$cohort_A$datapath, header=T, stringsAsFactors = F)
              B             = read.delim(input$cohort_B$datapath, header=T, stringsAsFactors = F)
              gene.set.path = ifelse(!input$customGS,get.gene.set(as.numeric(input$gs_dataset)),input$gene_set$datapath)
              geneset       = make_pathway_list(gene.set.path)
              fdr_th        = input$fdr
              bootstrapping = input$bootstrapping=="True"
              nsim          = input$nsim

              if(!is.null(A) & !is.null(B) & !is.null(geneset) ){
                MEGA(A, B, geneset, fdr_th, bootstrapping, nsim)
              }

      })


      # TAB results
       output$mega_results   = DT::renderDataTable({
         results = RUN()
         current$res = results
         DT::datatable(results, options = list(orderClasses = TRUE, pageLength = 50))
       })
       

       output$downloadData <- downloadHandler(

         filename = function() {
           "MEGA-RVs_results.tsv"
         },
         content = function(con) {
           write.table(current$res,con,col.names = T, row.names = F, quote = F, sep="\t")
         }
       )
  })
