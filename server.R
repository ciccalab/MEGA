source("MEGA.R")
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
    print(c)
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


      RUN = eventReactive(input$submit, {

              A             = read.delim(input$cohort_A$datapath, header=T)
              rn = A[,1]
              cn = colnames(A)
              A = as.matrix(A[,-1]); rownames(A) = rn; colnames(A) = cn[2:length(cn)]
              print(A)

              B             = read.delim(input$cohort_B$datapath, header=T)
              rn = B[,1]
              cn = colnames(B)
              B = as.matrix(B[,-1]); rownames(B) = rn; colnames(B) =  cn[2:length(cn)]
              print(B)
              geneset       = make_pathway_list(input$gene_set$datapath)
              print(geneset)
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

         DT::datatable(results, options = list(orderClasses = TRUE, pageLength = 100))

       })

  })
