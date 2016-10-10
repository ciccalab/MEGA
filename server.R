source("./functions/shiny_libs.R")

shinyInput <- function(FUN, len, id, ...) {
  inputs <- character(len)
  for (i in seq_len(len)) {
    inputs[i] <- as.character(FUN(paste0(id, i), ...))
  }
  inputs
}

options(shiny.maxRequestSize = 30*1024^2)
options(scipen=3)

shinyServer(function(input, output) {
      current = reactiveValues(res =NULL)


      RUN = eventReactive(input$submit, {
              A             = read.delim(input$cohort_A$datapath, header=T, stringsAsFactors = F)
              B             = read.delim(input$cohort_B$datapath, header=T, stringsAsFactors = F)
              gene.set.path = ifelse(!input$customGS,get.gene.set(as.numeric(input$gs_dataset)),input$gene_set$datapath)
              geneset       = read.gmt.file(gene.set.path)
              s.test        = input$stat.test
              fdr_th        = 0.1
              bootstrapping = input$bootstrapping=="True"
              montecarlo    = input$motecarlo=="True"
              nsim          = input$nsim
              genome        = input$genome
              
              if(!is.null(A) & !is.null(B) & !is.null(geneset) ){
                
                if (montecarlo)
                  load(paste("./RData/",genome,".gene.cds.length.RData",sep = ""))
                
                cpus = detectCores()
                if (!is.null(cpus))
                {
                  cpus = cpus - 1
                } else {
                  cpus = 2
                }
                
                MEGA(A, B, geneset, fdr_th, bootstrapping, nsim, s.test, montecarlo, gene.cds.length, cpus)
              }

      })


      # TAB results
      output$mega_results   = DT::renderDataTable({
        results = RUN()
        current$res = results
        DT::datatable(results, options = list(orderClasses = TRUE, pageLength = 15),selection= 'single',rownames = FALSE)
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
