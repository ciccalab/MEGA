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
