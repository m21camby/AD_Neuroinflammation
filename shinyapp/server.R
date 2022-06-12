.libPaths( c( "/shinyserver/libs/AD_Neuroinflammation/prod/site-library" , .libPaths() ) )
#.libPaths( c( "/shinyserver/libs/AD_Neuroinflammation/prod/site-library" , .libPaths() ) )
# change to .libPaths( c( "/shinyserver/libs/AD_Neuroinflammation/prod/site-library" , .libPaths() ) )
# later when Dan change

lib.path=Sys.readlink("/shinyserver/libs/AD_Neuroinflammation/prod/site-library")

library(shiny, lib.loc = lib.path)

server <- function(input, output, session){

  # 1. QC plots 
  updateSelectizeInput(session, 'QC_Features', choices=QC_Features, selected=QC_Features[[1]])

  output$QC1 <- renderPlot({

        output$QC1 <- renderPlot({
          QC1(input$QC_Features)
        },height = 400,width = 500)
  })

  output$Violin_plot <- renderPlot({
    Violin_plot(input$QC_Features)
  },height = 400,width = 500)


  # 2. cell type plots
  updateSelectizeInput(session, 'celltype_Features', choices=celltype_Features, selected=celltype_Features[[1]])

  output$UMAP_celltype <- renderPlot(height = 500,width = 600, {
    UMAP_celltype(input$celltype_Features)
  })


  # 3. gene plots
  updateSelectizeInput(session, 'genes_Features', choices=genes_Features, selected=genes_Features[[1]])

  output$UMAP_genes <- renderPlot(height = 800, width = 900, {
    UMAP_genes(input$genes_Features)
  })

  output$DE_genes <- renderPlot(height = 700, width = 900, {
    DE_genes(input$genes_Features)
  })

  # 4. Trajectory gene plots
  updateSelectizeInput(session, 'Oligo_Genes', choices=Oligo_Genes, selected=Oligo_Genes[[1]])

  output$pseudograph <- renderPlot(height = 400, width = 700, {
    pseudograph(input$Oligo_Genes)
  })

}

