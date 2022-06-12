shinyUI(navbarPage(title="AD Neuroinflammation", id="navpage",
                   selected="About",
    tabPanel("project overview",
             fluidRow(column(12, includeMarkdown("info.md"))),
                            img(src = "bimsb_logo", width = 250, height = 160), img(src = "charite.png", width = 250, height = 140)),

    tabPanel("QC", id="clustering",
                            fluidPage(titlePanel("Select features"),
                                      sidebarLayout(sidebarPanel(
    selectInput(inputId = "QC_Features", label = "please select features", choices = QC_Features), width=2),
    mainPanel(verticalLayout(plotOutput("QC1"), plotOutput("Violin_plot"), br(), br(), includeMarkdown("exp1.md")))))),

    tabPanel("Clustering",
            fluidPage(titlePanel("Select cell type"),
                     sidebarLayout(sidebarPanel(selectInput(inputId = "celltype_Features", label = "please select cell type", choices = celltype_Features),width=2),
            mainPanel(verticalLayout(plotOutput("UMAP_celltype"), br(), br(), br(), includeMarkdown("exp2.md")))))),

    tabPanel("Gene",
             fluidPage(titlePanel("Select gene"),
                       sidebarLayout(sidebarPanel(selectInput(inputId = "genes_Features", label = "please select gene", choices = genes_Features),width=2),
                      mainPanel(verticalLayout(plotOutput("UMAP_genes", inline = TRUE), plotOutput("DE_genes", inline = TRUE)),
                                column(width = 12, align="center"))))),
    tabPanel("Oligo trajectory",
             fluidPage(titlePanel("Select gene"),
                       sidebarLayout(sidebarPanel(selectInput(inputId = "Oligo_Genes", label = "please select gene", choices = Oligo_Genes),width=2),
                                     mainPanel(verticalLayout(plotOutput("pseudograph"), br(), br(), includeMarkdown("exp3.md"))))))

))

