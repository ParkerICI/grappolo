shinyUI(
    navbarPage("grappolo",
        tabPanel("Cluster data",
            uiOutput("clusteringUI")
        )
    )
)

