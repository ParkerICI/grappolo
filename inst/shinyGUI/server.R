

get_cluster_groups_table <- function(v, key) {
    tags$table(class = "table table-hover table-striped",
        tags$tr(tags$td(
            v[1],
            tags$button(class = "btn btn-xs btn-warning pull-right", onClick = sprintf("Shiny.onInputChange('clusteringui_remove_clustering_group', {'key':'%s', 'x':Math.random()})", key),
                tags$span(class = "glyphicon glyphicon-trash")
            )
        )),
        ifelse(length(v > 1),
            tagList(lapply(tail(v, n = -1), function(x) {tags$tr(tags$td(x))})),
            tagList()
        )
    )
}

render_clustering_ui <- function(working.directory, ...) {renderUI({
    fluidPage(
        fluidRow(
            column(12,
                selectInput("clusteringui_clustering_mode", "Clustering mode", choices = c("Single", "Pooled"), multiple = F),
                selectInput("clusteringui_file_for_markers", "Load marker names from file", choices = c("", list.files(path = working.directory, pattern = "*.fcs$")), width = "100%"),
                selectInput("clusteringui_markers", "Choose the markers for clustering", choices = c(""), multiple = T, width = "100%")
            )
        ),

        conditionalPanel(
            condition = "input.clusteringui_clustering_mode == 'Pooled'",
            fluidRow(
                column(6,
                    selectInput("clusteringui_files_list", label = "File list", choices = c(list.files(path = working.directory, pattern = "*.fcs$")),
                                selectize = F, multiple = T, width = "100%"),
                    actionButton("clusteringui_add_clustering_group", "Add clustering group")
                ),
                column(6,
                    uiOutput("clusteringui_clustering_groups_table")
                )
            ),
            fluidRow(
                column(12,
                    numericInput("clusteringui_downsample_to", "Pre-downsample data to (0 means no pre-downsampling)", value = 0)
                )
            )
        ),
        fluidRow(
            column(12,
                numericInput("clusteringui_num_clusters", "Number of clusters", value = 200, min = 1, max = 2000),
                numericInput("clusteringui_num_samples", "Number of samples (lower numbers lead to faster but less accurate results)", value = 50, min = 2),
                numericInput("clusteringui_asinh_cofactor", "Cofactor for asinh transformation", value = 5),
                numericInput("clusteringui_num_cores", "Number of CPU cores to use", value = 1),
                selectInput("clusteringui_negative_values", "Negative vaues", choices = c("truncate", "shift")),
                conditionalPanel(
                    condition = "input.clusteringui_negative_values == 'shift'",
                    numericInput("clusteringui_quantile_prob", "Quantile probability", value = 0.05, min = 0, max = 1, step = 0.01)
                ),
                actionButton("clusteringui_start", "Start clustering")
            )
        )
    )
})}

shinyServer(function(input, output, session) {
    working.directory <- dirname(file.choose())
    output$clusteringUI <- render_clustering_ui(working.directory, input, output, session)

    observe({
        if(!is.null(input$clusteringui_file_for_markers) && grepl("*.fcs$", input$clusteringui_file_for_markers)) {
            v <- grappolo:::get_fcs_col_names(file.path(working.directory, input$clusteringui_file_for_markers))
            updateSelectInput(session, "clusteringui_markers", choices = v)
        }
    })

    clusteringui_reactive_values <- reactiveValues(clustering_groups = NULL)

    output$clusteringui_clustering_groups_table = renderUI({
        dd <- clusteringui_reactive_values$clustering_groups
        return(tagList(mapply(get_cluster_groups_table, dd, names(dd), SIMPLIFY = F)))
    })

    observe({
        key <- input$clusteringui_remove_clustering_group$key
        if(!is.null(key) && key != "")
            isolate({
                clusteringui_reactive_values$clustering_groups[key] <- NULL
            })
    })

    observeEvent(input$clusteringui_add_clustering_group, {
        files_list <- isolate({input$clusteringui_files_list})
        clusteringui_reactive_values$clustering_groups <- c(clusteringui_reactive_values$clustering_groups,
                                                            setNames(list(files_list), files_list[1])
        )
    })

    observeEvent(input$clusteringui_start,
        isolate({
            showModal(modalDialog(
                sprintf("Clustering started with markers %s ",  paste(input$clusteringui_markers, collapse = " ")), br(),
                "please wait..."
            ))

            # This calls to force are necessary because the clustering functions will spawn
            # a different process, see https://github.com/rstudio/shiny/issues/2163

            num.cores <- force(input$clusteringui_num_cores)
            col.names <- force(input$clusteringui_markers)
            num.clusters <- force(input$clusteringui_num_clusters)
            asinh.cofactor <- force(input$clusteringui_asinh_cofactor)
            num.samples <- force(input$clusteringui_num_samples)
            downsample.to <- force(input$clusteringui_downsample_to)
            output.dir <- force(working.directory)
            negative.values <- force(input$clusteringui_negative_values)
            quantile.prob <- force(input$clusteringui_quantile_prob)

            if(input$clusteringui_clustering_mode == "Pooled") {
                input.files <- lapply(clusteringui_reactive_values$clustering_groups, function(s) {file.path(working.directory, s)})
                grappolo::cluster_fcs_files_groups(input.files,
                    num.cores = num.cores,
                    col.names = col.names,
                    num.clusters = num.clusters,
                    asinh.cofactor = asinh.cofactor,
                    num.samples = num.samples,
                    downsample.to = downsample.to,
                    negative.values = negative.values,
                    quantile.prob = quantile.prob,
                    output.dir = output.dir
                )
            } else {
                grappolo::cluster_fcs_files_in_dir(output.dir,
                    num.cores = num.cores,
                    col.names = col.names,
                    num.clusters = num.clusters,
                    asinh.cofactor = asinh.cofactor,
                    num.samples = num.samples,
                    output.dir = output.dir,
                    negative.values = negative.values,
                    quantile.prob = quantile.prob
                )
            }

            showModal(modalDialog(
                sprintf("Clustering completed. The output files can be found in %s", working.directory)

            ))
        })
    )

})
