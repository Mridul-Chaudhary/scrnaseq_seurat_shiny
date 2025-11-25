library(shiny)
library(Seurat)
library(Azimuth)
library(ggplot2)
library(dplyr)

# SERVER for scRNA-seq Analysis App
shinyServer(function(input, output, session) {

  # ====== Global Seurat Object (will be updated at each step) ======
  seurat_obj <- reactiveVal(NULL)

  # ====== Tab 1: Data Upload ======
  observeEvent(input$upload_data, {
    req(input$counts_file) # require user input
    # Add checks for valid uploaded files
    counts_path <- input$counts_file$datapath
    barcodes_path <- input$barcodes_file$datapath
    features_path <- input$features_file$datapath

    # ReadMTX (or compatible loader)
    # COMMENT: Connects to 'Data Import' chunk in Rmd
    counts <- tryCatch({
      ReadMtx(mtx = counts_path, cells = barcodes_path, features = features_path)
    }, error = function(e) { NULL })

    # Create Seurat object
    if (!is.null(counts)) {
      seurat_obj(CreateSeuratObject(
        counts = counts,
        project = input$project_name,
        min.cells = input$min_cells,
        min.features = input$min_features
      ))
      output$upload_summary <- renderPrint({
        paste("Data uploaded for project:", input$project_name, "\nDimensions:", dim(counts))
      })
      output$data_upload_status <- renderUI({
        tags$span(style="color:green", "Data uploaded and Seurat object created.")
      })
    } else {
      output$upload_summary <- renderPrint({
        "Data not uploaded: Check input file formats and try again."
      })
      output$data_upload_status <- renderUI({
        tags$span(style="color:red", "Upload error or wrong format.")
      })
    }
  })

  # ====== Tab 2: Quality Control & Filtering ======
  observeEvent(input$qc_go, {
    ki <- seurat_obj()
    req(ki)
    # Add mitochondrial percentage calculation
    ki[["percent.mt"]] <- PercentageFeatureSet(ki, pattern = "^MT-")

    # Filter by QC parameters (user sliders)
    filtered <- subset(ki,
      subset =
        nFeature_RNA > input$feature_range[1] & nFeature_RNA < input$feature_range[2] &
        nCount_RNA > input$count_range[1] & nCount_RNA < input$count_range[2] &
        percent.mt > input$mt_range[1] & percent.mt < input$mt_range[2]
    )
    seurat_obj(filtered)

    # Output plots
    output$qc_violin <- renderPlot({
      VlnPlot(filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
    })
    output$qc_scatter <- renderPlot({
      FeatureScatter(filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    })
  })

  # ====== Tab 3: Normalization & Feature Selection ======
  observeEvent(input$norm_go, {
    ki <- seurat_obj()
    req(ki)
    ki <- NormalizeData(ki, normalization.method = input$norm_method, scale.factor = input$scale_factor)
    ki <- FindVariableFeatures(ki, selection.method = input$variable_method, nfeatures = input$nfeatures)
    top10 <- head(VariableFeatures(ki), 10)

    seurat_obj(ki) # update

    output$variable_feature_plot <- renderPlot({
      VariableFeaturePlot(ki)
    })
    output$top_variable_features <- renderPrint({
      top10
    })
  })

  # ====== Tab 4: Scaling & PCA ======
  observeEvent(input$pca_go, {
    ki <- seurat_obj()
    req(ki)
    all.genes <- rownames(ki)
    ki <- ScaleData(ki, features = all.genes)
    ki <- RunPCA(ki, features = VariableFeatures(ki), npcs = input$num_pcs)
    seurat_obj(ki)

    output$pca_loadings <- renderPlot({
      VizDimLoadings(ki, dims = 1:2, reduction = "pca")
    })
    output$pca_plot <- renderPlot({
      DimPlot(ki, reduction = "pca") + NoLegend()
    })
    output$pca_heatmap <- renderPlot({
      DimHeatmap(ki, dims = 1:input$num_pcs, cells = 500, balanced = TRUE)
    })
    output$elbow_plot <- renderPlot({
      ElbowPlot(ki)
    })
  })

  # ====== Tab 5: Clustering & UMAP ======
  observeEvent(input$cluster_go, {
    ki <- seurat_obj()
    req(ki)
    # Clustering and UMAP
    ki <- FindNeighbors(ki, dims = 1:input$neighbor_dims)
    ki <- FindClusters(ki, resolution = input$cluster_resolution)
    ki <- RunUMAP(ki, dims = 1:input$neighbor_dims)
    seurat_obj(ki)

    output$umap_plot <- renderPlot({
      DimPlot(ki, reduction = "umap", label = TRUE)
    })
    output$num_clusters <- renderPrint({
      paste("Number of clusters: ", length(unique(Idents(ki))))
    })
    output$cluster_table <- renderPlot({
      # plot distribution of clusters
      cluster.counts <- table(Idents(ki))
      barplot(cluster.counts, main = "Cluster Distribution", col = "steelblue")
    })
  })

  # ====== Tab 6: Marker Discovery ======
  observeEvent(input$markers_go, {
    ki <- seurat_obj()
    req(ki)
    # Marker selection
    markers <- FindAllMarkers(ki, only.pos = input$only_positive)
    # Filter
    filtered.markers <- markers %>%
      filter(avg_log2FC > input$logfc_threshold & pct.1 > input$pct_threshold)

    output$marker_table <- renderDataTable({
      filtered.markers
    })
    output$marker_heatmap <- renderPlot({
      # Heatmap of top marker genes (top 10 per cluster)
      top10 <- filtered.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
      DoHeatmap(ki, features = top10$gene)
    })
  })

  # ====== Tab 7: Azimuth Annotation ======
  observeEvent(input$azimuth_go, {
    ki <- seurat_obj()
    req(ki)
    # Azimuth annotation; assumes reference name is valid in env or installed
    # COMMENT: Connects to 'Azimuth' section in Rmd
    ki <- RunAzimuth(ki, reference = input$azimuth_ref)
    seurat_obj(ki)

    # Plots by predicted subclass/class
    output$azimuth_subclass_plot <- renderPlot({
      DimPlot(ki, group.by = "predicted.subclass", label = TRUE, label.size = 3, repel = T) + NoLegend()
    })
    output$azimuth_class_plot <- renderPlot({
      DimPlot(ki, group.by = "predicted.class", label = TRUE, label.size = 3, repel = T) + NoLegend()
    })
  })

  # ====== Tab 8: Visualization & Export ======
  observeEvent(input$render_plot, {
    ki <- seurat_obj()
    req(ki)
    output$custom_plot <- renderPlot({
      switch(input$plot_type,
        "UMAP" = DimPlot(ki, reduction = "umap", label = TRUE),
        "PCA" = DimPlot(ki, reduction = "pca"),
        "QC" = VlnPlot(ki, features = c("nFeature_RNA", "nCount_RNA", "percent.mt")),
        "Markers" = {
          top.markers <- tryCatch({
            FindAllMarkers(ki, only.pos = TRUE) %>% group_by(cluster) %>% top_n(10, avg_log2FC)
          }, error = function(e) { NULL })
          if (!is.null(top.markers)) DoHeatmap(ki, features=top.markers$gene) else NULL
        },
        "Azimuth" = DimPlot(ki, group.by = "predicted.class", label = TRUE, label.size=3, repel=T)
      )
    })
    # Download handler for plot export
    output$download_plot <- downloadHandler(
      filename = function() {
        paste("scRNA_plot_", input$plot_type, ".png", sep="")
      },
      content = function(file) {
        png(file)
        ki <- seurat_obj()
        plot_fun <- switch(input$plot_type,
          "UMAP" = function() DimPlot(ki, reduction="umap", label=TRUE),
          "PCA" = function() DimPlot(ki, reduction="pca"),
          "QC" = function() VlnPlot(ki, features = c("nFeature_RNA", "nCount_RNA", "percent.mt")),
          "Markers" = function() {
            top.markers <- tryCatch({
              FindAllMarkers(ki, only.pos = TRUE) %>% group_by(cluster) %>% top_n(10, avg_log2FC)
            }, error=function(e){NULL})
            if (!is.null(top.markers)) DoHeatmap(ki, features=top.markers$gene) else plot.new()
          },
          "Azimuth" = function() DimPlot(ki, group.by = "predicted.class", label = TRUE, label.size=3, repel=T)
        )
        print(plot_fun())
        dev.off()
      }
    )
  })

  # COMMENT: All sections above closely follow Rmd notebook headers and keep Seurat/Azimuth-only logic.
  # To extend: add more UI/server tabs, refine inputs/outputs, connect to more optional metadata as test data is added.
})
