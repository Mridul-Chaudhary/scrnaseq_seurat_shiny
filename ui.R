library(shiny)
library(shinythemes)

# UI for scRNA-seq Seurat + Azimuth analysis app
shinyUI(
  navbarPage(
    title = "Single Cell RNA-seq Analysis",
    theme = shinytheme("flatly"),

    # --- Tab 1: Data Upload ---
    tabPanel("Data Upload", 
      sidebarLayout(
        sidebarPanel(
          fileInput("counts_file", "Count Matrix (MTX/CSV/RDS)", accept = c(".mtx", ".csv", ".rds")),
          fileInput("barcodes_file", "Barcodes File (TXT/CSV)", accept = c(".txt", ".csv")),
          fileInput("features_file", "Features/Gene List (TXT/CSV)", accept = c(".txt", ".csv")),
          textInput("project_name", "Project Name", value = "sc_project"),
          numericInput("min_cells", "Min Cells per Gene", value = 3, min = 1),
          numericInput("min_features", "Min Features per Cell", value = 200, min = 1),
          actionButton("upload_data", "Go"),
          helpText("Upload single-cell matrix/barcodes/genes as your starting point; adjust initial filtering parameters.") # COMMENT: Corresponds to 'Data Import' in Rmd/server
        ),
        mainPanel(
          verbatimTextOutput("upload_summary"),
          uiOutput("data_upload_status")
        )
      )
    ),

    # --- Tab 2: Quality Control ---
    tabPanel("QC & Filtering",
      sidebarLayout(
        sidebarPanel(
          sliderInput("feature_range", "Genes per cell", min = 0, max = 10000, value = c(200, 5000)),
          sliderInput("count_range", "UMIs per cell", min = 0, max = 50000, value = c(500, 20000)),
          sliderInput("mt_range", "Mitochondrial (%)", min = 0, max = 100, value = c(0, 10)),
          actionButton("qc_go", "Run QC"),
          helpText("Choose QC/filter thresholds for cell and gene. Run QC to visualize distributions.") # COMMENT: Corresponds to 'Quality Control' in Rmd/server
        ),
        mainPanel(
          plotOutput("qc_violin"),
          plotOutput("qc_scatter")
        )
      )
    ),

    # --- Tab 3: Normalization & Feature Selection ---
    tabPanel("Normalization",
      sidebarLayout(
        sidebarPanel(
          selectInput("norm_method", "Normalization Method", choices = c("LogNormalize", "CLR", "RC", "SCT"), selected="LogNormalize"),
          numericInput("scale_factor", "Scale Factor", value = 10000),
          selectInput("variable_method", "Variable Feature Selection", choices = c("vst", "mean.var.plot", "dispersion"), selected="vst"),
          numericInput("nfeatures", "Number of Variable Features", value = 2000, min = 100),
          actionButton("norm_go", "Run Normalization"),
          helpText("Set normalization and feature selection methods.") # COMMENT: Corresponds to 'Normalization and Feature Selection' in Rmd/server
        ),
        mainPanel(
          plotOutput("variable_feature_plot"),
          verbatimTextOutput("top_variable_features")
        )
      )
    ),

    # --- Tab 4: Scaling & PCA ---
    tabPanel("Scaling & PCA",
      sidebarLayout(
        sidebarPanel(
          numericInput("num_pcs", "Number of PCs", value = 10, min = 2, max = 50),
          actionButton("pca_go", "Run Scaling & PCA"),
          helpText("Select dimensionality parameters for PCA/feature scaling.") # COMMENT: Corresponds to 'Scaling and PCA' in Rmd/server
        ),
        mainPanel(
          plotOutput("pca_loadings"),
          plotOutput("pca_plot"),
          plotOutput("pca_heatmap"),
          plotOutput("elbow_plot")
        )
      )
    ),

    # --- Tab 5: Clustering & UMAP ---
    tabPanel("Clustering & UMAP",
      sidebarLayout(
        sidebarPanel(
          numericInput("cluster_resolution", "Resolution for Clustering", value = 0.5, min = 0.1, max = 5, step = 0.1),
          numericInput("neighbor_dims", "PCs for Neighbors/UMAP", value = 10, min = 2, max = 50),
          actionButton("cluster_go", "Run Clustering & UMAP"),
          helpText("Adjust clustering/UMAP settings and run both.") # COMMENT: Corresponds to 'Clustering' and 'UMAP Dimensionality Reduction' in Rmd/server
        ),
        mainPanel(
          plotOutput("umap_plot"),
          plotOutput("cluster_table"),
          verbatimTextOutput("num_clusters")
        )
      )
    ),

    # --- Tab 6: Marker Discovery ---
    tabPanel("Marker Discovery",
      sidebarLayout(
        sidebarPanel(
          checkboxInput("only_positive", "Only Positive Markers", value = TRUE),
          numericInput("logfc_threshold", "min log2FC", value = 0.25),
          numericInput("pct_threshold", "min pct.1", value = 0.25),
          actionButton("markers_go", "Find Markers"),
          helpText("Find cluster marker genes by changing thresholds.") # COMMENT: Corresponds to 'Marker Discovery' in Rmd/server
        ),
        mainPanel(
          dataTableOutput("marker_table"),
          plotOutput("marker_heatmap")
        )
      )
    ),

    # --- Tab 7: Azimuth Annotation ---
    tabPanel("Azimuth Annotation",
      sidebarLayout(
        sidebarPanel(
          textInput("azimuth_ref", "Reference Name/Path", value = "humancortexref"),
          actionButton("azimuth_go", "Run Azimuth"),
          helpText("Map to Azimuth reference for cell type annotation.") # COMMENT: Corresponds to 'Azimuth Annotation' in Rmd/server
        ),
        mainPanel(
          plotOutput("azimuth_subclass_plot"),
          plotOutput("azimuth_class_plot")
        )
      )
    ),

    # --- Tab 8: Visualization & Export ---
    tabPanel("Visualization",
      sidebarLayout(
        sidebarPanel(
          selectInput("plot_type", "Choose Plot Type", choices = c("UMAP", "PCA", "QC", "Markers", "Azimuth")),
          actionButton("render_plot", "Show Plot"),
          downloadButton("download_plot", "Download Plot"),
          helpText("Select plot to view and export.") # COMMENT: Connects to multiple outputs in server, general visualization/export
        ),
        mainPanel(
          plotOutput("custom_plot")
        )
      )
    )
  )
)
