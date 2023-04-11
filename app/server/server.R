#### server.R

# ______________________________________________________________________________
# DESCRIPTION OF BUTTONS
# 1. input$jump_to_settings_preprocessing_B
# 2. input$preprocessing_B
# 3. input$jump2overview_B
# 4. input$back2presettings_B
# 5. input$back2presettings2_B

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
names_colorblind_palette <- c("light blue", "light red", "mustard", "green", "violett", "magenta", 
                              "turquoise", "dirty yellow", "dark pink", "merlot", "sky", "light grey")

# ______________________________________________________________________________
# SERVER ---------------------------------------------------------
server <- function(input, output, session) {
  ## 0. Server settings ####
  nextflow_script_dir = paste0(normalizePath("./server/scripts_nextflow/"), "/") # will be shifted to another place soon
  
  observeEvent(run_dir_selected(), {
    req(docker())
    print(docker()$run_dir)
    if (dir.exists(docker()$run_dir)){
      if (dir.exists(paste0(docker()$run_dir, "bam_merged_files"))){
        showModal(modalDialog(size = "l",
                              title = "Output directory is already present",
                              "The preprocessing of sequencing has started. Please go to next tab to select the analysis type of choice!",
                              easyClose = TRUE,
                              footer = NULL
        ))
      }
    }
  })
  
  ## 0.1 Volumes
  # Extract volumes for file choose function
  volumes = getVolumes()()
  
  
  ### 0.2 Cores for DEA ####
  # Set number of cores of default cores used for R process (from input in DEA tab)
  observe({register(MulticoreParam(input$cores))})
  
  
  ## RAM in use (currently not embedded into ui)
  output$memuse.out <- renderText({
    invalidateLater(5000, session)
    as.character(Sys.meminfo()[[2]])
  })
  
 
  metadata_raw <- data.frame(rbind(data.frame(Samples=rep("doubleclick to insert Samplename",1),	Condition=rep("doubleclick to insert Condition",1),	Replicate=rep("doubleclick to insert Replicate",1)),
    data.frame(Samples=rep("...",119),	Condition=rep("...",119),	Replicate=rep("...",119))))


  metadata_reac <- metadata_raw 

  

  ##0.3 Editable metadata file
  output$editable.metadata.table.out <- DT::renderDataTable({
    metadata_temp <- DT::datatable(
      metadata_raw,
      caption = 'Please make sure that samplenames fit names of sequencing output. If barcoded samples use barcode01-barcode96',
      editable = "cell",
      options = list(scrollY = 200,
                     scrollX = 200,
                     scroller = T,
                     fixedColumns = TRUE
                    ),
      selection = 'none',
      rownames = F)
    
  })

  observeEvent(input$editable.metadata.table.out_cell_edit, {
      info <- input$editable.metadata.table.out_cell_edit
      print(info)
      i <- as.numeric(info$row)
      j <- as.numeric(info$col+1)
      v <- as.character(info$value)
      metadata_reac[i, j] <<- v
    })

  output$saveMetadata = downloadHandler(
        filename = function(){paste0("metadata_", Sys.Date(), ".tsv")},
        content = function(fname) {
          download_metadata =  metadata_reac[which(metadata_reac$Condition != "..."),]
          print(download_metadata)
          write.table(download_metadata,fname,sep = "\t",row.names=F,col.names=T,quote = FALSE)
        }
    )


  ## 1. Preprocessing #######
  # Process variable for background processes executed with processx
  myProcess <- NULL
  
  # Not really needed
  alignment = reactiveVal()
  alignment(0)
  
  ### Buttons Actions ####
  # Dashboard --> Preprocessing 
  observeEvent(input$jump_to_settings_preprocessing_B, {
    newtab <- switch(input$tabs,
                     "start_preprocessing_page" = "pre_settings",
                     "pre_settings" = "start_preprocessing_page"
    )
    updateTabItems(session, "tabs", newtab)
  })
  
  # Preprocessing --> Run Overview 
  observeEvent(input$preprocessing_B,{
    print("change tabs to run overview")
    showModal(modalDialog(size = "l",
                          title = "Preprocessing started",
                          "The preprocessing of sequencing has started. You will be redirected once the preprocessing has finished!",
                          easyClose = TRUE,
                          footer = NULL
    ))
    shinyjs::disable("preprocessing_B")  # disable preprocessing Button, when clicked
    updateTabItems(session, "tabs", "run_overview")
    updateTabsetPanel(session, "main_quality_control",
                      selected = "read_length_dist_panel")
  }, 
  priority = 10)
  
  observeEvent(input$jump2overview_B, {
    updateTabsetPanel(session, "tabset1",
                      selected = "panel2")
  })
  
  observeEvent(input$back2presettings_B, {
    updateTabsetPanel(session, "tabset1",
                      selected = "panel1")
  })
  observeEvent(input$back2presettings2_B, {
    updateTabsetPanel(session, "tabset1",
                      selected = "panel2")
  })
  observeEvent(input$run_dea_B, {
    updateTabsetPanel(session, "dea.box",
                      selected = "dea.tab")
  })
  observeEvent(input$run_dt_preprocess, {
    updateTabsetPanel(session, "dea.box",
                      selected = "deu.tab")
  })
  
  observeEvent(input$missing_values_B, {
    showModal(modalDialog(size = "l",
                          title = "",
                          div(tags$b("Missing input data!", style = "color: red; font-size:100%;")),
                          easyClose = TRUE,
                          footer = NULL
    ))
  })
  
  observeEvent(input$jump_to_settings_overview_B, {
    updateTabsetPanel(session, "tabset1",
                      selected = "panel3")
  })
  
  
  # Confirm whole config 
  output$confirm_preprocess_settings <- renderUI({
    if (is.null(table_of_settings_transformed())){
      return()
    } else {
      actionButton(inputId = "preprocessing_B" , label = "Start", class = "btn btn-primary", 
                   style='padding:20px; 
                   font-size:200%; 
                   color: white; 
                   background-color: #d64045; 
                   border-radius: 5px;')
      
    }
  })
  
  # Start preprocessing is disabled when clicked
  # listen to changed settings
  listenChangedInput <- reactive({
    list(input$cores,
         input$config_file,
         input$preprocess,
         input$barcoded,
         input$DRS,
         input$fastq_files, 
         input$metadata_file, 
         input$genome_fasta_file, 
         input$transcriptome_fasta_file, 
         input$gtf_file, input$run_dir, input$bed_file,
         input$design_column, input$feature_A, input$feature_B, 
         input$input_data)
  })
  # when input changed, enable Start Preprocessing again
  observeEvent(listenChangedInput(), {
    #input$preprocess <- 0
    enable("preprocessing_B")
    updateActionButton(session, "preprocessing_B", "Start")
  })

  transform_path <- function(path){
    if (dir.exists("/NanopoReaTA_linux_docker")){
      if (!startsWith(path,"/NanopoReaTA_linux_docker")){
      transformed_path = paste0("/NanopoReaTA_linux_docker",path)
      }
      else {
        transformed_path=path
      }
    }
    else if (dir.exists("/NanopoReaTA_windows_docker")){
      if (!startsWith(path,"/NanopoReaTA_windows_docker")){ 
        transformed_path = gsub("\\\\", "/", as.character(path)) 
        transformed_path = sub(".*?\\/","",transformed_path)
        transformed_path = paste0("/NanopoReaTA_windows_docker/",transformed_path)
      }
      else {
        transformed_path = path
      }
    }
    else {
    transformed_path = path
    }
    return(transformed_path)
  }

  metadata_file_selected <- reactiveVal(NULL)
  fastq_files_selected <- reactiveVal(NULL)
  genome_fasta_file_selected <- reactiveVal(NULL)
  transcriptome_fasta_file_selected <- reactiveVal(NULL)
  bed_file_selected <- reactiveVal(NULL)
  gtf_file_selected <- reactiveVal(NULL)
  run_dir_selected <- reactiveVal(NULL)
  

  docker <- reactive({
    req(run_dir_selected())
    req(metadata_file_selected())
    req(fastq_files_selected())
    req(genome_fasta_file_selected())
    req(transcriptome_fasta_file_selected())
    req(gtf_file_selected())
    req(bed_file_selected())
    req(run_dir_selected())
    #if (run_dir_selected() != "" | metadata_file_selected() != "" |  fastq_files_selected() != "" | genome_fasta_file_selected() != "" | transcriptome_fasta_file_selected() != "" | gtf_file_selected() != "" | run_dir_selected() != "" | bed_file_selected() != ""){
    list(metadata_file = transform_path(metadata_file_selected()),
    fastq_files = transform_path(fastq_files_selected()),
    genome_fasta_file = transform_path(genome_fasta_file_selected()),
    transcriptome_fasta_file = transform_path(transcriptome_fasta_file_selected()),
    gtf_file = transform_path(gtf_file_selected()),
    bed_file = transform_path(bed_file_selected()),
    run_dir = transform_path(run_dir_selected())
    )
    #}
  })
  
  
  
  #Create button to confirm preprocessing config
  output$jump2overview_B.out <- renderUI({
    if (is.null(metadata_file_selected())){
      return()
    } else {
      #print(metadata_file_selected())
      if (run_dir_selected() == "" | metadata_file_selected() == "" |  fastq_files_selected() == ""){
        actionButton(inputId = "missing_values_B" ,
                     label = "",
                     icon = icon("arrow-alt-circle-right"),
                     class = "btn btn-primary",
                     style='font-size:200%; color: white; background-color: #gray; border-radius: 10px')
      } else {
        actionButton(inputId = "jump2overview_B",
                     label = "",
                     icon = icon("arrow-alt-circle-right"),
                     class = "btn btn-primary",
                     style='font-size:200%; color: white; background-color: #83c5be; border-radius: 10px')
      }
    }
  })
  
  ## INPUT FILES
  
  shinyFileChoose(input, "config_file", roots = volumes, session = session)
  shinyFileChoose(input, "metadata_file", roots = volumes, session = session)
  shinyDirChoose(input, "fastq_files", roots = volumes, session = session)
  shinyFileChoose(input, "genome_fasta_file", roots = volumes, session = session)
  shinyFileChoose(input, "transcriptome_fasta_file", roots = volumes, session = session)
  shinyFileChoose(input, "gtf_file", roots = volumes, session = session)
  shinyFileChoose(input, "bed_file", roots = volumes, session = session)
  shinyDirChoose(input, "run_dir", roots = volumes, session = session)
  
  config_selected <- eventReactive(input$config_file, {
    req(input$config_file)
    print(input$config_file)
    y = parseFilePaths(volumes, input$config_file)
    gsub("\\/\\/", "\\/", as.character(y$datapath))
  })
  
  observeEvent(input$metadata_file, {
    req(input$metadata_file)
    print(input$metadata_file)
    y = parseFilePaths(volumes, input$metadata_file)
    metadata_file_selected(gsub("\\/\\/", "\\/", as.character(y$datapath)))
  })
  
  observeEvent(input$fastq_files, {
    req(input$fastq_files)
    print(input$fastq_files)
    y = parseDirPath(getVolumes(), input$fastq_files)
    print(as.character(y))
    fastq_files_selected(paste0(gsub("\\/\\/", "\\/", as.character(y)), "/"))
  })
  
  observeEvent(input$genome_fasta_file, {
    req(input$genome_fasta_file)
    print(input$genome_fasta_file)
    y = parseFilePaths(volumes, input$genome_fasta_file)
    print(as.character(y))
    
    genome_fasta_file_selected(gsub("\\/\\/", "\\/", as.character(y$datapath)))
  })
  
  observeEvent(input$transcriptome_fasta_file, {
    req(input$transcriptome_fasta_file)
    print(input$transcriptome_fasta_file)
    y = parseFilePaths(volumes, input$transcriptome_fasta_file)
    transcriptome_fasta_file_selected(gsub("\\/\\/", "\\/", as.character(y$datapath)))
  })
  
  observeEvent(input$gtf_file, {
    req(input$gtf_file)
    print(input$gtf_file)
    y = parseFilePaths(volumes, input$gtf_file)
    gtf_file_selected(gsub("\\/\\/", "\\/", as.character(y$datapath)))
  })
  
  observeEvent(input$bed_file, {
    req(input$bed_file)
    print(input$bed_file)
    y = parseFilePaths(volumes, input$bed_file)
    print(y)
    bed_file_selected(gsub("\\/\\/", "\\/", as.character(y$datapath)))
  })
  
  observeEvent(input$run_dir, {
    req(input$run_dir)
    print(input$run_dir)
    y = parseDirPath(getVolumes(), input$run_dir)
    print(as.character(y))
    run_dir_selected(paste0(gsub("\\/\\/", "\\/", as.character(y)), "/"))
  })
  
  output$conig_file_out <- renderText({
    if(is.null(config_selected())){
      x = "No file loaded"
    } else {
      x = config_selected()
    }
    x
  })
  
  output$metadata_file_out <- renderText({
    metadata_file_selected()
  })
  
  output$fastq_file_out <- renderText({
    fastq_files_selected()
  })
  
  output$genome_fasta_file_out <- renderText({
    genome_fasta_file_selected()
  })
  
  output$transcriptome_fasta_file_out <- renderText({
    transcriptome_fasta_file_selected()
  })
  
  output$gtf_file_out <- renderText({
    gtf_file_selected()
  })
  
  output$bed_file_out <- renderText({
    req(bed_file_selected())
    bed_file_selected()
  })
  
  output$run_dir_out <- renderText({
    req(run_dir_selected())
    run_dir_selected()
  })
  
  
  # Load yaml file into R
  # Loaded this way to get values, that can be changed and saved
  table_of_settings_old <- reactiveValues()
  config_loaded <- eventReactive(input$config_file, {
    req(config_selected())
    print("Looks fine")
    print(config_selected())
    if (!is.null(config_selected())){
      table_of_settings_old <<- read_yaml(file = config_selected())
      return(1)
    } else {
      showModal(modalDialog(size = "l",
                            title = "Configuration file not loaded!",
                            "Please load a configuration file or fill out open fields manually!",
                            easyClose = TRUE,
                            footer = NULL
      ))
      return(0)
    }
  })
  

  
  observeEvent(config_loaded(), {
    req(config_loaded())
    if (config_loaded() == 1){
      # display message, when something is wrong with the parameters or fix simple things directly
      print(table_of_settings_old)
      ### METADATA
      if (!is.null(table_of_settings_old$metadata) | length(table_of_settings_old$metadata) != 0){
        if (!file.exists(transform_path(table_of_settings_old$metadata))){
          showModal(modalDialog(size = "l",
                                title = "File not found!", 
                                "Please make sure the metadata file exists and is of type 'txt' or 'tsv'", 
                                easyClose = T, 
                                footer = NULL))
          metadata_file_selected("") 
        } else if (!endsWith(transform_path(table_of_settings_old$metadata), ".txt") & !endsWith(transform_path(table_of_settings_old$metadata), ".tsv")){
          showModal(modalDialog(size = "l",
                                title = "Wrong file type!", 
                                "Please make sure the metadata file is of type 'txt' or 'tsv'", 
                                easyClose = T, 
                                footer = NULL))
          metadata_file_selected("") 
          
        } else {
          metadata_file_selected(table_of_settings_old$metadata)
        }
      }
      
      if (!is.null(table_of_settings_old$barcoded) | length(table_of_settings_old$barcoded) != 0){
        updateRadioButtons(session, 'barcoded', selected = table_of_settings_old$barcoded)
      } else {
        updateRadioButtons(session, "barcoded", selected = 1)
      }
      
      if (!is.null(table_of_settings_old$DRS) | length(table_of_settings_old$DRS) != 0){
        updateRadioButtons(session, 'DRS', selected = table_of_settings_old$DRS)
      } else {
        updateRadioButtons(session, "DRS", selected = 1)
      }
      
      
      if (!is.null(table_of_settings_old$general_folder)){
        print("yes")
        if (!file.exists(transform_path(table_of_settings_old$general_folder))){
          showModal(modalDialog(size = "l",
                                title = "Folder not found!", 
                                "Please make sure the path to the main directory exists and contains the specified folder.", 
                                easyClose = T, 
                                footer = NULL))
          fastq_files_selected("Not found")
        } else if (length(table_of_settings_old$general_folder) != 0 ){
          if (!endsWith(table_of_settings_old$general_folder, "/")) {
            table_of_settings_old$general_folder <<- paste0(table_of_settings_old$general_folder, "/")
            print("yes")
          }
          fastq_files_selected(table_of_settings_old$general_folder)
        } else {
          fastq_files_selected("Not found")
        }  
      } else {
        fastq_files_selected("Not found")
      }
      
      if (!is.null(table_of_settings_old$genome_fasta) | length(table_of_settings_old$genome_fasta) != 0){
        if (!file.exists(transform_path(table_of_settings_old$genome_fasta))){
          showModal(modalDialog(size = "l",
                                title = "File not found!", 
                                "Please make sure the genome (fasta) file exists and is in the correct folder!", 
                                easyClose = T, 
                                footer = NULL))
          genome_fasta_file_selected("Not found") 
        } else {
          genome_fasta_file_selected(table_of_settings_old$genome_fasta)
        } 
      } else {
        genome_fasta_file_selected("Not found")
      }
      
      if (!is.null(table_of_settings_old$transcriptome_fasta) | length(table_of_settings_old$transcriptome_fasta) != 0){
        if (!file.exists(transform_path(table_of_settings_old$transcriptome_fasta))){
          showModal(modalDialog(size = "l",
                                title = "File not found!", 
                                "Please make sure the transcriptome (fasta) file exists and is in the correct folder!", 
                                easyClose = T, 
                                footer = NULL))
          transcriptome_fasta_file_selected("Not found")
        } else {
          transcriptome_fasta_file_selected(table_of_settings_old$transcriptome_fasta)
        }
      } else {
        transcriptome_fasta_file_selected("Not found")
      }
      
      if (!file.exists(transform_path(table_of_settings_old$genome_gtf))){
        showModal(modalDialog(size = "l",
                              title = "File not found!", 
                              "Please make sure the gtf file exists and is in the correct folder!", 
                              easyClose = T, 
                              footer = NULL))
        gtf_file_selected("Not found")
      } else if (!endsWith(table_of_settings_old$genome_gtf, "gtf")){
        showModal(modalDialog(size = "l",
                              title = "Wrong file type!", 
                              "Please make sure the annotation file is of type 'gtf'", 
                              easyClose = T, 
                              footer = NULL))
        gtf_file_selected("Not found")
      } else if (!is.null(table_of_settings_old$genome_gtf) | length(table_of_settings_old$genome_gtf) != 0){
        
        gtf_file_selected(table_of_settings_old$genome_gtf)
      } else {
        gtf_file_selected("Not found")
      }
      
      if (!is.null(table_of_settings_old$bed_file) | length(table_of_settings_old$bed_file) != 0){
        if (!file.exists(transform_path(table_of_settings_old$bed_file))){
          showModal(modalDialog(size = "l",
                                title = "File not found!", 
                                "Please make sure the bed file exists and is in the correct folder!", 
                                easyClose = T, 
                                footer = NULL))
          bed_file_selected("Not found")
        } else if (!endsWith(table_of_settings_old$bed_file, "bed")){
          showModal(modalDialog(size = "l",
                                title = "Wrong file type!", 
                                "Please make sure the annotation file is of type 'bed'", 
                                easyClose = T, 
                                footer = NULL))
          bed_file_selected("Not found")
        } else {
          bed_file_selected(table_of_settings_old$bed_file)
        }
      } else {
        bed_file_selected("Not found")
      }
      
      
      if (!is.null(table_of_settings_old$run_dir) | length(table_of_settings_old$run_dir) != 0){
        if (!endsWith(table_of_settings_old$run_dir, "/")) {
          table_of_settings_old$run_dir <<- paste0(table_of_settings_old$run_dir, "/")
        }
        if (!file.exists(transform_path(table_of_settings_old$run_dir))){
          showModal(modalDialog(size = "l",
                                title = "Folder not found!", 
                                "Please make sure the output directory exists and lays the specified path.", 
                                easyClose = T, 
                                footer = NULL))
          run_dir_selected("Not found")
        } else {
          run_dir_selected(table_of_settings_old$run_dir)
        }
      } else {
        run_dir_selected("Not found")
      }
    }
  })
  
  # Save current settings in reactive variable
  table_of_settings <- reactive({
    #print(input$jump2overview_B)
    x = list("threads" = as.integer(input$cores),
             "barcoded" = as.integer(input$barcoded),
             "DRS" = as.integer(input$DRS),
             "metadata_file" = metadata_file_selected(),
             "fastq_files" = fastq_files_selected(),
             "genome_fasta_file" = genome_fasta_file_selected(),
             "transcriptome_fasta_file" = transcriptome_fasta_file_selected(),
             "gtf_file" = gtf_file_selected(),
             "bed_file" = bed_file_selected(),
             "run.dir" = run_dir_selected())
    print(unlist(x))
    x
  })


  #list(metadata_file = transform_windows(metadata_file_selected()),
  #  fastq_files = transform_windows(fastq_files_selected()),
  #  genome_fasta_file = transform_windows(genome_fasta_file_selected()),
  #  transcriptome_fasta_file = transform_windows(transcriptome_fasta_file_selected()),
  #  gtf_file =  transform_windows(gtf_file_selected()),
  #  bed_file =  transform_windows(bed_file_selected()),
  #  run_dir = transform_windows(run_dir_selected())

  table_of_settings_transformed <- reactive({
    req(docker())
    #print(input$jump2overview_B)
    x = list("threads" = as.integer(input$cores),
             "barcoded" = as.integer(input$barcoded),
             "DRS" = as.integer(input$DRS),
             "metadata_file" = docker()$metadata_file,
             "fastq_files" = docker()$fastq_files,
             "genome_fasta_file" = docker()$genome_fasta_file,
             "transcriptome_fasta_file" = docker()$transcriptome_fasta_file,
             "gtf_file" = docker()$gtf_file,
             "bed_file" = docker()$bed_file,
             "run.dir" = docker()$run_dir)
    print(unlist(x))
    x
  })
  
  output$metadata.tables.out <- renderUI({
    metadata_table_box(input, output, session)
  })
  
  output$design_column_in <- renderUI({
    if (is.null(metadata())) return()
    selectInput(inputId = 'design_column', 
                label = 'Select condition column', 
                choices = colnames(metadata()),selected = colnames(metadata())[2])#,
                # options = list(`actions-box` = TRUE))
  })
  
  output$feature_column_A <- renderUI({
    if (is.null(metadata())) return()
    selectInput(inputId = 'feature_A', 
                label = 'Select condition A', 
                choices = unique(metadata()[,input$design_column]))# ,
                # options = list(`actions-box` = TRUE))
  })
  
  output$feature_column_B <- renderUI({
    if (is.null(metadata())) return()
    x = unique(metadata()[,input$design_column])
    
    selectInput(inputId = 'feature_B', 
                label = 'Select condition B', 
                choices = x[!(x %in% input$feature_A)])#,
                # options = list(`actions-box` = TRUE))
  })
  

  output$color_feature_A <- renderUI({
    if (is.null(metadata())) return()
    choices_df = data.frame(
      names = names_colorblind_palette, 
      hex = safe_colorblind_palette
    )
    selectInput("color_A", 
                label = "Select color A", 
                choices = c("#88CCEE", "#DDCC77", "#CC6677", "#117733", "#332288", "#AA4499", 
                            "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888"),
                selected = "#88CCEE")
  })
  
  output$color_feature_B <- renderUI({
    if (is.null(metadata())) return()
    choices_df = data.frame(
      names = names_colorblind_palette, 
      hex = safe_colorblind_palette
    )
    selectInput(inputId = "color_B", 
            label = "Select color B", 
            choices = c("#88CCEE", "#DDCC77", "#CC6677", "#117733", "#332288", "#AA4499", 
                            "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888"),
            selected = "#DDCC77")
  })

  condi_cols = reactive({
    req(metadata())
    condi_cols = c(input$color_A, input$color_B)
    names(condi_cols) = c(input$feature_A, input$feature_B)
    condi_cols
  })
  
  output$design.matrix.out <- renderUI({
    design_matrix_box(input, output, session, F)
  })
  metadata <- reactive({
      if (table_of_settings_transformed()$metadata != ""){
        data.table::fread(table_of_settings_transformed()$metadata, data.table = F)
      } else {
        return (NULL)
      }
      
  })
  
  output$metadata.out <- renderDT({
    if (is.null(metadata())){
      h5("No metadata")
    } else {
      metadata() #< reactive metadata value
    }
  })
  
  # Print table of settings
  output$table_of_settings_df <- DT::renderDataTable({
    req(table_of_settings_transformed())
    config = table_of_settings_transformed()
    config$conditions = "No multiple conditions"
    
    if (is.null(metadata())) return()
    if (!is.null(input$design_column)){
      if (!(input$design_column == "")){
        config$conditions = paste0(input$design_column, ": ", input$feature_A, " vs. ", input$feature_B)
        
      }
    }
    config_out = data.frame(Configurations = unlist(config), row.names = names(config))
    DT::datatable(
      config_out,
      extensions = c('Buttons', 'Scroller'),
      options = list(scrollY = 500,
                     scrollX = 500,
                     deferRender = TRUE,
                     scroller = F,
                     paging = TRUE,
                     dom = 't',
                     fixedColumns = TRUE),
      rownames = T)
  })
  
  table_of_genes <- reactiveVal()

  ### Start nextflow pipeline ####
  process_running <- reactiveValues(x = NULL)

  observeEvent(input$preprocessing_B, {
    req(table_of_settings_transformed())
    req(metadata())
    req(run_dir_selected())
    req(alignment())
    req(docker())
    if (input$preprocessing_B == 1 & !is.null(metadata())){
      if (!dir.exists(docker()$run_dir)) dir.create(docker()$run_dir)
      dir.create(paste0(docker()$run_dir, "nextflow_work_dir"))

      # Save config file
      config = table_of_settings()
      names(config) = c("threads", "barcoded", "DRS", "metadata", "general_folder", "genome_fasta", "transcriptome_fasta", "genome_gtf", "bed_file", "run_dir")
      currentwd = dirname(getwd())
      print(currentwd)
      write_yaml(x = config, file = paste0(docker()$run_dir, "config.yaml"))
      print(list.files(paste0(currentwd, "/server/bash_scripts"), ".sh"))

      config = table_of_settings_transformed()
      names(config) = c("threads", "barcoded", "DRS", "metadata", "general_folder", "genome_fasta", "transcriptome_fasta", "genome_gtf", "bed_file", "run_dir")
      config[["script_dir"]] = nextflow_script_dir
      write_yaml(x = config, file = paste0(docker()$run_dir, "cf_transformed.yaml"))
      print(list.files(paste0(currentwd, "/server/bash_scripts"), ".sh"))
      print(config)

      

      if(input$preprocess == 1){
        
        if(!is.null(myProcess)){
          if(myProcess$is_alive()){
             myProcess$kill()
          }
        }
        
        myProcess <<- processx::process$new(command = paste0(currentwd, "/app/server/bash_scripts/run_nextflow.sh"), 
                                            args = c(paste0(config[["script_dir"]], "new_np_pipe.nf"), paste0(docker()$run_dir, "cf_transformed.yaml"), paste0(docker()$run_dir, "nextflow_work_dir")), 
                                            echo_cmd = T,
                                            stderr = paste0(docker()$run_dir,"error.log"), 
                                            stdout = paste0(docker()$run_dir,"index.log"))

        write.table(c(1), file = paste0(docker()$run_dir, "process_running.txt"), col.names = F, row.names = F, quote = F)
      } else {
        myProcess <<- processx::process$new(command = "echo", args = "Preprocessing skipped", echo_cmd = T,stderr =paste0(docker()$run_dir,"error.log"), stdout = paste0(docker()$run_dir,"index.log"))
        write.table(c(0), file = paste0(docker()$run_dir, "process_running.txt"), col.names = F, row.names = F, quote = F)

      }
      if (!myProcess$is_alive()) print("DEAD")
    }
    alignment(1)
    print(alignment())
    
  }, priority = 1)
  
  ##### Load gtf file ####
  observeEvent({input$tabs}, {
    # Load gtf file
    req(docker()$gtf_file)
    req(docker())
    if (is.null(table_of_genes())){
      if(input$tabs == "run_overview"){
        print("load gtf started")
        if (file.exists(paste0(docker()$run_dir, "converted_gtf.csv"))){
          g = read.csv(paste0(docker()$run_dir, "converted_gtf.csv"), header = T)
        } else {
          g = getGeneSymbolFromGTF(docker()$gtf_file, ".")
        }
        table_of_genes(unique(g[,c("gene_id", "gene_name")]))
        
        print("loading finished!")
        write.csv(g[, c("gene_id", "gene_name", "transcript_id", "transcript_name")], paste0(docker()$run_dir, "converted_gtf.csv"))
      }
    }
  }, priority = 0)
  
  ## 1.1 Output from Preprocessing ####
  
  ### Load counts file ####
  countsfile = reactiveValues(df = NULL, df_trans = NULL)

  observeEvent({input$preprocessing_B + input$run_dea_B + input$submit_gene_selection}, {
    print("In Gene Count")
    print(paste0("Current tab: ", input$tabs))
    print(input$preprocessing_B)
    if (sum(c(input$preprocessing_B) > 0)){
      print("###################################")
      print("##     Load/Update Counts file   ##")
      print("###################################")
      print(paste0(Sys.time(), ": Starting count file loading/updating"))

      req(metadata())
      counts_file = paste0(table_of_settings_transformed()$run.dir, "merged_all.csv")

      print(sprintf("Count file does exists?: %s", file.exists(counts_file)))
      if (!file.exists(counts_file)) {
        showModal(modalDialog(size = "l",
                              title = "No gene count file",
                              "No gene count file found! Wait until the preprocessing pipeline has generated a gene count file! This can take a few minutes",
                              easyClose = TRUE,
                              footer = NULL
        ))
      } else {
        x = fread(counts_file, data.table = F)
        geneNames = x$Geneid
        x$Geneid <- NULL
        y = sapply(x, as.numeric)
        row.names(y) = geneNames
        nulls = sapply(0,function(x)rowSums(y==x))
        num_nulls = length(which(nulls[,1] == 0))
        if (num_nulls == 0){
          showModal(modalDialog(size = "l",
                                title = "Warning",
                                "The gene counts file contains too many 0. All samples from your metadata table must match the sample names used for barcoding!",
                                easyClose = TRUE,
                                footer = NULL
          ))
        }
        countsfile$df <- as.data.frame(y)
        print("Counts file loaded successfully!")
      }
    } 
  }, 
  priority = 2
  )
  
  
  observeEvent({input$preprocessing_B + input$run_dt_preprocess}, {
    print("In Transcript Count")
    req(input$preprocessing_B)
    req(input$run_dt_preprocess)
    req(docker())
    req(metadata())
    print(paste0("Current tab: ", input$tabs))
    print("###################################")
    print("##     Load/Update Transcript Counts file   ##")
    print("###################################")
    print(paste0(Sys.time(), ": Starting transcript count file loading/updating"))
    
    counts_file = paste0(table_of_settings_transformed()$run.dir, "salmon_merged_absolute.csv")
    
    print(sprintf("Salmon Count file does exists?: %s", file.exists(counts_file)))
    if (file.exists(counts_file)) {
      x = fread(counts_file, data.table = T)
      geneNames = x$Name
      x = x %>% dplyr::select(-1:-4)
      x = as.data.frame(x)
      y = sapply(x, as.numeric)
      row.names(y) = geneNames
      nulls = sapply(0,function(x)rowSums(y==x))
      num_nulls = length(which(nulls[,1] == 0))
      if (num_nulls == 0){
        showModal(modalDialog(size = "l",
                              title = "Warning",
                              "The gene counts file contains too many 0. All samples from your metadata table must match the sample names used for barcoding!",
                              easyClose = TRUE,
                              footer = NULL
        ))
      }
      countsfile$df_trans <- as.data.frame(y)
      print("Transcript counts file loaded successfully!")

    }
  }, 
  priority = 2
  )


  ## 3. Run Overview ####
  
  ### Read length distribution ####
  # store inputs from reactiveFileReader later
  overview_variables <- reactiveValues(
    readLengthTable = NULL
  )

  observeEvent(ignoreNULL = TRUE,
               eventExpr = {
                 input$preprocessing_B
               },
               handlerExpr = {
                 overview_variables$readLengthTable <- reactivePoll(10000, session,
                                                                       checkFunc = function(){
                                                                         req(docker()$run_dir)
                                                                        inputSizeDir = paste0(docker()$run_dir, "ReadLengthFolder")
                                                                        if (dir.exists(inputSizeDir)){
                                                                          filesList <- list.files(inputSizeDir, full.names = TRUE)
                                                                          file.info(filesList[1])$mtime[1]
                                                                        } else {
                                                                          print("Directory does not exists!")
                                                                        }
                                                                       },
                                                                       valueFunc = function(){
                                                                        req(docker()$run_dir)
                                                                        inputSizeDir = paste0(docker()$run_dir, "ReadLengthFolder")
                                                                        if (dir.exists(inputSizeDir)){
                                                                          filesList <- list.files(inputSizeDir, full.names = TRUE)
                                                                          if (length(filesList) == 0) return()

                                                                          names(filesList) = gsub("_read_lengths_pass.txt", "", basename(filesList))

                                                                            # Load files and add pass/fail column
                                                                          readLengths <- lapply(filesList, read.table, header = T)
                                                                          readLengths
                                                                        } else {
                                                                          print("Dir does not exists!")
                                                                          return()
                                                                        }
                                                                      }
                                                                    )
                 }
  )

  readLengthTab <- reactive({
    req(docker()$run_dir)
    
    req(overview_variables$readLengthTable())
    print("### Overview Length reactive")
    print(length(overview_variables$readLengthTable()))
    if (is.null(overview_variables$readLengthTable())) {
      return()
    } else {
      readLengths = overview_variables$readLengthTable()
      # Filter out lists of no reads
      readLengths = readLengths[unlist(lapply(readLengths, function(x) nrow(x) > 0))]
      # Create data.frame of list
      readLengths_df = data.table::rbindlist(readLengths, idcol = "Sample")
      upper_bound <- quantile(readLengths_df$Length, 0.99)
      readLengths_df_filt = readLengths_df[readLengths_df$Length < upper_bound,]
      readLengths_df_filt$Sample = gsub("_.*", "", readLengths_df_filt$Sample)
      readLengths_df_filt

    }
  })

  output$samplewise_read_length.out <- renderPlot({
    req(readLengthTab())
    print("Samplewise read length plot")
    print(length(readLengthTab()))
    if(is.null(readLengthTab())) return()
    g = createLengthPlots(readLengthTab(), metadata(), input$design_column, c(input$feature_A, input$feature_B), c(input$color_A, input$color_B))
    g[[1]]
  }, bg = "transparent")

  output$groupwise_read_length.out <- renderPlot({
    req(readLengthTab())
    if(is.null(readLengthTab())) return()
    g = createLengthPlots(readLengthTab(), metadata(), input$design_column, c(input$feature_A, input$feature_B), c(input$color_A, input$color_B))
    g[[2]]
  }, bg = "transparent")

  # downloads

  output$samplewise_read_length.down <- downloadHandler(
    filename = function(){
      paste0("Read_length_distribution_samplewise_", Sys.Date(), ".", input$samplewise_read_length.type)
    },
    content = function(file){
      if(is.null(readLengthTab())) return()
      plot = samplewise_read_length.download(readLengthTab(), metadata(), input$design_column, c(input$feature_A, input$feature_B))
      ggsave(file, plot = plot, device = input$inner_var_sample.type, height = 5, width = 14, units = "in")
    })
  output$groupwise_read_length.down <- downloadHandler(
    filename = function(){
      paste0("Read_length_distribution_groupwise_", Sys.Date(), ".", input$groupwise_read_length.type)
    },
    content = function(file){
      if(is.null(readLengthTab())) return()
      plot = groupwise_read_length.download(readLengthTab(), metadata(), input$design_column, c(input$feature_A, input$feature_B), c(input$color_A, input$color_B))
      ggsave(file, plot = plot, device = input$groupwise_read_length.type, height = 5, width = 14, units = "in")
    })

  output$process_time.down <- downloadHandler(
    filename = function(){
      paste0("Process_time_tools_over_time_", Sys.Date(), ".", input$process_time.type)
    },
    content = function(file){
      req(process_time$df())
      if (!is.null(process_time$df())){
        print(head(process_time$df()))
        x = process_time$df()
        # replace _ with whitespace in tool names for a cleaner legend
        x = x %>% mutate(Tool = gsub("_", " ", Tool))
        
        # get the x-tick interval based on the number of iterations in the data
        max_it <- max(x$Iteration)
        xticks <- 1:max_it
        if (max_it <= 10) { # every tick
          xticks <- xticks
        } else if (max_it <= 50) { # every 5th tick
          xticks <- xticks[xticks %% 5 == 0]
        } else if (max_it <= 100) { # every 10th tick
          xticks <- xticks[xticks %% 10 == 0]
        } else if (max_it <= 1000) { # every 50th tick
          xticks <- xticks[xticks %% 50 == 0]
        } else if (max_it <= 10000) { # every 100th tick
          xticks <- xticks[xticks %% 100 == 0]
        } else { # every 1000th tick
          xticks <- xticks[xticks %% 1000 == 0]
        }
        plot = ggplot(x, aes(x = Iteration, y = Time, fill = Tool)) +
          geom_bar(stat = "identity") +
          scale_x_continuous(breaks = xticks) + # label each x tick individually + # label each x tick individually
          ylab("Time [s]") + # add [s] to y-axis label
          theme_bw() +
          scale_fill_manual("Preprocessing steps", values = c("#CC6677", "#d8c66c", "#117733", "#88CCEE", "#AA4499", "#24011e","#092a36")) +
          theme(
          panel.grid.minor = element_blank(), # remove minor grid lines
          panel.grid.major.x = element_blank(),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          axis.text = element_text(angle = 45, hjust = 1, size = 17),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 23),
          axis.title = element_text(size = 23))
        ggsave(file, plot = plot, device = input$process_time.type, height = 10, width = 10, units = "in")
      }
  })
  
  ### Gene Expression Variability ####
  
  # store inputs from reactiveFileReader later
  overview_variables_infer <- reactiveValues(
    inner_var_sample.table = NULL, 
    expGenes_counted.table = NULL
  )
  
  observeEvent(input$preprocessing_B,{
    overview_variables_infer$inner_var_sample.table <- reactivePoll(10000, session,
                                   checkFunc = function(){
                                     req(docker()$run_dir)
                                     if (dir.exists(docker()$run_dir)) {
                                        inner_var=paste0(docker()$run_dir, "inner_variability_per_sample.csv")
                                        if (file.exists(inner_var)){   
                                          file.info(inner_var)$mtime[1]
                                        }
                                     } else {
                                       print("Inner variability file does not exists!")
                                     }
                                   },
                                   valueFunc = function(){
                                     req(docker()$run_dir)
                                     if (!dir.exists(docker()$run_dir)) return(NULL)
                                     inner_var=paste0(docker()$run_dir, "inner_variability_per_sample.csv")
                                     
                                     if (file.exists(inner_var)) {
                                      read.table(inner_var, sep="\t", header=T)
                                     } else {
                                      NULL
                                     }
                                   })
    overview_variables_infer$expGenes_counted.table <- reactivePoll(10000, session,
                                                              checkFunc = function(){
                                                                req(docker()$run_dir)
                                                                if (dir.exists(docker()$run_dir)) {
                                                                  genes_counted=paste0(docker()$run_dir, "exp_genes_counted_per_sample.csv")
                                                                    if (file.exists(genes_counted)){   
                                                                      file.info(genes_counted)$mtime[1]
                                                                    }
                                                                } else {
                                                                  print("Counts file does not exists!")
                                                                }
                                                              },
                                                              valueFunc = function(){
                                                                req(docker()$run_dir)
                                                                if (!dir.exists(docker()$run_dir)) return(NULL)
                                                                genes_counted=paste0(docker()$run_dir, "exp_genes_counted_per_sample.csv")
                                                                
                                                                if (file.exists(genes_counted)) {
                                                                  read.table(genes_counted, sep="\t", header=T)
                                                                } else {
                                                                  NULL
                                                                }
                                                              })
  },
  priority = -1
  )
  

  output$inner_var_sample.out <- renderPlot({
    if (is.null(overview_variables_infer$inner_var_sample.table())){ 
      inner_var_plot_per_sample(data.frame())
    }
    else{
    inner_var_plot_per_sample(overview_variables_infer$inner_var_sample.table())
    }
  }, bg = "transparent")
      
  output$inner_var_condition.out <- renderPlot({
    req(metadata())
    if (is.null(overview_variables_infer$inner_var_sample.table())){
        inner_var_plot_per_condition(data.frame(), metadata(), c(input$color_A, input$color_B))
    } 
    else{
    inner_var_plot_per_condition(overview_variables_infer$inner_var_sample.table(), metadata(), c(input$color_A, input$color_B))
    }
  }, bg = "transparent")
  
  output$total_genes_sample.out <- renderPlot({
    req(overview_variables_infer$expGenes_counted.table())

    if (is.null(overview_variables_infer$expGenes_counted.table())){
        total_genes_counted_plot_per_sample(data.frame())
      }
    else{
      total_genes_counted_plot_per_sample(overview_variables_infer$expGenes_counted.table())
    }
    }, bg = "transparent")
      
  output$total_genes_condition.out <- renderPlot({
    req(metadata())
    if(is.null(overview_variables_infer$expGenes_counted.table())){
      total_genes_counted_plot_per_condition(data.frame(), metadata(), c(input$color_A, input$color_B))
    }
    else{
      total_genes_counted_plot_per_condition(overview_variables_infer$expGenes_counted.table(), metadata(), c(input$color_A, input$color_B))
    }
  
  }, bg = "transparent")
    
  # downloads 
  

  output$inner_var_sample.down <- downloadHandler(
   filename = function(){
    paste0("Expression_variability_over_time_samplewise_", Sys.Date(), ".", input$inner_var_sample.type)
    }, 
   content = function(file){
    if(is.null(overview_variables_infer$inner_var_sample.table())) return()
    plot = inner_var_plot_per_sample.download(overview_variables_infer$inner_var_sample.table())
    ggsave(file, plot = plot, device = input$inner_var_sample.type, height = 5, width = 14, units = "in")
   })



  output$inner_var_condition.down <- downloadHandler(
    filename = function(){
      paste0("Expression_variability_over_time_conditionwise_", Sys.Date(), ".", input$inner_var_condition.type)
    }, 
    content = function(file){
      if(is.null(overview_variables_infer$inner_var_sample.table()) | is.null(metadata())) return()
      plot = inner_var_plot_per_condition.download(overview_variables_infer$inner_var_sample.table(), metadata(), c(input$color_A, input$color_B))
      ggsave(file, plot = plot, device = input$inner_var_condition.type, height = 5, width = 14, units = "in")
    })

  

  output$total_genes_sample.down <- downloadHandler(
    filename = function(){
      paste0("Gene_expression_counts_over_time_samplewise_", Sys.Date(), ".", input$total_genes_sample.type)
    }, 
    content = function(file){
      if(is.null(overview_variables_infer$expGenes_counted.table())) return()
      plot = total_genes_counted_plot_per_sample.download(overview_variables_infer$expGenes_counted.table())
      ggsave(file, plot = plot, device = input$total_genes_sample.type, height = 5, width = 14, units = "in")
    })


 
  output$total_genes_condition.down <- downloadHandler(
    filename = function(){
      paste0("Gene_expression_counts_over_time_conditionwise_", Sys.Date(), ".", input$total_genes_condition.type)
    }, 
    content = function(file){
      if(is.null(overview_variables_infer$expGenes_counted.table()) | is.null(metadata())) return()
      plot = total_genes_counted_plot_per_condition.download(overview_variables_infer$expGenes_counted.table(), metadata(), c(input$color_A, input$color_B))
      ggsave(file, plot = plot, device = input$total_genes_condition.type, height = 5, width = 14, units = "in")
    })

  
  
  
  
  ## 2. Quality Control ####
  ### Buttons ####
  
  output$submit_gene_selection.out <- renderUI({
    if (!is.null(table_of_genes())){
         actionButton("submit_gene_selection", "Submit genes",
                      class = "btn btn-primary",
                      style='font-size:200%; color: white;
                             background-color: #e76f51; 
                             border-radius: 5px')
      } else {
        actionButton("submit_gene_selection_disabled", "Submit genes",
                     class = "btn btn-primary",
                     style='font-size:200%; color: white;
                            background-color: gray; 
                            border-radius: 5px')
    }
  })
  
  ### Gene-wise analysis #####
  table_of_normCounts <- reactiveValues(df_norm = NULL)
  
  observeEvent({input$run_dea_B + input$submit_gene_selection}, {
    # req(metadata())
    print(paste0("Norm current tab: ", input$tabs))
    req(table_of_settings_transformed())
    req(countsfile$df)
    print(!is.null(countsfile$df))
    if (is.null(countsfile$df)){
      showModal(modalDialog(size = "l",
                            title = "No countsfile loaded!",
                            "The countsfile has not been loaded so the selected analysis is not possible, yet. 
                            Please proceed here when it has been loaded as indicated in 'Run Overview'. ",
                            easyClose = TRUE,
                            footer = NULL
      ))
    }
    print(!is.null(metadata()))
    if (!is.null(countsfile$df) & !is.null(metadata())){
      if (is.null(input$design_column)){
        print("input$design_column is null")
        return()
      }
      print(input$design_column)
      table_of_normCounts$df_norm <- createDDS(counts.file = countsfile$df,
                meta.file = metadata(),
                condition.col = input$design_column,
                first.level = input$feature_A,
                ref.level = input$feature_B)
    } else {
      print(paste0(Sys.time(), ": Analysis not started!"))
      table_of_normCounts$df_norm <- NULL
    }
  }, 
  priority = 3)
  
   
  output$counts_plots_box <- renderUI({
    input$submit_gene_selection
    # cat(unlist(unlist(input$submit_gene_selection)))
    if (is.null(input$submit_gene_selection)) return()
    if (input$submit_gene_selection <= 0){
      box(title = "Counts plots", status = "primary", 
          solidHeader = T, 
          collapsible = T, 
          collapsed = T, 
          width = 12,
          # cat(input$submit_gene_selection),
          fluidRow(
            column(12, radioButtons("disp", "Display",
                                    choices = c("Boxplot" = "boxplot",
                                                "Dotplot" = "dotplot",
                                                "Violin plot" = "violinplot"),
                                    selected = "boxplot", inline = T)),
            column(12, 
                   div(style="display:inline-block; vertical-align:top", downloadButton('down_tea')), 
                   div(style="display:inline-block; width:80px; vertical-align:top", selectInput("tea_down.type", NULL, choices = c("png", "pdf")))),
            column(12, plotOutput("teaPlot", height = "500px", click = clickOpts(id="plot_click"),
                                  hover = hoverOpts(id = "plot_hover", delayType = "throttle"),
                                  brush = brushOpts(id = "plot_brush"))  %>% withSpinner(color="#0dc5c1"))
          ))
    } else {
      # if (input$submit_gene_selection == 0) return()
      box(title = "Counts plots", status = "primary", 
          solidHeader = T, 
          collapsible = T, 
          collapsed = F, 
          width = 12,
          # cat(input$submit_gene_selection),
          fluidRow(
            column(12, radioButtons("disp", "Display",
                                    choices = c("Boxplot" = "boxplot",
                                                "Dotplot" = "dotplot",
                                                "Violin plot" = "violinplot"),
                                    selected = "boxplot", inline = T)),
            column(12, 
                   div(style="display:inline-block; vertical-align:top", downloadButton('down_tea')), 
                   div(style="display:inline-block; width:80px; vertical-align:top", selectInput("tea_down.type", NULL, choices = c("png", "pdf")))),
            column(12, plotOutput("teaPlot", height = "500px", click = clickOpts(id="plot_click"),
                                  hover = hoverOpts(id = "plot_hover", delayType = "throttle"),
                                  brush = brushOpts(id = "plot_brush"))  %>% withSpinner(color="#0dc5c1"))
          ))
    }
  })
  
  
  genes.list <- reactive({
    req(table_of_genes())
    gtf_genes = unique(table_of_genes())
    print(input$table_of_genes_df_rows_selected)
    gtf_genes[input$table_of_genes_df_rows_selected, ]
  })
  tea_res = reactiveValues(df_res = NULL)
  observeEvent({input$submit_gene_selection}, {
    print(paste0("TEA current tab: ", input$tabs))

    cat(input$submit_gene_selection)
    # print("before ifs")
    if (is.null(input$submit_gene_selection)) return()
    if (!(input$tabs == "qualitycontrol")) return()
    # print("after if")
    req(table_of_genes())
    req(countsfile$df)
    if (is.null(countsfile$df)){
      showModal(modalDialog(size = "l",
                            title = "No countsfile loaded!",
                            "The countsfile has not been loaded so the selected analysis is not possible, yet. 
                            Please proceed here when it has been loaded as indicated in 'Run Overview'. ",
                            easyClose = TRUE,
                            footer = NULL
      ))
    }
    if (is.null(genes.list())){
      print("genes list is null")
      return()
    }
    req(genes.list())
    if (is.null(table_of_normCounts$df_norm)) {
      print("table of norm counts is null")
      return()}
    if (length(table_of_normCounts$df_norm)== 0) {
      print("table of norm counts is 0")
      return()
    }
    req(table_of_normCounts$df_norm)
    print(head(table_of_normCounts$df_norm[[1]]))
    print(head(metadata()))
    print(paste0(Sys.time(), ": Starting TEA"))
    print(genes.list()[,1])
    tea_res$df_res <- TEA(counts = table_of_normCounts$df_norm[[1]],
        norm_counts = table_of_normCounts$df_norm[[2]], 
        genes.list = genes.list(), 
        metadata = table_of_normCounts$df_norm[[3]],
        pvalue = 0.05, 
        output.dir = paste0(table_of_settings_transformed()$run.dir, "TEA"), condi_cols())
  }, priority = 2)
  
  output$table_of_genes_df <- DT::renderDataTable({
    req(table_of_genes())
    config = table_of_genes()
    
    config_out = config
    config_out = unique(config_out)
    DT::datatable(
      config_out,
      extensions = c('Buttons', 'Scroller'),
      options = list(scrollX = 100, scrollY = 400,
                     deferRender = TRUE,
                     scroller = TRUE,
                     paging = TRUE,
                     dom = 'lBfrtip',
                     fixedColumns = TRUE),
      rownames = T)
  })
  
  output$selectedGenes <- DT::renderDataTable({
    req(genes.list())
    config = genes.list()
    config_out = config
    DT::datatable(
      config_out, 
      extensions = c('Buttons', 'Scroller'),
      options = list(scrollX = 100, scrollY = 400,
                     deferRender = TRUE,
                     scroller = TRUE,
                     paging = TRUE,
                     dom = 't',
                     fixedColumns = TRUE),
      rownames = T)
    
  })
  
    output$teaPlot <- renderPlot({
    #if (is.null(tea_res())) return()
    req(tea_res$df_res)
    if (input$disp == "boxplot"){
      return(tea_res$df_res[["Boxplot"]])
    }
    if (input$disp == "dotplot"){
      return(tea_res$df_res[["Dotplot"]])
    }
    if (input$disp == "violinplot"){
      return(tea_res$df_res[["Violinplot"]])
    }
  }, bg="transparent")
  
  # Save button for TEA plot
  observe({
    req(tea_res$df_res)
    if (input$disp == "boxplot"){
      output$down_tea = downloadHandler(
        filename = function() {
          paste0("Gene_counts_boxplot_", Sys.Date(), ".", input$tea_down.type)
        },
        content = function(file) {
          plot <- tea_res$df_res[["Boxplot.down"]]
          ggsave(file, plot, device = input$tea_down.type, bg="white", height = 8, width = 20, units = "in")
        }
      )
    }
    if (input$disp == "dotplot"){
      output$down_tea = downloadHandler(
        filename = function() {
          paste0("Gene_counts_dotplot_", Sys.Date(), ".", input$tea_down.type)
          
        },
        content = function(file) {
          plot <- tea_res$df_res[["Dotplot.down"]]
          ggsave(file, plot, device = input$tea_down.type, bg="white", height = 8, width = 20, units = "in")
        }
      )
    }
    if (input$disp == "violinplot"){
      output$down_tea = downloadHandler(
        filename = function() {
          
          paste0("Gene_counts_violinplot_", Sys.Date(), ".", input$tea_down.type)
          
        },
        content = function(file) {
          plot <- tea_res$df_res[["Violinplot.down"]]
          ggsave(file, plot, device = input$tea_down.type, bg="white", height = 8, width = 20, units = "in")
        }
      )
    }
    
  })
  ## Gene Body Coverage ####

  output$table_of_genes_df_gC <- DT::renderDataTable({
    req(table_of_genes())
    config = table_of_genes()
    
    config_out = config
    config_out = unique(config_out)
    DT::datatable(
      config_out,
      extensions = c('Buttons', 'Scroller'),
      options = list(scrollX = 100, scrollY = 200,
                     deferRender = TRUE,
                     scroller = TRUE,
                     paging = TRUE,
                     dom = 'lBfrtip',
                     fixedColumns = TRUE),
                     selection = list(mode = "single", target = "row"),
      rownames = T)
  })

   genes.list_gC <- reactive({
    req(table_of_genes())
    gtf_genes = unique(table_of_genes())
    print(input$table_of_genes_df_gC_rows_selected)
    gtf_genes[input$table_of_genes_df_gC_rows_selected, ]
  })

  output$selectedGene_gC <- DT::renderDataTable({
    req(genes.list_gC())
    config = genes.list_gC()
    config_out = config
    DT::datatable(
      config_out, 
      extensions = c('Buttons', 'Scroller'),
      options = list(scrollX = 100, scrollY = 200,
                     deferRender = TRUE,
                     scroller = TRUE,
                     paging = TRUE,
                     dom = 'lBfrtip',
                     fixedColumns = TRUE),
      selection = list(mode = "single"
                      , target = "row"),
      rownames = T)
    
  })

  output$submit_gene_selection_gC.out <- renderUI({
    if (!is.null(table_of_genes())){
         actionButton("submit_gene_selection_gC", "Submit gene",
                      class = "btn btn-primary",
                      style='font-size:200%; color: white;
                             background-color: #e76f51; 
                             border-radius: 5px')
      } else {
        actionButton("submit_gene_selection_gC_disabled", "Submit gene",
                     class = "btn btn-primary",
                     style='font-size:200%; color: white;
                            background-color: gray; 
                            border-radius: 5px')
    }
  })
  
  ### Start genebodycoverage pipeline ####
  gB_results = reactiveVal(NULL)
  observeEvent(input$submit_gene_selection_gC, {
    req(table_of_settings_transformed())
    req(input$submit_gene_selection_gC)
    print("HERE")
    if (!file.exists(paste0(docker()$run_dir,"g_percentiles.json"))) return()
    if (!is.null(input$submit_gene_selection_gC)){

      # Save config file
      currentwd = dirname(getwd())
      print(currentwd)
      print(genes.list_gC())
      bamFileList = paste0(docker()$run_dir, "bam_genome_merged/")
      print(bamFileList)
      converted_gtf = paste0(docker()$run_dir, "converted_gtf.csv")
      print(converted_gtf)
      myProcess <<- processx::process$new(command = paste0(currentwd, "/app/server/bash_scripts/run_genebodycoverage.sh"), 
                                          args = c(paste0(currentwd, "/app/server/python_scripts/get_geneBody_coverage.py"), bamFileList, genes.list_gC()$gene_id, converted_gtf, docker()$run_dir), echo_cmd = T,stderr = "error_GC.log", stdout = "index_GC.log")
      myProcess$wait()
      
      df = paste0(docker()$run_dir, "samples.geneBodyCoverage.txt")
      if (!file.exists(df)) print("something went wrong")
      geneBodyCov = read.table(df, header = T)
      print(geneBodyCov)
      print("I'm about to plot")
      geneBodyCov.plots = geneBodyCov.plot(geneBodyCov, genes.list_gC(), metadata(), condi_cols())
      gB_results(geneBodyCov.plots)
      print("GB is done!")
    }
  }, priority = 2)
  
  output$geneBodyCoveragePlot <- renderUI({
    input$submit_gene_selection_gC
    # cat(unlist(unlist(input$submit_gene_selection)))
    if (is.null(input$submit_gene_selection_gC)) return()
    if (input$submit_gene_selection_gC <= 0){
      box(title = "Gene Body Coverage", status = "primary", 
          solidHeader = T, 
          collapsible = T, 
          collapsed = T, 
          width = 12,
          # cat(input$submit_gene_selection),
          fluidRow(
            column(12, radioButtons("disp", "Display",
                                    choices = c("Barplot" = "barplot",
                                                "Dotplot" = "dotplot",
                                                "Line plot" = "lineplot"),
                                    selected = "lineplot", inline = T)),
            column(12, 
                   div(style="display:inline-block; vertical-align:top", downloadButton('down_geneBod')), 
                   div(style="display:inline-block; width:80px; vertical-align:top", selectInput("geneBod_down.type", NULL, choices = c("png", "pdf")))),
            column(12, plotOutput("geneBodyPlot", height = "500px", click = clickOpts(id="plot_click"),
                                  hover = hoverOpts(id = "plot_hover", delayType = "throttle"),
                                  brush = brushOpts(id = "plot_brush"))  %>% withSpinner(color="#0dc5c1"))
          ))
    } else {
      # if (input$submit_gene_selection == 0) return()
      box(title = "Gene Body Coverage", status = "primary", 
          solidHeader = T, 
          collapsible = T, 
          collapsed = F, 
          width = 12,
          # cat(input$submit_gene_selection),
          fluidRow(
            column(12, 
                   div(style="display:inline-block; vertical-align:top", downloadButton('down_geneBod')), 
                   div(style="display:inline-block; width:80px; vertical-align:top", selectInput("geneBod_down.type", NULL, choices = c("png", "pdf")))),
            column(12, plotOutput("geneBodyPlot", height = "500px", click = clickOpts(id="plot_click"),
                                  hover = hoverOpts(id = "plot_hover", delayType = "throttle"),
                                  brush = brushOpts(id = "plot_brush"))  %>% withSpinner(color="#0dc5c1"))
          ))
    }
  })

  output$geneBodyPlot <- renderPlot({   
    req(gB_results())
    ggarrange(gB_results()[[1]], gB_results()[[2]], ncol = 2)
  }, bg="transparent")
  
  # Save button for Coverage plot
  output$down_geneBod = downloadHandler(
    filename = function() {
        paste0("GeneBodyCoverage_",genes.list_gC()$gene_name, "_", Sys.Date(), ".", input$geneBod_down.type)
  
    },
    content = function(file) {
      g1 = gB_results()[[1]] + theme_bw() +
          theme(
        legend.title = element_text(size = 20, color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        axis.text.y = element_text(angle = 45, hjust = 1, size = 17, color = "black"),
        axis.text.x = element_text(size = 17, color = "black"), # remove rotation from x tick labels
        plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "black"),
        axis.title = element_text(size = 23, color = "black"))

      g2 = gB_results()[[2]] + theme_bw() +
          theme(
        legend.title = element_text(size = 20, color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        axis.text.y = element_text(angle = 45, hjust = 1, size = 17, color = "black"),
        axis.text.x = element_text(size = 17, color = "black"), # remove rotation from x tick labels
        plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "black"),
        axis.title = element_text(size = 23, color = "black"))

      g_out = ggarrange(g1, g2, ncol = 2)
     ggsave(file, plot = g_out, device = input$geneBod_down.type, bg="white", height = 10, width = 25, units = "in")
    }
  )


  # Infer experiment #####
  infer_plots <- reactive({
    req(table_of_settings_transformed())
    req(countsfile$df)
    if (is.null(countsfile$df)){
      showModal(modalDialog(size = "l",
                            title = "No countsfile loaded!",
                            "The countsfile has not been loaded so the selected analysis is not possible, yet. 
                            Please proceed here when it has been loaded as indicated in 'Run Overview'. ",
                            easyClose = TRUE,
                            footer = NULL
      ))
    }
    if (!is.null(countsfile$df) & !is.null(metadata())){
      if (is.null(input$design_column) | is.null(input$design_column) | is.null(input$design_column)){
        return()
      }
      p1 = total_genes_counted_plot_per_condition(countsfile$df, metadata())
      p2 = total_genes_counted_plot_per_sample(countsfile$df)
      list(p1, p2)
    } else {
      print(paste0(Sys.time(), ": Analysis not started!"))
      return()
      
    }
  })
  
  
  output$gene_counts_develop <- renderPlot({
    if (is.null(countsfile$df) | is.null(metadata())) return()
    req(countsfile$df)
    req(metadata())
    infer_plots()[[1]]
    
  }, bg="transparent")
  output$gene_counts_develop2 <- renderPlot({
    if (is.null(countsfile$df)) return()
    req(countsfile$df)
    infer_plots()[[2]]
  }, bg="transparent")
  
  
  
  ## 3. DEA ####
  ### Buttons & Actions ####
  #Confirm button for DEA settings
  output$submit_dge <- renderUI({
    if (!(is.null(input$design_column) | is.null(input$feature_A) | is.null(input$feature_B))){
      actionButton(inputId = "run_dea_B" , label = "Gene Expression Analysis", class = "btn btn-primary",
                   style = "white-space: normal;display:inline-block;width:50%;text-align: center;
                                                        color: #fff; height: 70pt;
                                                        font-size:180%;
                                                        background-color: #e76f51;")
    } else {
      actionButton(inputId = "run_dea_B_disabled" , label = "Gene Expression Analysis", class = "btn btn-primary",
                   style = "white-space: normal;display:inline-block;width:50%;text-align: center;
                                                        color: #fff; height: 70pt;
                                                        font-size:180%;
                                                        background-color: gray;")
    }
  })
  
  output$submit_dt_preprocess <- renderUI({
    if (!(is.null(input$design_column) | is.null(input$feature_A) | is.null(input$feature_B))){
      actionButton(inputId = "run_dt_preprocess" , label = "Transcript Expression Analysis", class = "btn btn-primary",
                   style = "white-space: normal;display:inline-block;width:50%;text-align: center;
                                                        color: #fff; height: 70pt;
                                                        font-size:180%;
                                                        background-color: #e76f51;")
    } else {
      actionButton(inputId = "run_dt_preprocess_disabled" , label = "Transcript Expression Analysis", class = "btn btn-primary",
                   style = "white-space: normal;display:inline-block;width:50%;text-align: center;
                                                        color: #fff; height: 70pt;
                                                        font-size:180%;
                                                        background-color: gray;")
      
    }
  })

  #output$run_dte <- renderUI({
  #  if (!(is.null(input$design_column) | is.null(input$feature_A) | is.null(input$feature_B))){
  #    actionButton(inputId = "run_dte_B" , label = "Transcript Expression", icon("arrow-alt-circle-right"), class = "btn btn-primary",
  #                 style='font-size:200%; color: white; background-color: #669bbc; border-radius: 10px')
  #  } else {
  #    disabled(actionButton(inputId = "run_dte_B" , label = "Transcript Expression", icon("arrow-alt-circle-right"), class = "btn btn-primary",
  #                          style='font-size:200%; color: white; background-color: #669bbc; border-radius: 10px'))
  #    
  #  }
  #})
    
  #output$submit_dte <- renderUI({
  #  if (!(is.null(input$design_column) | is.null(input$feature_A) | is.null(input$feature_B))){
  #    actionButton(inputId = "run_dte_B" , label = "Transcript Expression", icon("arrow-alt-circle-right"), class = "btn btn-primary",
  #                 style='font-size:200%; color: white; background-color: #669bbc; border-radius: 10px')
  #  } else {
  #    disabled(actionButton(inputId = "run_dte_B" , label = "Transcript Expression", icon("arrow-alt-circle-right"), class = "btn btn-primary",
  #                          style='font-size:200%; color: white; background-color: #669bbc; border-radius: 10px'))
  #    
  #  }
  #})
  
  output$submit_dtu <- renderUI({
    if (!(is.null(input$design_column) | is.null(input$feature_A) | is.null(input$feature_B))){
      actionButton(inputId = "run_dtu_B" , label = "Submit gene selection", class = "btn btn-primary",
                      style='font-size:200%; color: white;
                             background-color: #e76f51; 
                             border-radius: 5px')
    } else {
      actionButton(inputId = "run_dtu_B_disabled" , label = "Submit gene selection", class = "btn btn-primary",
                     style='font-size:200%; color: white;
                            background-color: gray; 
                            border-radius: 5px')
    }
  })
  
  observeEvent(input$run_dea_B, {
    updateTabsetPanel(session, "dea.box",
                      selected = "dea.tab")
  }, priority = 5)
  
  
  observeEvent(input$run_dt_preprocess, {
    updateTabsetPanel(session, "dea.box",
                      selected = "deu.tab")
  }, priority = 5)

  #observeEvent(input$run_dte_B, {
  #  updateTabsetPanel(session, "dea.box",
  #                    selected = "dte.tab")
  #}, priority = 5)
  
  
  
  #observeEvent(input$jump_to_settings_overview_B, {
  #  if (is.null(metadata())) return()
  #  if (!is.null(input$design_column)){
  #    if (!(input$design_column == "")){
  #      output$design_column_in2 <- renderUI({
  #        selectInput(inputId = 'design_column', 
  #                    label = 'Select condition column', 
  #                    choices = colnames(metadata()),selected = input$design_column)#,
  #                    #options = list(`actions-box` = TRUE))
  #      })
  #      output$feature_column_A2 <- renderUI({
  #        selectInput(inputId = 'feature_A2', 
  #                    label = 'Select condition A', 
  #                    choices = unique(metadata()[,input$design_column]), selected = input$feature_A)#, 
  #                    #options = list(`actions-box` = TRUE))
  #      })
        
  #      output$feature_column_B2 <- renderUI({
  #        if (is.null(metadata())) return()
  #        x = unique(metadata()[,input$design_column])
  #        selectInput(inputId = 'feature_B2', 
  #                    label = 'Select condition B', 
  #                    choices = x[!(x %in% input$feature_A)], selected = input$feature_B)#,
  #                    #options = list(`actions-box` = TRUE))
  #      })
  #      
  #    }
  #  }
  #})

  
  # Dialog is only run if table_of_genes is loaded
  observeEvent(input$run_dea_B,{
    req(table_of_genes())
    showModal(modalDialog(size = "l",
                          title = "Differential Expression Analysis started",
                          "The DEA has started!",
                          easyClose = TRUE,
                          footer = NULL
    ))
    
    if (input$run_dea_B == 1){
      print(head(table_of_genes()))
    }
  })
  ### Select Colors####
    
  output$select_color.heat <- renderUI({
    fluidRow(
      box(
      column(12, div(style = "font-size: 13px", selectInput("main_color.heat", "Select Colorpalette", choices = c("RdBu", "PiYG", "PRGn", "PuOr", "Blue-Yellow",
                                                                                                                  "Teal", "Sunset", "Viridis")))),
      #column(6, div(style = "font-size: 13px", radioButtons("color_one.heat", "First Group Color", choices = c("#88CCEE", "#DDCC77", "#CC6677", "#117733", "#332288", "#AA4499", 
      #                                                                                                         "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")))), 
      #column(6, div(style = "font-size: 13px", radioButtons("color_two.heat", "Second Group Color", choices = c("#88CCEE", "#DDCC77", "#CC6677", "#117733", "#332288", "#AA4499", 
      #                                                                                                          "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888"), selected = "#DDCC77"))), 
      width = 12, solidHeader = T, collapsible = T, collapsed = T, status = "primary", title = div("Customize Colors", style = "font-size: 16px")
    )
    )
  })

  ######## TRANSCRIPT ANALYSIS #########
  preProcTrans <- reactiveValues(pre_list = NULL)
  observeEvent({input$run_dt_preprocess}, {
    req(countsfile$df_trans)
    print(input$tabs)
    print(input$dea.box)
    if (input$run_dt_preprocess > 0){
      preProcTrans$pre_list <- DRIM_seq_prep(
                                    table = countsfile$df_trans,
                                    run.dir = docker()$run_dir, 
                                    condition_col = input$design_column, 
                                    first.level = input$feature_A, 
                                    ref.level = input$feature_B, 
                                    samps = metadata(), 
                                    gtf_file = docker()$gtf_file,
                                    cores = 4
      )
      print(names(preProcTrans$pre_list))
    } else {
      print(paste0(Sys.time(), ": Analysis not started!"))
      
    }
    }, priority  = 0)
  

  #############DTE Analysis##########################
  
  DTE_run <- reactiveValues(df_res_dte = NULL)
  
  observeEvent({input$run_dt_preprocess}, {
    req(preProcTrans$pre_list)
    if (is.null(preProcTrans$pre_list)) return()
    gtf_csv = paste0(docker()$run_dir, "converted_gtf.csv")
    gtf_table <- read.table(gtf_csv, ",",header = T) 
    print("Until here it works!")
    print(names(preProcTrans$pre_list))
 
    DTE_run$df_res_dte <- DTE_general(
      d_list = preProcTrans$pre_list, 
                                  condition_col = input$design_column, 
                                  first.level = input$feature_A, 
                                  ref.level = input$feature_B, 
                                  samps = metadata(), 
                                  gtf_table = gtf_table,
                                  cores = 4
    )
    #save_rds(DTE_run$df_res_dte, table_of_settings_transformed()$run.dir)
    
  })
    
    
  
  ##################################################
  
  
  ###########DTU Analysis##########################  
  DTU_run <- reactiveValues(df_res_dtu = NULL)

  observeEvent(input$run_dtu_B, {
    print("NOW I AM IN DTU ANALYSIS")
    req(preProcTrans$pre_list)
    print(input$tabs)
    if (is.null(preProcTrans$pre_list)) return()

    s = input$dtu_tab_rows_selected
    print("Input DTU rows")
    print(s)
    if (is.null(table_of_genes()$gene_id[input$dtu_tab_rows_selected])){
      print("Please select gene of interest in datatable")
    } else{
      gtf_csv = paste0(docker()$run_dir, "converted_gtf.csv")
      gtf_table <- read.table(gtf_csv, ",",header = T) 

      DTU_run$df_res_dtu <- DTU_special(
        d_list = preProcTrans$pre_list,  
                                        condition_col = input$design_column, 
                                        first.level = input$feature_A, 
                                        ref.level = input$feature_B,
                                        goi_id = table_of_genes()$gene_id[s],
                                        gtf_tab = gtf_table,
                                        cores = 4
      )
      
      #save_rds(DTU_run$df_res_dtu, table_of_settings_transformed()$run.dir)
    }
  })
  
  
  
  
  
  #################################################
  dea_res_preprocess <- reactiveValues(df_res = NULL)
  
  observeEvent({input$run_dea_B}, {
    req(table_of_settings_transformed())
    req(countsfile$df)
    req(table_of_genes())
    req(docker()$gtf_file)

    if (is.null(countsfile$df)){
      showModal(modalDialog(size = "l",
                            title = "No countsfile loaded!",
                            "The countsfile has not been loaded so the selected analysis is not possible, yet. 
                            Please proceed here when it has been loaded as indicated in 'Run Overview'. ",
                            easyClose = TRUE,
                            footer = NULL
      ))
    }
    if (input$tabs == "dea_main" & !(is.null(countsfile$df)) & !is.null(metadata())){
      
      cat(unlist(metadata()))
      dea_res_preprocess$df_res <- run_preprocessing_dea(meta.file = metadata(),
                            counts.file = countsfile$df, 
                            condition.col = input$design_column,
                            first.level = input$feature_A, 
                            ref.level = input$feature_B, 
                            pvalue = as.numeric(input$pvalue), 
                            gtf_file = table_of_genes()
                            )
      save_rds(dea_res_preprocess$df_res, table_of_settings_transformed()$run.dir)
    } else {
      print(paste0(Sys.time(), ": Analysis not started!"))
      
    }

  })
  
  dea_res <- reactive({
    req(table_of_settings_transformed())
    req(table_of_genes())
    req(dea_res_preprocess$df_res)    
    run_dea(metadata = dea_res_preprocess$df_res$metadata, 
            counts = dea_res_preprocess$df_res$counts,
            condition.col = input$design_column, 
            first.level = input$feature_A, 
            ref.level = input$feature_B, 
            rld = dea_res_preprocess$df_res$rld, 
            dds = dea_res_preprocess$df_res$dds, 
            res_df = dea_res_preprocess$df_res$res_df,
            # pvalue = as.numeric(input$pvalue), 
            # output.dir = table_of_settings_transformed()$run.dir, 
            gtf_file = table_of_genes(), 
            condi_cols = condi_cols()
    )
  })
  
  output$text <- renderText({
    req(dea_res())
    input$preprocessing_B
    if (length(dea_res()) == 0){
      "No file loaded!"
    } else {
      paste0(Sys.time(), " ", paste0(names(dea_res()), collapse = ","))
      
    }
  })
  
  ### Download buttons ####
  #!!!!!!MUCH ROOM FOR IMPROVEMENT!!!!!!!
  # Save button for PCA plot
  output$down_pca = downloadHandler(
    
    filename = function() {
        paste0("pca_", Sys.Date(), ".", input$down_pca.type)
    },
    content = function(file) {
      req(dea_res_preprocess$df_res)
      plot <- createPCA(dea_res_preprocess$df_res$rld, 
              first.level = input$feature_A, 
              ref.level = input$feature_B, 
              condi_cols()) + theme_bw() +
              theme(legend.title = element_text(size = 20),
                legend.text = element_text(size = 20),
                axis.text = element_text(angle = 45, hjust = 1, size = 17),
                plot.title = element_text(hjust = 0.5, face = "bold", size = 23),
                axis.title = element_text(size = 23))

      ggsave(file, plot = plot, device = input$down_pca.type, width = 10, height = 10, bg = "white")
    }
  )
  
  # Save button for volcano plot
  output$down_vol = downloadHandler(
    
    filename = function() {
      paste0("DEA_volcano_plot_", Sys.Date(), ".", input$down_vol.type)
    },
    content = function(file) {
      req(dea_res_preprocess$df_res)
      plot <- createVolcano(dea_res_preprocess$df_res$res_df, condi_cols()) + theme_bw() +
              theme(#legend.title = element_text(size = 20),
                legend.title = element_blank(), # remove legend title
                legend.text = element_text(size = 20),
                axis.text = element_text(angle = 45, hjust = 1, size = 17),
                plot.title = element_text(hjust = 0.5, face = "bold", size = 23),
                axis.title = element_text(size = 23))

      ggsave(file, plot = plot, device = input$down_vol.type, bg = "white", width = 10, height = 7)
    }
  )
  
  # Save button for Sam2sam plot
  output$down_s2s = downloadHandler(
    
    filename = function() {
      paste0("DEA_sam2sam_", Sys.Date(), ".", input$down_s2s.type)
    },
    content = function(file) {
      req(dea_res_preprocess$df_res)

      if (input$down_s2s.type == "pdf"){
        pdf(file = file, width = 10, height = 10)
      } 
      if (input$down_s2s.type == "png"){
        png(file = file, width = 10, height = 20, units = "cm", res = 300)
      }
      draw(createSam2Sam(dea_res_preprocess$df_res$rld)[[2]])
      dev.off()
      # save_rds(dea_res(), table_of_settings_transformed()$run.dir)
      
    }
  )
  
  # Save button for heatmap
  output$down_heat = downloadHandler(
    
    filename = function() {
      paste0("DEA_expression_heatmap_", Sys.Date(), ".", input$down_heat.type)
    },
    content = function(file) {
      req(dea_res_preprocess$df_res)
      req(input$main_color.heat)
      color_plot <- condi_cols()
      plot <- createHeatmap(dea_res_preprocess$df_res$dds,
                           dea_res_preprocess$df_res$rld, 
                           condi_col = color_plot, 
                           genes = dea_res_preprocess$df_res$res_df[1:20, ], 
                           main_color = input$main_color.heat)
      if (input$down_heat.type == "pdf"){
        pdf(file = file,width= 20, height = 20)
      } 
      if (input$down_heat.type == "png"){
        png(file = file,width= 20, height = 20, units = "cm", res = 300)
      }
      print(plot$heat.down)
      dev.off()
    }
  )

  ## Download button for DEX Volcano
  output$down_vol_dex = downloadHandler(
    
    filename = function() {
      paste0("DEX_volcano_plot_", Sys.Date(), ".", input$down_vol_dex.type)
    },
    content = function(file) {
      req(DTE_run$df_res_dte$volcano_plot)
      plot = DTE_run$df_res_dte$volcano_plot + 
        theme_bw() +
        theme(legend.title = element_blank(),
                legend.text = element_text(size = 20),
                axis.text = element_text(angle = 45, hjust = 1, size = 17),
                plot.title = element_text(hjust = 0.5, face = "bold", size = 23),
                axis.title = element_text(size = 23))

      ggsave(file, plot = plot, device = input$down_vol_dex.type, bg = "white", width = 10, height = 7)
      # save_rds(dea_res(), table_of_settings_transformed()$run.dir)
      
    }
  )

    ## Download button for DEX Volcano
  output$down_dtu_boxplot = downloadHandler(
    
    filename = function() {
      req(input$dtu_tab_rows_selected)
      paste0("DTU_", table_of_genes()$gene_id[input$dtu_tab_rows_selected], "_", Sys.Date(), ".", input$down_dtu_boxplot.type)
    },
    content = function(file) {
      req(DTU_run$df_res_dtu$bp)
      plot = DTU_run$df_res_dtu$bp + 
        theme_bw() +
        theme(legend.title = element_text(size = 20),
                legend.text = element_text(size = 20),
                axis.text = element_text(angle = 45, hjust = 1, size = 8),
                plot.title = element_text(hjust = 0.5, face = "bold", size = 23),
                axis.title = element_text(size = 23))

      ggsave(file, plot = plot, device = input$down_dtu_boxplot.type, bg = "white", width = 10, height = 7)
      # save_rds(dea_res(), table_of_settings_transformed()$run.dir)
      
    }
  )
  
  ## Sequencing run table ####
  seqRunInfos = reactiveValues(df = NULL) 
  
  observeEvent(input$preprocessing_B,ignoreNULL = T,{
              seqRunInfos$df <- reactivePoll(10000, session,
                checkFunc = function(){
                  req(docker()$run_dir)
                  mapStats=paste0(docker()$run_dir, "mapping_stats.txt")
                  if (file.exists(mapStats)){   
                    file.info(mapStats)$mtime[1]
                  } else {
                    print("Not processed yet")
                  }
                },
                valueFunc = function(){
                  
                  df_out = data.frame("Samples" = metadata()$Samples, "mapped_reads" = 0, "gene_counts" = 0, "transcript_counts" = 0)
                  print(df_out)

                  # Num Q reads and mapped reads
                  if (file.exists(paste0(docker()$run_dir, "mapping_stats.txt"))){
                    print("Mapping statistics exists.")
                    mapStats = read.table(paste0(docker()$run_dir, "mapping_stats.txt"), header = T, sep = "\t")
                    print(head(mapStats))
                    colnames(mapStats) = c("Samples","mapped_reads", "X")
                    mapStats = mapStats[match(df_out$Samples, mapStats$Samples),]
                    df_out$mapped_reads = mapStats$mapped_reads
                  }

                  ## Num gene counts 
                  if (!is.null(countsfile$df)){
                    print(head(countsfile$df))
                    print(colSums(countsfile$df))
                    df_out$gene_counts =  as.integer(colSums(countsfile$df[c(df_out$Samples)]))
                  } else {
                    print("Counts file not present")
                  }
                  
                  ## Num gene counts 
                  if (!is.null(countsfile$df_trans)){
                    print(head(countsfile$df_trans))
                    print(colSums(countsfile$df_trans))
                    print(as.integer(colSums(countsfile$df_trans)[c(df_out$Samples)]))
                    df_out$transcript_counts =  as.integer(colSums(countsfile$df_trans)[c(df_out$Samples)])
                  } else {
                    print("Transcript counts file not present")
                  }
                  df_out
                })
               },
               priority = -1
  )


  ### Tables and Plots ####
  # Show table of seq. run infos
  output$table_seq_run <- renderDT({
    input$preprocessing_B
    print(seqRunInfos$df())
    #req(seqRunInfos)
    DT::datatable(
      seqRunInfos$df(),
      extensions = c('Buttons', 'Scroller'),
      options = list(scrollY = 200,
                     scrollX = 500,
                     deferRender = TRUE,
                     scroller = TRUE,
                     paging = TRUE,
                     dom = 't',
                     fixedColumns = TRUE), 
      rownames = F)
  })


  process_time = reactiveValues(df = NULL)
  
  observeEvent(input$preprocessing_B,ignoreNULL = T,{
              process_time$df <- reactivePoll(10000, session,
                checkFunc = function(){
                  req(docker()$run_dir)
                  if (file.exists(paste0(docker()$run_dir, "processing_time_table.csv"))){
                    file.info(paste0(docker()$run_dir, "processing_time_table.csv"))$mtime[1]
                  }
                },
                valueFunc = function(){
                  if (file.exists(paste0(docker()$run_dir, "processing_time_table.csv"))){
                    read.table(paste0(docker()$run_dir, "processing_time_table.csv"), header = T, sep = ",")
                  } else {
                    NULL
                  }
                })
               },
               priority = -1
  )
  ## Plot Process time per Iteration
   output$process_time.out <- renderPlot({
    req(process_time$df())
    if (!is.null(process_time$df())){
      print(head(process_time$df()))
      x = process_time$df()
      # replace _ with whitespace in tool names for a cleaner legend
      x = x %>% mutate(Tool = gsub("_", " ", Tool))
      
      # get the x-tick interval based on the number of iterations in the data
      max_it <- max(x$Iteration)
      xticks <- 1:max_it
      if (max_it <= 10) { # every tick
        xticks <- xticks
      } else if (max_it <= 50) { # every 5th tick
        xticks <- xticks[xticks %% 5 == 0]
      } else if (max_it <= 100) { # every 10th tick
        xticks <- xticks[xticks %% 10 == 0]
      } else if (max_it <= 1000) { # every 50th tick
        xticks <- xticks[xticks %% 50 == 0]
      } else if (max_it <= 10000) { # every 100th tick
        xticks <- xticks[xticks %% 100 == 0]
      } else { # every 1000th tick
        xticks <- xticks[xticks %% 1000 == 0]
      }

      ggplot(x, aes(x = Iteration, y = Time, fill = Tool)) +
        geom_bar(stat = "identity") +
        scale_x_continuous(breaks = xticks) + # label each x tick individually
        ylab("Time [s]") + # add [s] to y-axis label
        theme_bw() +
         scale_fill_manual("Preprocessing steps", values = c("#CC6677", "#d8c66c", "#117733", "#88CCEE", "#AA4499", "#24011e","#092a36")) +
        theme(
        panel.grid.minor = element_blank(), # remove minor grid lines
        panel.grid.major.x = element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        # panel.grid.major = element_blank(), # get rid of major grid
        # panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        #legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
        legend.title = element_text("", size = 20, color = "white"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.text = element_text(size = 20, color = "white"),
        axis.text = element_text(angle = 45, hjust = 1, size = 17, color = "white"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 23, color = "white"),
        axis.title = element_text(size = 23, color = "white"))

    }
  }, bg="transparent")
  
  
  
  # Show table of DEGs 
  output$degs_tab <- renderDT({
    input$preprocessing_B
    req(dea_res_preprocess$df_res)
    res = dea_res_preprocess$df_res$res_df
    degs = res[res$padj < 0.05,]
    DT::datatable(
      degs,
      extensions = c('Buttons', 'Scroller'),
      options = list(scrollY = 200,
                     scrollX = 500,
                     deferRender = TRUE,
                     scroller = TRUE,
                     paging = TRUE,
                     dom = 'lBfrtip',
                     fixedColumns = TRUE), 
      rownames = F)
  })
  
  # Show table of DETs
  output$dte_tab <- renderDT({
    input$preprocessing_B
    req(DTE_run$df_res_dte$dxr_df)
    res = DTE_run$df_res_dte$dxr_df
    #print(res)
    dtes = res[res$padj < 0.05,]
    dtes = data.frame(gene_name = as.character(dtes$gene_name),group_id = as.character(dtes$groupID), feature_id = as.character(dtes$featureID),exon_base_mean = as.numeric(as.character(dtes$exonBaseMean)),p_adjusted = as.numeric(as.character(dtes$padj)), dispersion = as.numeric(as.character(dtes$dispersion)))
    colnames(dtes) = c("names","gencode","feature ID", "Exon base mean", "p-adjusted", "dispersion")
    DT::datatable(
      dtes,
      extensions = c('Buttons', 'Scroller'),
      options = list(scrollY = 200,
                     scrollX = 500,
                     deferRender = TRUE,
                     scroller = TRUE,
                     paging = TRUE,
                     dom = 'lBfrtip',
                     fixedColumns = TRUE), 
      rownames = F)
  })
  
  output$dtu_tab <- DT::renderDataTable({
    input$preprocessing_B
    #s = input$dtu_tab_row_selected
    req(table_of_genes())
    config = table_of_genes()
    
    config_out = config
    DT::datatable(
      config_out,
      extensions = c('Buttons', 'Scroller'),
      options = list(scrollX = 100, scrollY = 400,
                     deferRender = TRUE,
                     pageLength = 5,
                     scroller = TRUE,
                     paging = TRUE,
                     dom = 'lBfrtip',
                     fixedColumns = TRUE),
      selection = list(mode = "single"
                       , target = "row"),
      rownames = T)
  })
  
  output$selectedGene_dtu <- DT::renderDataTable({
    req(table_of_genes())
    s = table_of_genes()[input$dtu_tab_rows_selected,]
    config = s
    config_out = config
    DT::datatable(
      config_out, 
      extensions = c('Buttons', 'Scroller'),
      options = list(scrollX = 100, scrollY = 400,
                     deferRender = TRUE,
                     scroller = TRUE,
                     paging = TRUE,
                     dom = 't',
                     fixedColumns = TRUE),
      rownames = T)
    
  })
  
  # Section One ---------------------------------
  
  # Plot volcano plot stored in RDS object
  output$volcano_plot <- renderPlot({
    req(dea_res_preprocess$df_res)
    createVolcano(dea_res_preprocess$df_res$res_df, condi_cols())
  }, bg="transparent")
  
  
  output$pca_plot <- renderPlot({
    print("PCA is made")
    req(dea_res_preprocess$df_res)
    createPCA(dea_res_preprocess$df_res$rld, 
              first.level = input$feature_A, 
              ref.level = input$feature_B, 
              condi_cols())   
  }, bg="transparent")
  
  output$s2s_plot <- renderPlot(bg = "transparent",{
    print("Heatmap")

    req(dea_res_preprocess$df_res)
    createSam2Sam(dea_res_preprocess$df_res$rld)[[1]]
    
  })
  
  output$heat_plot <- renderPlot(bg = "transparent", {
    print("Heatmap")
    req(dea_res_preprocess$df_res)
    req(input$main_color.heat)
    # does not load correctly after being downloaded, all values still exist but nothing is rendered
    # the plot can be seen when printed in RStudio, though  
    plot = createHeatmap(dea_res_preprocess$df_res$dds,
                  dea_res_preprocess$df_res$rld, 
                  condi_col = condi_cols(), 
                  genes = dea_res_preprocess$df_res$res_df[1:20, ], 
                  main_color = input$main_color.heat)
    plot$heat
  })
  
  output$volcano_dex_output <- renderPlot(bg = "transparent",{
    print("DEX Volcano plot")
    req(DTE_run$df_res_dte$volcano_plot)
    DTE_run$df_res_dte$volcano_plot
    
  })
  
  output$dtu_boxplot <- renderPlot(bg = "transparent",{
    print("DTU Boxplot")
    req(DTU_run$df_res_dtu)
    DTU_run$df_res_dtu$bp
  })
    
  output$dtu_boxplot_UI <- renderUI({
      if (!is.null(DTU_run$df_res_dtu)) {
        if (is.null(DTU_run$df_res_dtu$bp)) {
          h3(p(em(strong("No transcripts found!"))))
        } else {
          column(12, 
          div(style = "display: inline-block; vertical-align: top", downloadButton("down_dtu_boxplot")),
          div(style = "display: inline-block; vertical-align: top; width: 80px", selectInput("down_dtu_boxplot.type", NULL, choices = c("png", "pdf"))),
          plotOutput("dtu_boxplot") %>% withSpinner(color = "#0dc5c1"))
        }
      }
  })

  output$dea_results.out <- renderUI({
    if (is.null(input$run_dea_B)) return()
    if(input$run_dea_B > 0){
      fluidRow(
          column(12, box(
            title = "Differential expressed genes",
            width = 12,
            status = "primary",
            solidHeader = T,
            collapsible = T,
            collapsed = F,
            DT::dataTableOutput("degs_tab") %>% withSpinner(color = "#0dc5c1")
          )),
          column(
            12,
            tabBox(
              width = 12, title = "Plots",
              # title = "Quality Box",
              # The id lets us use input$tabset1 on the server to find the current tab
              id = "tabset.plots",
              tabPanel(
                title = "PCA", value = "pca.tab",
                fluidRow(column(
                  12,
                  div(style = "display:inline-block; vertical-align: top", downloadButton("down_pca")),
                  div(style = "display:inline-block; vertical-align: top; width: 80px", selectInput("down_pca.type", NULL, choices = c("png", "pdf"))),
                  # div(style = "display:inline-block; vertical-align: top; margin-top: 8px", radioButtons("custom_color.pca", NULL, choices = c("default colors" = F, "custom colors" = T), inline = T)),
                  #uiOutput("select_color.pca"),
                  plotOutput("pca_plot") %>% withSpinner(color = "#0dc5c1")
                ))
              ),
              tabPanel(
                title = "Volcano", value = "volcano.tab",
                fluidRow(column(
                  12,
                  div(style = "display: inline-block; vertical-align: top", downloadButton("down_vol")),
                  div(style = "display: inline-block; vertical-align: top; width: 80px", selectInput("down_vol.type", NULL, choices = c("png", "pdf"))),
                  # div(style = "display:inline-block; vertical-align: top; margin-top: 8px", radioButtons("custom_color.volcano", NULL, choices = c("default colors" = F, "custom colors" = T), inline = T)),
                  uiOutput("select_color.volcano"),
                  plotOutput("volcano_plot") %>% withSpinner(color = "#0dc5c1")
                ))
              ),
              tabPanel(
                title = "Sample2Sample", value = "s2s.tab",
                fluidRow(column(
                  12,
                  div(style = "display:inline-block; vertical-align: top", downloadButton("down_s2s")),
                  div(style = "display:inline-block; vertical-align: top; width: 80px", selectInput("down_s2s.type", NULL, choices = c("png", "pdf"))),
                  plotOutput("s2s_plot") %>% withSpinner(color = "#0dc5c1")
                ))
              ),
              tabPanel(
                title = "Heatmap", value = "heatmap.tab",
                fluidRow(column(
                  12,
                  div(style = "display:inline-block; vertical-align:top", downloadButton("down_heat")),
                  div(style = "display:inline-block; vertical-align:top; width:80px", selectInput("down_heat.type", NULL, choices = c("png", "pdf"))),
                  uiOutput("select_color.heat"),
                  plotOutput("heat_plot") %>% withSpinner(color = "#0dc5c1")
                ))
              )
            ) # tabBox Plots
          )
        ) # fluidRowcolumn
    }
  })

  ##Stop preprocessing
  observeEvent(input$preprocessing_B,{
                 process_running$x <- reactivePoll(10000, session,
                                  checkFunc = function(){
                                    req(docker()$run_dir)
                                    process_run = paste0(docker()$run_dir, "process_running.txt")
                                    if (dir.exists(docker()$run_dir)){
                                      if (file.exists(process_run)){
                                        file.info(process_run)$mtime[1]
                                      }
                                    }
                                  },
                                  valueFunc = function(){
                                    req(docker()$run_dir)
                                    process_run = paste0(docker()$run_dir, "process_running.txt")
                                    if (dir.exists(docker()$run_dir)){
                                      if (file.exists(process_run)){
                                        print("Reading process x")
                                        print(read.table(process_run, header = F)[1,1])
                                        read.table(process_run, header = F)[1,1]
                                      }
                                    }
                                }
                              )
                 }
  )


  observeEvent(input$preprocessing_B, {
    req(process_running$x())
      print("Here process running: ")
      print(process_running$x())
      if(process_running$x() == 1){
              output$stop_preprocessing <- renderUI({
                        actionButton(inputId = "stop_process" , label = "Stop preprocessing", class = "btn btn-primary",
                        style='font-size:100%; color: white;
                        background-color: #e76f51; 
                        border-radius: 5px')
                  })
              output$resume_preprocessing <- renderUI({
                    actionButton(inputId = "resume_process_disabled" , label = "Resume preprocessing", class = "btn btn-primary",
                    style='font-size:100%; color: white;
                            background-color: gray; 
                            border-radius: 5px')
              })
        } 

  })
  

  observeEvent(input$stop_process, {
    req(input$preprocessing_B)
    req(input$stop_process)
    req(process_running$x())
    write.table(c(0), file = paste0(docker()$run_dir, "process_running.txt"), col.names = F, row.names = F, quote = F)
    print("Nextflow pipeline is stopped!")
    output$stop_preprocessing <- renderUI({ 
      
      actionButton(inputId = "stop_process_disabled" , label = "Stop preprocessing", class = "btn btn-primary",
                     style='font-size:100%; color: white;
                            background-color: gray; 
                            border-radius: 5px')
      })
    output$resume_preprocessing <- renderUI({
      actionButton(inputId = "resume_process_disabled" , label = "Resume preprocessing", class = "btn btn-primary",
                     style='font-size:100%; color: white;
                            background-color: gray; 
                            border-radius: 5px')
      })
  }, priority = 0)

### Resume preprocessing
  autoInvalidate <- reactiveTimer(4000)
  observe({
  autoInvalidate()
    req(input$preprocessing_B)
    req(process_running$x())
    if (file.exists(paste0(docker()$run_dir, "process_running.txt"))){
      restart = read.table(paste0(docker()$run_dir, "process_running.txt"), header = F)[1,1]
    }
    if (!is.null(input$preprocessing_B) & as.integer(restart) == 2){
      if (input$preprocessing_B > 0){
        output$resume_preprocessing <- renderUI({
          actionButton(inputId = "resume_process" , label = "Resume preprocessing", class = "btn btn-primary",
                                style='font-size:100%; color: white;
                                      background-color: #e76f51; 
                                      border-radius: 5px')
        })
        }
      }
  })

  observeEvent(input$resume_process, {
    write.table(c(1), file = paste0(docker()$run_dir, "process_running.txt"), col.names = F, row.names = F, quote = F)
    print("Nextflow pipeline is resumed!")
    output$resume_preprocessing <- renderUI({
      actionButton(inputId = "resume_process_disabled" , label = "Resume preprocessing", class = "btn btn-primary",
                     style='font-size:100%; color: white;
                            background-color: gray; 
                            border-radius: 5px')
        }) 
      output$stop_preprocessing <- renderUI({ 
      actionButton(inputId = "stop_process" , label = "Stop preprocessing", class = "btn btn-primary",
                                style='font-size:100%; color: white;
                                      background-color: #e76f51; 
                                      border-radius: 5px')
        })   
  })


  
  # Kill all processes when app is closed
  onStop(function(){
    cat(sprintf("Session %s was closed\n", session$token))
    if(!is.null(myProcess)){
      if(myProcess$is_alive()){
        myProcess$kill()
        unlink("work", recursive = TRUE) 
      }
    }
    
  })
  
}
