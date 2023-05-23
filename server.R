extract_which_chr <- function(data_in){
  unique(data_in$CHR)[order(match(unique(data_in$CHR), c(paste0("",1:22), "X", "Y")))]
}

get_cumulative_length <- function(chrom_lengths){
  cumulative_length <- 0
  
  if(length(chrom_lengths) > 1){
    cumulative_length <- head(c(0,cumsum(x = unname(chrom_lengths))), -1)
    names(cumulative_length) <- names(chrom_lengths)
  }
  return(cumulative_length)
}

get_x_breaks <- function(chrom_lengths){
  cumulative_length <- get_cumulative_length(chrom_lengths)
  x_breaks <-cumulative_length+round(chrom_lengths/2)
  names(x_breaks)=gsub('chr', '', names(x_breaks))
  if(length(chrom_lengths) == 21){
    names(x_breaks)[20]=''
  }
  if(length(chrom_lengths) > 21){
    names(x_breaks)[20]=''
    names(x_breaks)[22]=''
  }
  return(x_breaks)
}

add_cumulative_pos <- function(data_in, chrom_lengths){
  cumulative_length <- get_cumulative_length(chrom_lengths)
  
  tmp <- Map(function(x,y){x$cumulative_pos <- x$POS + y; return(x)}, 
             split(data_in, data_in$CHR)[names(chrom_lengths)], get_cumulative_length(chrom_lengths)) 
  return(do.call(rbind, tmp))
}  

add_color <- function(data_in, color1='black', color2='grey'){
  if ('color'%in%colnames(data_in)){
    user_color <- data_in$color
  } else {
    user_color <- rep(NA, nrow(data_in))
  }
  data_in$color <- ifelse(data_in$CHR %in% extract_which_chr(data_in)[c(TRUE, FALSE)], color1, color2)
  data_in$color <- ifelse(is.na(user_color), data_in$color, user_color)
  return(data_in)
}


server <- function(input, output, session) {

  # UI - GENERAL ---------------------------------------------------------------
  data <- reactiveValues(df = NULL)
  
  # READ INPUT FILE ------------------------------------------------------------
  observeEvent(input$inputFile, {
    req(input$inputFile)
    
    # Read File
    df <- fread(input$inputFile$datapath, sep = "\t", header = T, data.table = F)
    
    # Convert P-values
    df$P_Value_Marginal    <- -log10(df$P_Value_Marginal)
    df$P_Value_Interaction <- -log10(df$P_Value_Interaction)
    df$P_Value_Joint       <- -log10(df$P_Value_Joint)
    df$robust_P_Value_Marginal    <- -log10(df$robust_P_Value_Marginal)
    df$robust_P_Value_Interaction <- -log10(df$robust_P_Value_Interaction)
    df$robust_P_Value_Joint       <- -log10(df$robust_P_Value_Joint)
    
    # Get interactions including marginal columns
    interactions <- c("Beta_Marginal", "Beta_G", colnames(df)[grepl("^Beta_G-", colnames(df))])
    interactions <- gsub("Beta_", "", interactions)
    int_colnames <- c("Marginal", "Main", gsub("[-]", " x ", interactions[grepl("^G[-]", interactions)]))
    data$beta_columns <- paste0("Beta_", interactions[-1])
    data$se_columns   <- paste0("SE_Beta_", interactions[-1])
    data$robust_se_columns <- paste0("robust_", data$se_columns)
    data$int_colnames <- int_colnames
    
    covariances  <- unlist(lapply(combn(interactions[-1], 2, simplify = FALSE), FUN = function(x) {paste0(x[1], "_", x[2])}))
    cov_rownames <- unlist(lapply(combn(int_colnames[-1], 2, simplify = FALSE), FUN = function(x) {paste0(x[1], "_", x[2])}))
    robust_cov_rownames <- gsub("[_]", ", ", paste0("Cov<sub>R</sub>(", cov_rownames, ")"))
    cov_rownames <- gsub("[_]", ", ", paste0("Cov(", cov_rownames, ")"))
    data$covariances  <- paste0("Cov_Beta_", covariances)
    data$robust_covariances <- paste0("robust_", data$covariances)
    data$cov_rownames <- cov_rownames
    data$robust_cov_rownames <- robust_cov_rownames
    
    # Categorical interactions
    cat_interactions <- gsub("G[-]", "", interactions[-c(1,2)])
    cat_interactions <- colnames(df)[grepl(paste0("^N[_]", cat_interactions), colnames(df))]
    cat_interactions <- gsub("N[_]", "", cat_interactions)
    cat_n  <- paste0("<center>N<br>", gsub("[_]", " - ", cat_interactions), "</center>")
    cat_af <- paste0("<center>AF<br>", gsub("[_]", " - ", cat_interactions), "</center>")
    cat_colnames <- unlist(lapply(1:length(cat_n), FUN = function(x) c(cat_n[x], cat_af[x])))
    var_colnames <- c("-log<SUB>10</SUB>(p)", "SNP ID", "CHROM", "POS", "<center>NON-EFFECT<br>ALLELE</center>", "<center>EFFECT<br>ALLELE</center>", "<center>N<br>SAMPLES</center>", "AF", cat_colnames)
    
    cat_n  <- paste0("N_", cat_interactions)
    cat_af <- paste0("AF_", cat_interactions)
    cat_interactions <- unlist(lapply(1:length(cat_n), FUN = function(x) c(cat_n[x], cat_af[x])))
    data$cat_interactions <- cat_interactions
    data$var_colnames <- var_colnames
    
    # Manhattan plot data
    chrom_lengths <- chrom_lengths_hg38[extract_which_chr(df)]
    df <- add_cumulative_pos(df, chrom_lengths)
    df <- add_color(df, color1 = "black", color2 = "grey")
    
    x_breaks <- get_x_breaks(chrom_lengths_hg38)
    names(x_breaks)[20]="20"
    names(x_breaks)[22]="22"
    data$x_breaks <- x_breaks
    
    color_map <- unique(df$color)
    names(color_map) <- unique(df$color)
    data$color_map <- color_map
    
    data$df <- df
  })
  
  
  
  # UI - GWIS Panels -----------------------------------------------------------
  selectInputs <- reactive({
    list(input$gwis_choice, input$se_choice)
  })
  
  observeEvent(selectInputs(), {
    if (selectInputs()[[1]] == "marginal" & selectInputs()[[2]] == "modelbased") {
      show("gwis_mb_marginal_panel")
      hide("gwis_rb_marginal_panel")
      hide("gwis_mb_interaction_panel")
      hide("gwis_rb_interaction_panel")
      hide("gwis_mb_joint_panel")
      hide("gwis_rb_joint_panel")
    } else if (selectInputs()[[1]] == "marginal" & selectInputs()[[2]] == "robust") {
      hide("gwis_mb_marginal_panel")
      show("gwis_rb_marginal_panel")
      hide("gwis_mb_interaction_panel")
      hide("gwis_rb_interaction_panel")
      hide("gwis_mb_joint_panel")
      hide("gwis_rb_joint_panel")
    } else if (selectInputs()[[1]] == "interaction" & selectInputs()[[2]] == "modelbased") {
      hide("gwis_mb_marginal_panel")
      hide("gwis_rb_marginal_panel")
      show("gwis_mb_interaction_panel")
      hide("gwis_rb_interaction_panel")
      hide("gwis_mb_joint_panel")
      hide("gwis_rb_joint_panel")
    } else if (selectInputs()[[1]] == "interaction" & selectInputs()[[2]] == "robust") {
      hide("gwis_mb_marginal_panel")
      hide("gwis_rb_marginal_panel")
      hide("gwis_mb_interaction_panel")
      show("gwis_rb_interaction_panel")
      hide("gwis_mb_joint_panel")
      hide("gwis_rb_joint_panel")
    } else if (selectInputs()[[1]] == "joint" & selectInputs()[[2]] == "modelbased") {
      hide("gwis_mb_marginal_panel")
      hide("gwis_rb_marginal_panel")
      hide("gwis_mb_interaction_panel")
      hide("gwis_rb_interaction_panel")
      show("gwis_mb_joint_panel")
      hide("gwis_rb_joint_panel")
    }  else if (selectInputs()[[1]] == "joint" & selectInputs()[[2]] == "robust") {
      hide("gwis_mb_marginal_panel")
      hide("gwis_rb_marginal_panel")
      hide("gwis_mb_interaction_panel")
      hide("gwis_rb_interaction_panel")
      hide("gwis_mb_joint_panel")
      show("gwis_rb_joint_panel")
    }
  })
  
  
  
  # UI - Manhattan Box ---------------------------------------------------------
  output$mb_marginal_manhattan_box <- renderUI({
    manhattan_box("mb_marginal_manhattan_plot")
  })
  output$rb_marginal_manhattan_box <- renderUI({
    manhattan_box("rb_marginal_manhattan_plot")
  })
  
  output$mb_interaction_manhattan_box <- renderUI({
    manhattan_box("mb_interaction_manhattan_plot")
  })
  output$rb_interaction_manhattan_box <- renderUI({
    manhattan_box("rb_interaction_manhattan_plot")
  })
  
  output$mb_joint_manhattan_box <- renderUI({
    manhattan_box("mb_joint_manhattan_plot")
  })
  output$rb_joint_manhattan_box <- renderUI({
    manhattan_box("rb_joint_manhattan_plot")
  })
  
  

  # Manhattan Plots ------------------------------------------------------------
  output$mb_marginal_manhattan_plot <- renderPlot({
    plot_manhattan(data$df, data$x_breaks, data$color_map, "Marginal", FALSE)
  })
  output$rb_marginal_manhattan_plot <- renderPlot({
    plot_manhattan(data$df, data$x_breaks, data$color_map, "Marginal", TRUE)
  })
  
  output$mb_interaction_manhattan_plot <- renderPlot({
    plot_manhattan(data$df, data$x_breaks, data$color_map, "Interaction", FALSE)
  })
  output$rb_interaction_manhattan_plot <- renderPlot({
    plot_manhattan(data$df, data$x_breaks, data$color_map, "Interaction", TRUE)
  })
  
  output$mb_joint_manhattan_plot <- renderPlot({
    plot_manhattan(data$df, data$x_breaks, data$color_map, "Joint", FALSE)
  })
  output$rb_joint_manhattan_plot <- renderPlot({
    plot_manhattan(data$df, data$x_breaks, data$color_map, "Joint", TRUE)
  })
  
  
  
  # UI - QQ Box ---------------------------------------------------------------- 
  output$mb_marginal_qq_box <- renderUI({
    qq_box("mb_marginal_qq_plot")
  })
  output$rb_marginal_qq_box <- renderUI({
    qq_box("rb_marginal_qq_plot")
  })
  
  output$mb_interaction_qq_box <- renderUI({
    qq_box("mb_interaction_qq_plot")
  })
  output$rb_interaction_qq_box <- renderUI({
    qq_box("rb_interaction_qq_plot")
  })
  
  output$mb_joint_qq_box <- renderUI({
    qq_box("mb_joint_qq_plot")
  })
  output$rb_joint_qq_box <- renderUI({
    qq_box("rb_joint_qq_plot")
  })
  
  
  
  # QQ Plot --------------------------------------------------------------------
  output$mb_marginal_qq_plot <- renderPlot({
    plot_qq(data$df, "P_Value_Marginal")
  })
  output$rb_marginal_qq_plot <- renderPlot({
    plot_qq(data$df, "robust_P_Value_Marginal")
  })
  
  output$mb_interaction_qq_plot <- renderPlot({
    plot_qq(data$df, "P_Value_Interaction")
  })
  output$rb_interaction_qq_plot <- renderPlot({
    plot_qq(data$df, "robust_P_Value_Interaction")
  })
  
  output$mb_joint_qq_plot <- renderPlot({
    plot_qq(data$df, "P_Value_Joint")
  })
  output$rb_joint_qq_plot <- renderPlot({
    plot_qq(data$df, "robust_P_Value_Joint")
  })
  
  
  
  # UI - Variant Table Box -----------------------------------------------------
  observeEvent(input$mb_marginal_manhattan_plot_click, {
    output$mb_marginal_variants_table <- renderUI({
      variantTable_box("mb_marginal_manhattan_plot_table")
    })
    
    data$mb_marginal_nearest_points <- nearPoints(data$df, 
                                                     input$mb_marginal_manhattan_plot_click,
                                                     xvar = "cumulative_pos", 
                                                     yvar = "P_Value_Marginal")
    
    output$mb_marginal_manhattan_plot_table <- DT::renderDT({
      variant_table(data$mb_marginal_nearest_points, "P_Value_Marginal", data$var_colnames, data$cat_interactions)
    })
  })
  observeEvent(input$rb_marginal_manhattan_plot_click, {
    output$rb_marginal_variants_table <- renderUI({
      variantTable_box("rb_marginal_manhattan_plot_table")
    })
    
    data$rb_marginal_nearest_points <- nearPoints(data$df, 
                                                     input$rb_marginal_manhattan_plot_click,
                                                     xvar = "cumulative_pos", 
                                                     yvar = "robust_P_Value_Marginal")
    
    output$rb_marginal_manhattan_plot_table <- DT::renderDT({
      variant_table(data$rb_marginal_nearest_points, "robust_P_Value_Marginal", data$var_colnames, data$cat_interactions)
    })
  })
  
  observeEvent(input$mb_interaction_manhattan_plot_click, {
    output$mb_interaction_variants_table <- renderUI({
      variantTable_box("mb_interaction_manhattan_plot_table")
    })
    
    data$mb_interaction_nearest_points <- nearPoints(data$df, 
                                               input$mb_interaction_manhattan_plot_click,
                                               xvar = "cumulative_pos", 
                                               yvar = "P_Value_Interaction")
    
    output$mb_interaction_manhattan_plot_table <- DT::renderDT({
      variant_table(data$mb_interaction_nearest_points, "P_Value_Interaction", data$var_colnames, data$cat_interactions)
    })
  })
  observeEvent(input$rb_interaction_manhattan_plot_click, {
    output$rb_interaction_variants_table <- renderUI({
      variantTable_box("rb_interaction_manhattan_plot_table")
    })
    
    data$rb_interaction_nearest_points <- nearPoints(data$df, 
                                               input$rb_interaction_manhattan_plot_click,
                                               xvar = "cumulative_pos", 
                                               yvar = "robust_P_Value_Interaction")
    
    output$rb_interaction_manhattan_plot_table <- DT::renderDT({
      variant_table(data$rb_interaction_nearest_points, "robust_P_Value_Interaction", data$var_colnames, data$cat_interactions)
    })
  })
  
  observeEvent(input$mb_joint_manhattan_plot_click, {
    output$mb_joint_variants_table <- renderUI({
      variantTable_box("mb_joint_manhattan_plot_table")
    })

    data$mb_joint_nearest_points <- nearPoints(data$df, 
                                               input$mb_joint_manhattan_plot_click,
                                               xvar = "cumulative_pos", 
                                               yvar = "P_Value_Joint")
    
    output$mb_joint_manhattan_plot_table <- DT::renderDT({
      variant_table(data$mb_joint_nearest_points, "P_Value_Joint", data$var_colnames, data$cat_interactions)
    })
  })
  observeEvent(input$rb_joint_manhattan_plot_click, {
    output$rb_joint_variants_table <- renderUI({
      variantTable_box("rb_joint_manhattan_plot_table")
    })
    
    data$rb_joint_nearest_points <- nearPoints(data$df, 
                                               input$rb_joint_manhattan_plot_click,
                                               xvar = "cumulative_pos", 
                                               yvar = "robust_P_Value_Joint")
    
    output$rb_joint_manhattan_plot_table <- DT::renderDT({
      variant_table(data$rb_joint_nearest_points, "robust_P_Value_Joint", data$var_colnames, data$cat_interactions)
    })
  })
  
  
  
  # UI - Summary Statistics Tables ---------------------------------------------
  observeEvent(input$mb_marginal_manhattan_plot_table_rows_selected, {
    row <- input$mb_marginal_manhattan_plot_table_rows_selected
    
    output$mb_marginal_ss_tables <- renderUI({
      ssTable_box(data$mb_marginal_nearest_points$SNPID[row], "mb_marginal_ss_table")
    })
    
    ss_tables(output, "mb", "marginal", data$mb_marginal_nearest_points, row, data$int_colnames, data$beta_columns, data$se_columns, data$covariances, data$cov_rownames)
  })
  observeEvent(input$rb_marginal_manhattan_plot_table_rows_selected, {
    row <- input$rb_marginal_manhattan_plot_table_rows_selected
    
    output$rb_marginal_ss_tables <- renderUI({
      ssTable_box(data$rb_marginal_nearest_points$SNPID[row], "rb_marginal_ss_table")
    })
    
    ss_tables(output, "rb", "marginal", data$rb_marginal_nearest_points, row, data$int_colnames, data$beta_columns, data$robust_se_columns, data$robust_covariances, data$robust_cov_rownames)
  })
  
  observeEvent(input$mb_interaction_manhattan_plot_table_rows_selected, {
    row <- input$mb_interaction_manhattan_plot_table_rows_selected
    
    output$mb_interaction_ss_tables <- renderUI({
      ssTable_box(data$mb_interaction_nearest_points$SNPID[row], "mb_interaction_ss_table")
    })
    
    ss_tables(output, "mb", "interaction", data$mb_interaction_nearest_points, row, data$int_colnames, data$beta_columns, data$se_columns, data$covariances, data$cov_rownames)
  })
  observeEvent(input$rb_interaction_manhattan_plot_table_rows_selected, {
    row <- input$rb_interaction_manhattan_plot_table_rows_selected
    
    output$rb_interaction_ss_tables <- renderUI({
      ssTable_box(data$rb_interaction_nearest_points$SNPID[row], "rb_interaction_ss_table")
    })
    
    ss_tables(output, "rb", "interaction", data$rb_interaction_nearest_points, row, data$int_colnames, data$beta_columns, data$robust_se_columns, data$robust_covariances, data$robust_cov_rownames)
  })
  
  observeEvent(input$mb_joint_manhattan_plot_table_rows_selected, {
    row <- input$mb_joint_manhattan_plot_table_rows_selected
    
    output$mb_joint_ss_tables <- renderUI({
      ssTable_box(data$mb_joint_nearest_points$SNPID[row], "mb_joint_ss_table")
    })
    
    ss_tables(output, "mb", "joint", data$mb_joint_nearest_points, row, data$int_colnames, data$beta_columns, data$se_columns, data$covariances, data$cov_rownames)
  })
  observeEvent(input$rb_joint_manhattan_plot_table_rows_selected, {
    row <- input$rb_joint_manhattan_plot_table_rows_selected
    
    output$rb_joint_ss_tables <- renderUI({
      ssTable_box(data$rb_joint_nearest_points$SNPID[row], "rb_joint_ss_table")
    })
    
    ss_tables(output, "rb", "joint", data$rb_joint_nearest_points, row, data$int_colnames, data$beta_columns, data$robust_se_columns, data$robust_covariances, data$robust_cov_rownames)
  })
  

}


  