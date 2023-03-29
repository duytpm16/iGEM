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

  # UI - GENERAL --------------------------------------------------------------
  
  # use action buttons as tab selectors
    update_all <- function(x) {
    updateSelectInput(session, "tab",
                      choices = c("", "GWAS", "Diagnostics", "MARGINAL", "INTERACTION", "JOINT"),
                      label = "",
                      selected = x
    )
  }
  
  observeEvent(input$gwas, {
    update_all("GWAS")
  })
  observeEvent(input$gwas_marginal, {
    update_all("MARGINAL")
  })
  observeEvent(input$gwas_interaction, {
    update_all("INTERACTION")
  })
  observeEvent(input$gwas_joint, {
    update_all("JOINT")
  })
  observeEvent(input$diagnostics, {
    update_all("Diagnostics")
  })
  
  # update confirm button
  
  observeEvent(input$confirm, {
    updateButton(
      session, 
      inputId = "confirm", 
      label = "CONFIRM SELECTION", 
      icon = icon("bar-chart-o"), 
      style = "primary")
  })
  
  # hide the underlying selectInput in sidebar for better design
  observeEvent("", {
    hide("tab")
  })
  
  
  # DYNAMIC RENDER RULES ----------------------------------------------------
  
  observeEvent("", {
    show("gwas_panel")
    hide("diagnostics_panel")
    hide("gwas_mb_marginal_panel")
    hide("gwas_mb_interaction_panel")
    hide("gwas_mb_joint_panel")
  }, once = TRUE)
  
  observeEvent(input$gwas, {
    show("gwas_panel")
    hide("diagnostics_panel")
    hide("gwas_mb_marginal_panel")
    hide("gwas_mb_interaction_panel")
    hide("gwas_mb_joint_panel")
  })
  
  observeEvent(input$diagnostics, {
    show("diagnostics_panel")
    hide("patients_panel")
  })
  
  observeEvent(input$gwas_marginal, {
    show("gwas_panel")
    show("gwas_mb_marginal_panel")
    hide("gwas_mb_interaction_panel")
    hide("gwas_mb_joint_panel")
    hide("diagnostics_panel")
  })
  
  observeEvent(input$gwas_interaction, {
    show("gwas_panel")
    show("gwas_mb_interaction_panel")
    hide("gwas_mb_marginal_panel")
    hide("gwas_mb_joint_panel")
    hide("diagnostics_panel")
    
  })
  observeEvent(input$gwas_joint, {
    show("gwas_panel")
    show("gwas_mb_joint_panel")
    hide("gwas_mb_marginal_panel")
    hide("gwas_mb_interaction_panel")
    hide("diagnostics_panel")
  })
  
  
  # show active button with color
  
  observeEvent(input$tab, {
    x <- input$tab
    updateButton(session, "GWAS", style = {
      if (x == "GWAS") {
        paste("warning")
      } else {
        paste("success")
      }
    })
    
    updateButton(session, "diagnostics", style = {
      if (x == "Diagnostics") {
        paste("warning")
      } else {
        paste("success")
      }
    })
  })
  
  # READ INPUT FILE ------------------------------------------------------------
  data <- reactiveValues(df = NULL)
  observeEvent(input$inputFile, {
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
    data$int_colnames <- int_colnames
    
    covariances  <- unlist(lapply(combn(interactions[-1], 2, simplify = FALSE), FUN = function(x) {paste0(x[1], "_", x[2])}))
    cov_rownames <- unlist(lapply(combn(int_colnames[-1], 2, simplify = FALSE), FUN = function(x) {paste0(x[1], "_", x[2])}))
    cov_rownames <- gsub("[_]", ", ", paste0("cov(", cov_rownames, ")"))
    data$covariances  <- paste0("Cov_Beta_", covariances)
    data$cov_rownames <- cov_rownames
    
    # Categorical interactions
    cat_interactions <- gsub("G[-]", "", interactions[-c(1,2)])
    cat_interactions <- colnames(df)[grepl(paste0("^N_", cat_interactions), colnames(df))]
    cat_interactions <- gsub("N[_]", "", cat_interactions)
    cat_n  <- paste0("<center>N<br>", gsub("[_]", " - ", cat_interactions), "</center>")
    cat_af <- paste0("<center>AF<br>", gsub("[_]", " - ", cat_interactions), "</center>")
    cat_colnames <- unlist(lapply(1:length(cat_n), FUN = function(x) c(cat_n[x], cat_af[x])))
    var_colnames <- c("SNP ID", "CHROM", "POS", "<center>NON-EFFECT<br>ALLELE</center>", "<center>EFFECT<br>ALLELE</center>", "<center>N<br>SAMPLES</center>", "AF", cat_colnames)
    
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
  
  # UI - MB - MARGINAL - 1 ----------------------------------------------------------
  
  output$mb_marginal_manhattan_box <- renderUI({
      manhattan_box("mb_marginal_manhattan_plot")
  })
  
  observeEvent(input$gwas_marginal, {
    output$mb_marginal_manhattan_plot <- renderPlot({
      plot_manhattan(data$df, data$x_breaks, data$color_map, "Marginal")
    })
  })
  
  
  # UI - MB - MARGINAL - 2 -------------------------------------------------------

  output$mb_marginal_qq_box <- renderUI({
    box(
      title = p("Quantile-Quantile Plot", style = 'font-size:21px;'),
      status = "primary",
      collapsible = FALSE,
      solidHeader = FALSE,
      width = 12,
      withSpinner(
        plotOutput("mb_marginal_qq_plot", height = 300)
      )
    )
  })
  
  observeEvent(input$gwas_marginal, {
    output$mb_marginal_qq_plot <- renderPlot({
      ggplot(data$df, aes(x = P_Value_Marginal, y = P_Value_Marginal)) +
        
        geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
        geom_point() +
        theme(panel.background = element_blank(),
              panel.grid       = element_line(color = "grey97"),
              axis.line        = element_line(linewidth = 0.6)) +
        ylab(expression(paste('Observed ', -log[10](italic(p))))) +
        xlab(expression(paste('Expected ',-log[10](italic(p)))))
    })
  })
  

  # UI - MB - MARGINAL - 3 -------------------------------------------------------

  observeEvent(input$mb_marginal_manhattan_plot_click, {
    output$mb_marginal_variants_table <- renderUI({
        variant_box("mb_marginal_manhattan_plot_table")
    })
    
    data$mb_mrg_nearest_points <- nearPoints(data$df, input$mb_marginal_manhattan_plot_click,
                                             xvar = "cumulative_pos", yvar = "P_Value_Marginal")
    
    output$mb_marginal_manhattan_plot_table <- DT::renderDT({
        variant_table(data$mb_mrg_nearest_points, data$var_colnames, data$cat_interactions)
    })
  })
  
  
  # UI - MB - MARGINAL - 4 -------------------------------------------------------
  
  observeEvent(input$mb_marginal_manhattan_plot_table_rows_selected, {
    row <- input$mb_marginal_manhattan_plot_table_rows_selected
    
    output$mb_marginal_ss_tables <- renderUI({
      ss_box(data$mb_mrg_nearest_points$SNPID[row], "mb_marginal_manhattan_plot_table")
    })
    
    mb_ss_tables(output, "marginal", data$mb_mrg_nearest_points, row, data$beta_columns, data$se_columns, data$int_colnames, data$covariances, data$cov_rownames)
  })
  
  
  # UI - MB - INTERACTION - 1 ----------------------------------------------------------
  
  output$mb_interaction_manhattan_box <- renderUI({
    manhattan_box("mb_interaction_manhattan_plot")
  })
  
  observeEvent(input$gwas_interaction, {
    output$mb_interaction_manhattan_plot <- renderPlot({
      plot_manhattan(data$df, data$x_breaks, data$color_map, "Interaction")
    })
  })
  
  
  # UI - MB - INTERACTION - 2 -------------------------------------------------------
  
  output$mb_interaction_qq_box <- renderUI({
    box(
      title = p("Quantile-Quantile Plot", style = 'font-size:21px;'),
      status = "primary",
      collapsible = FALSE,
      solidHeader = FALSE,
      width = 12,
      withSpinner(
        plotOutput("mb_interaction_qq_plot", height = 300)
      )
    )
  })
  
  observeEvent(input$gwas_interaction, {
    output$mb_interaction_qq_plot <- renderPlot({
      ggplot(data$df, aes(x = P_Value_Interaction, y = P_Value_Interaction)) +
        
        geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
        geom_point() +
        theme(panel.background = element_blank(),
              panel.grid       = element_line(color = "grey97"),
              axis.line        = element_line(linewidth = 0.6)) +
        ylab(expression(paste('Observed ', -log[10](italic(p))))) +
        xlab(expression(paste('Expected ',-log[10](italic(p)))))
    })
  })
  
  
  # UI - MB - INTERACTION - 3 -------------------------------------------------------
  
  observeEvent(input$mb_interaction_manhattan_plot_click, {
    output$mb_interaction_variants_table <- renderUI({
        variant_box("mb_interaction_manhattan_plot_table")
    })
    
    data$mb_int_nearest_points <- nearPoints(data$df, input$mb_interaction_manhattan_plot_click,
                                             xvar = "cumulative_pos", yvar = "P_Value_Interaction")
    
    output$mb_interaction_manhattan_plot_table <- DT::renderDT({
        variant_table(data$mb_int_nearest_points, data$var_colnames, data$cat_interactions)
    })
  })
  
  
  # UI - MB - INTERACTION - 4 -------------------------------------------------------
  
  observeEvent(input$mb_interaction_manhattan_plot_table_rows_selected, {
    row <- input$mb_interaction_manhattan_plot_table_rows_selected
    
    output$mb_interaction_ss_tables <- renderUI({
      ss_box(data$mb_int_nearest_points$SNPID[row], "mb_interaction_manhattan_plot_table")
    })
    
    mb_ss_tables(output, "interaction", data$mb_int_nearest_points, row, data$beta_columns, data$se_columns, data$int_colnames, data$covariances, data$cov_rownames)
  })
  
  
  # UI - MB - JOINT - 1 ----------------------------------------------------------
  
  output$mb_joint_manhattan_box <- renderUI({
    manhattan_box("mb_joint_manhattan_plot")
  })
  
  observeEvent(input$gwas_joint, {
    output$mb_joint_manhattan_plot <- renderPlot({
      plot_manhattan(data$df, data$x_breaks, data$color_map, "Joint")
    })
  })
  
  
  # UI - MB - JOINT - 2 -------------------------------------------------------
  
  output$mb_joint_qq_box <- renderUI({
    box(
      title = p("Quantile-Quantile Plot", style = 'font-size:21px;'),
      status = "primary",
      collapsible = FALSE,
      solidHeader = FALSE,
      width = 12,
      withSpinner(
        plotOutput("mb_joint_qq_plot", height = 300)
      )
    )
  })
  
  observeEvent(input$gwas_joint, {
    output$mb_joint_qq_plot <- renderPlot({
      ggplot(data$df, aes(x = P_Value_Joint, y = P_Value_Joint)) +
        
        geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
        geom_point() +
        theme(panel.background = element_blank(),
              panel.grid       = element_line(color = "grey97"),
              axis.line        = element_line(linewidth = 0.6)) +
        ylab(expression(paste('Observed ', -log[10](italic(p))))) +
        xlab(expression(paste('Expected ',-log[10](italic(p)))))
    })
  })
  
  
  # UI - MB - JOINT - 3 -------------------------------------------------------
  
  observeEvent(input$mb_joint_manhattan_plot_click, {
    output$mb_joint_variants_table <- renderUI({
      variant_box("mb_joint_manhattan_plot_table")
    })
    
    data$mb_jnt_nearest_points <- nearPoints(data$df, input$mb_joint_manhattan_plot_click,
                                             xvar = "cumulative_pos", yvar = "P_Value_Joint")
    
    output$mb_joint_manhattan_plot_table <- DT::renderDT({
      variant_table(data$mb_jnt_nearest_points, data$var_colnames, data$cat_interactions)
    })
  })
  
  
  # UI - MB - JOINT - 4 -------------------------------------------------------
  
  observeEvent(input$mb_joint_manhattan_plot_table_rows_selected, {
    row <- input$mb_joint_manhattan_plot_table_rows_selected
    
    output$mb_joint_ss_tables <- renderUI({
      ss_box(data$mb_jnt_nearest_points$SNPID[row], "mb_joint_manhattan_plot_table")
    })
    
    mb_ss_tables(output, "joint", data$mb_jnt_nearest_points, row, data$beta_columns, data$se_columns, data$int_colnames, data$covariances, data$cov_rownames)
  })
}
  

  
