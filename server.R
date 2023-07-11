server <- function(input, output, session) {

  # UI - GENERAL ---------------------------------------------------------------
  data <- reactiveValues(df = NULL)
  
  mh_data <- reactiveValues(mb_marginal    = NULL,
                            rb_marginal    = NULL,
                            mb_interaction = NULL,
                            rb_interaction = NULL,
                            mb_joint       = NULL,
                            rb_joint       = NULL)
  
  mh_nvar <- reactiveValues(mb_marginal_var_count    = NULL,
                            rb_marginal_var_count    = NULL,
                            mb_interaction_var_count = NULL,
                            rb_interaction_var_count = NULL,
                            mb_joint_var_count       = NULL,
                            rb_joint_var_count       = NULL)
  
  qq_data <- reactiveValues(mb_marginal    = NULL,
                            rb_marginal    = NULL,
                            mb_interaction = NULL,
                            rb_interaction = NULL,
                            mb_joint       = NULL,
                            rb_joint       = NULL)
  
  columnExists <- reactiveValues(p_value_marginal           = FALSE,
                                 p_value_interaction        = FALSE,
                                 p_value_joint              = FALSE,
                                 robust_p_value_marginal    = FALSE,
                                 robust_p_value_interaction = FALSE,
                                 robust_p_value_joint       = FALSE)
  
  # READ INPUT FILE ------------------------------------------------------------
  observeEvent(input$inputFile, {
    req(input$inputFile)
    
    ## Read File
    df   <- fread(input$inputFile$datapath, sep = "\t", header = T)
    nvar <- nrow(df)
    coln <- colnames(df)
    gc(verbose = FALSE)
    
    setorderv(df, c("CHR", "POS"))
    chrom_lengths     <- chrom_lengths_hg38[df[,kit::funique(CHR)]]
    cum_chrom_lengths <- get_cumulative_length(chrom_lengths)
    df[, cumulative_pos := POS + cum_chrom_lengths[CHR]]
    df[, index := 1:nvar]
    gc(verbose = FALSE)
  
    ## P-values------------------------------------------------------------------
    index <- c()
    # P-values
    if ("P_Value_Marginal" %in% coln) {
      columnExists$p_value_marginal = TRUE
      
      qq_data$mb_marginal_lambda <- format(round(median(qchisq(1-df$P_Value_Marginal,1)) / qchisq(0.5, 1), 2), nsmall = 2)
      
      df[, P_Value_Marginal := -log10(P_Value_Marginal)]
      
      qq_data$mb_marginal <- fastqq::drop_dense(df$P_Value_Marginal, -log10(stats::ppoints(nvar)))
      
      subDF <- subset_data(df[,c("CHR", "cumulative_pos", "P_Value_Marginal")], "P_Value_Marginal", nvar)
  
      index <- c(index, subDF$index)
      mh_data$mb_marginal <- subDF[!subDF$duplicated, c("index", "CHR", "POS", "LOGP")]
      mh_nvar$mb_marginal_var_count <- dplyr::count(mh_data$mb_marginal, CHR)
      gc(verbose = FALSE)
    }
    
    if ("robust_P_Value_Marginal" %in% coln) {
      columnExists$robust_p_value_marginal = TRUE
      
      qq_data$rb_marginal_lambda <- format(round(median(qchisq(1-df$robust_P_Value_Marginal,1)) / qchisq(0.5, 1), 2), nsmall = 2)
      
      df[, robust_P_Value_Marginal := -log10(robust_P_Value_Marginal)]
      
      qq_data$rb_marginal <- fastqq::drop_dense(df$robust_P_Value_Marginal, -log10(stats::ppoints(nvar)))
      
      subDF <- subset_data(df[,c("CHR", "cumulative_pos", "robust_P_Value_Marginal")], "robust_P_Value_Marginal", nvar)
      
      index <- c(index, subDF$index)
      mh_data$rb_marginal <- subDF[!subDF$duplicated, c("index", "CHR", "POS", "LOGP")]
      mh_nvar$rb_marginal_var_count <- dplyr::count(mh_data$rb_marginal, CHR)
      gc(verbose = FALSE)
    }
    
    if ("P_Value_Interaction" %in% coln) {
      columnExists$p_value_interaction = TRUE
      
      qq_data$mb_interaction_lambda <- format(round(median(qchisq(1-df$P_Value_Interaction, 1)) / qchisq(0.5, 1), 2), nsmall = 2)
      
      df[, P_Value_Interaction := -log10(P_Value_Interaction)]
      
      qq_data$mb_interaction <- fastqq::drop_dense(df$P_Value_Interaction, -log10(stats::ppoints(nvar)))
      
      subDF <- subset_data(df[,c("CHR", "cumulative_pos", "P_Value_Interaction")], "P_Value_Interaction", nvar)
      
      index <- c(index, subDF$index)
      mh_data$mb_interaction <- subDF[!subDF$duplicated, c("index", "CHR", "POS", "LOGP") ]
      mh_nvar$mb_interaction_var_count <- dplyr::count(mh_data$mb_interaction, CHR)
      gc(verbose = FALSE)
    }
    
    if ("robust_P_Value_Interaction" %in% coln) {
      columnExists$robust_p_value_interaction = TRUE
      
      qq_data$rb_interaction_lambda <- format(round(median(qchisq(1-df$robust_P_Value_Interaction, 1)) / qchisq(0.5, 1), 2), nsmall = 2)
      
      df[, robust_P_Value_Interaction := -log10(robust_P_Value_Interaction)]
      
      qq_data$rb_interaction <- fastqq::drop_dense(df$robust_P_Value_Interaction, -log10(stats::ppoints(nvar)))
      
      subDF <- subset_data(df[,c("CHR", "cumulative_pos", "robust_P_Value_Interaction")], "robust_P_Value_Interaction", nvar)
      
      index <- c(index, subDF$index)
      mh_data$rb_interaction <- subDF[!subDF$duplicated, c("index", "CHR", "POS", "LOGP") ]
      mh_nvar$rb_interaction_var_count <- dplyr::count(mh_data$rb_interaction, CHR)
      gc(verbose = FALSE)
    }
    
    if ("P_Value_Joint" %in% coln) {
      columnExists$p_value_joint = TRUE
      
      qq_data$mb_joint_lambda <- format(round(median(qchisq(1-df$P_Value_Joint, 1)) / qchisq(0.5, 1), 2), nsmall = 2)
      
      df[, P_Value_Joint := -log10(P_Value_Joint)]
      
      qq_data$mb_joint <- fastqq::drop_dense(df$P_Value_Joint, -log10(stats::ppoints(nvar)))
      
      subDF <- subset_data(df[,c("CHR", "cumulative_pos", "P_Value_Joint")], "P_Value_Joint", nvar)
      
      index <- c(index, subDF$index)
      mh_data$mb_joint <- subDF[!subDF$duplicated, c("index", "CHR", "POS", "LOGP") ]
      mh_nvar$mb_joint_var_count <- dplyr::count(mh_data$mb_joint, CHR)
      gc(verbose = FALSE)
    }
    
    if ("robust_P_Value_Joint" %in% coln) {
      columnExists$robust_p_value_joint = TRUE
      
      qq_data$rb_joint_lambda <- format(round(median(qchisq(1-df$robust_P_Value_Joint, 1)) / qchisq(0.5, 1), 2), nsmall = 2)
      
      df[, robust_P_Value_Joint := -log10(robust_P_Value_Joint)]
      
      qq_data$rb_joint <- fastqq::drop_dense(df$robust_P_Value_Joint, -log10(stats::ppoints(nvar)))
      
      subDF <- subset_data(df[,c("CHR", "cumulative_pos", "robust_P_Value_Joint")], "robust_P_Value_Joint", nvar)
      
      index <- c(index, subDF$index)
      mh_data$rb_joint <- subDF[!subDF$duplicated, c("index", "CHR", "POS", "LOGP") ]
      mh_nvar$rb_joint_var_count <- dplyr::count(mh_data$rb_joint, CHR)
      gc(verbose = FALSE)
    }
    
    index <- index[!fduplicated(index)]
    index <- index[order(index)]
    
    df <- df[index, ]
    df$new_index <- 1:nrow(df)
    for (x in names(mh_data)) {
      mh_data[[x]]$index <- df$new_index[match(mh_data[[x]]$index, df$index)]
    }
    gc(verbose = FALSE)
    
    
    # Get interactions including marginal columns
    interactions <- c("Beta_Marginal", "Beta_G", coln[grepl("^Beta_G-", coln)])
    interactions <- gsub("Beta_", "", interactions)
    int_colnames <- c("Marginal", "Main", gsub("G[-]", "G x ", interactions[grepl("^G[-]", interactions)]))
    
    ## Betas and SE-------------------------------------------------------------
    data$mb_beta <- paste0("Beta_", interactions[-1])
    data$rb_beta <- paste0("Beta_", interactions[-1])
    data$mb_se   <- paste0("SE_Beta_", interactions[-1])
    data$rb_se   <- paste0("robust_SE_Beta_", interactions[-1])
    data$int_colnames <- int_colnames[-1]
    
    ## Covariance----------------------------------------------------------------
    covariances  <- unlist(lapply(combn(interactions[-1], 2, simplify = FALSE), FUN = function(x) {paste0(x[1], "_", x[2])}))
    cov_rownames <- unlist(lapply(combn(int_colnames[-1], 2, simplify = FALSE), FUN = function(x) {paste0(x[1], ", ", x[2])}))
    data$mb_covs <- paste0("Cov_Beta_", covariances)
    data$rb_covs <- paste0("robust_Cov_Beta_", covariances)
    data$mb_cov_rownames <- paste0("Cov(", cov_rownames, ")")
    data$rb_cov_rownames <- paste0("Cov<sub>R</sub>(", cov_rownames, ")")
    
    ## Categorical---------------------------------------------------------------
    cat_interactions <- gsub("G[-]", "", interactions[-c(1,2)])
    cat_interactions <- coln[grepl(paste0("^N[_]", cat_interactions, collapse = "|"), coln)]
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
    data$x_breaks <- get_x_breaks(chrom_lengths_hg38)
  
    data$df <- as.data.frame(df)
  })
  
  
  
  # Inputs ---------------------------------------------------------------------
  selectInputs <- reactive({
    list(input$gwis_choice, input$se_choice)
  })
  
  mh_sigThreshold <- reactive({
    input$mh_sigThreshold
  })
  
  mh_sigColor <- reactive({
    input$mh_sigColor
  })
  
  
  
  # GWIS Panels ----------------------------------------------------------------
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
  
  
  
  # Manhattan Box --------------------------------------------------------------
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
  
  
  
  # Manhattan Colors------------------------------------------------------------
  mh_chrColor <- reactive({
    dfs <- list()
    colors <- strsplit(input$mh_chrColor, split = ";")[[1]]
    
    if (!is.null(mh_nvar$mb_marginal_var_count)) {
      dfs$mb_marginal <- get_chr_colors(mh_nvar$mb_marginal_var_count, colors)
    }
    
    if (!is.null(mh_nvar$rb_marginal_var_count)) {
      dfs$rb_marginal <- get_chr_colors(mh_nvar$rb_marginal_var_count, colors)
    }
    
    if (!is.null(mh_nvar$mb_interaction_var_count)) {
      dfs$mb_interaction <- get_chr_colors(mh_nvar$mb_interaction_var_count, colors)
    }
    
    if (!is.null(mh_nvar$rb_interaction_var_count)) {
      dfs$rb_interaction <- get_chr_colors(mh_nvar$rb_interaction_var_count, colors)
    }
    
    if (!is.null(mh_nvar$mb_joint_var_count)) {
      dfs$mb_joint <- get_chr_colors(mh_nvar$mb_joint_var_count, colors)
    }
    
    if (!is.null(mh_nvar$rb_joint_var_count)) {
      dfs$rb_joint <- get_chr_colors(mh_nvar$rb_joint_var_count, colors)
    }
    
    return(dfs)
  })
  

  
  # Manhattan Plot -------------------------------------------------------------
  output$mb_marginal_manhattan_plot <- renderPlot({
    manhattan_plot(mh_data$mb_marginal, data$x_breaks, mh_sigThreshold(), mh_sigColor(), mh_chrColor()$mb_marginal)
  })
  
  output$rb_marginal_manhattan_plot <- renderPlot({
    manhattan_plot(mh_data$rb_marginal, data$x_breaks, mh_sigThreshold(), mh_sigColor(), mh_chrColor()$rb_marginal)
  })
  
  output$mb_interaction_manhattan_plot <- renderPlot({
    manhattan_plot(mh_data$mb_interaction, data$x_breaks, mh_sigThreshold(), mh_sigColor(), mh_chrColor()$mb_interaction)
  })
  
  output$rb_interaction_manhattan_plot <- renderPlot({
    manhattan_plot(mh_data$rb_interaction, data$x_breaks, mh_sigThreshold(), mh_sigColor(), mh_chrColor()$rb_interaction)
  })
  
  output$mb_joint_manhattan_plot <- renderPlot({
    manhattan_plot(mh_data$mb_joint, data$x_breaks, mh_sigThreshold(), mh_sigColor(), mh_chrColor()$mb_joint)
  })

  output$rb_joint_manhattan_plot <- renderPlot({
    manhattan_plot(mh_data$rb_joint, data$x_breaks, mh_sigThreshold(), mh_sigColor(), mh_chrColor()$rb_joint)
  })
  
  
  
  # Manhattan Tooltip-----------------------------------------------------------
  output$mb_marginal_manhattan_plot_hover_info <- renderUI({
    manhattan_tooltip(input[["mb_marginal_manhattan_plot_hover"]], mh_data$mb_marginal)
  })

  output$rb_marginal_manhattan_plot_hover_info <- renderUI({
    manhattan_tooltip(input[["rb_marginal_manhattan_plot_hover"]], mh_data$rb_marginal)
  })

  output$mb_interaction_manhattan_plot_hover_info <- renderUI({
    manhattan_tooltip(input[["mb_interaction_manhattan_plot_hover"]], mh_data$mb_interaction)
  })

  output$rb_interaction_manhattan_plot_hover_info <- renderUI({
    manhattan_tooltip(input[["rb_interaction_manhattan_plot_hover"]], mh_data$rb_interaction)
  })

  output$mb_joint_manhattan_plot_hover_info <- renderUI({
      manhattan_tooltip(input[["mb_joint_manhattan_plot_hover"]], mh_data$mb_joint)
  })
  
  output$rb_joint_manhattan_plot_hover_info <- renderUI({
    manhattan_tooltip(input[["rb_joint_manhattan_plot_hover"]], mh_data$rb_joint)
  })
  
  
  
  # QQ Box ---------------------------------------------------------------------
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
    qq_plot(qq_data$mb_marginal, qq_data$mb_marginal_lambda)
  })
  
  output$rb_marginal_qq_plot <- renderPlot({
    qq_plot(qq_data$rb_marginal, qq_data$rb_marginal_lambda)
  })
  
  output$mb_interaction_qq_plot <- renderPlot({
    qq_plot(qq_data$mb_interaction, qq_data$mb_interaction_lambda)
  })
  
  output$rb_interaction_qq_plot <- renderPlot({
    qq_plot(qq_data$rb_interaction, qq_data$rb_interaction_lambda)
  })
  
  output$mb_joint_qq_plot <- renderPlot({
    qq_plot(qq_data$mb_joint, qq_data$mb_joint_lambda)
  })
  
  output$rb_joint_qq_plot <- renderPlot({
    qq_plot(qq_data$rb_joint, qq_data$rb_joint_lambda)
  })
  
  
  
  # Variant Table Box ----------------------------------------------------------
  output$mb_marginal_variants_table <- renderUI({
    variantTable_box("mb_marginal_manhattan_plot_table")
  })
  
  output$rb_marginal_variants_table <- renderUI({
    variantTable_box("rb_marginal_manhattan_plot_table")
  })
  
  output$mb_interaction_variants_table <- renderUI({
    variantTable_box("mb_interaction_manhattan_plot_table")
  })
  
  output$rb_interaction_variants_table <- renderUI({
    variantTable_box("rb_interaction_manhattan_plot_table")
  })
  
  output$mb_joint_variants_table <- renderUI({
    variantTable_box("mb_joint_manhattan_plot_table")
  })
  
  output$rb_joint_variants_table <- renderUI({
    variantTable_box("rb_joint_manhattan_plot_table")
  })
  
  
  
  # Variant Table Nearest Points------------------------------------------------
  observeEvent(input$mb_marginal_manhattan_plot_click, {
    data$mb_marginal_nearest_points <- nearPoints(mh_data$mb_marginal, input$mb_marginal_manhattan_plot_click, xvar = "POS", yvar = "LOGP")
  })
  
  observeEvent(input$rb_marginal_manhattan_plot_click, {
    data$rb_marginal_nearest_points <- nearPoints(mh_data$rb_marginal, input$rb_marginal_manhattan_plot_click, xvar = "POS", yvar = "LOGP")
  })
  
  observeEvent(input$mb_interaction_manhattan_plot_click, {
    data$mb_interaction_nearest_points <- nearPoints(mh_data$mb_interaction, input$mb_interaction_manhattan_plot_click, xvar = "POS", yvar = "LOGP")
  })
  
  observeEvent(input$rb_interaction_manhattan_plot_click, {
    data$rb_interaction_nearest_points <- nearPoints(mh_data$rb_interaction, input$rb_interaction_manhattan_plot_click, xvar = "POS", yvar = "LOGP")
  })
  
  observeEvent(input$mb_joint_manhattan_plot_click, {
    data$mb_joint_nearest_points <- nearPoints(mh_data$mb_joint, input$mb_joint_manhattan_plot_click, xvar = "POS", yvar = "LOGP")
  })
  
  observeEvent(input$rb_joint_manhattan_plot_click, {
    data$rb_joint_nearest_points <- nearPoints(mh_data$rb_joint, input$rb_joint_manhattan_plot_click, xvar = "POS", yvar = "LOGP")
  })
  
  
  
  # Variant Table---------------------------------------------------------------
  output$mb_marginal_manhattan_plot_table <- DT::renderDT({
    req(data$df[data$mb_marginal_nearest_points$index, ])
    variantTable(data$df[data$mb_marginal_nearest_points$index, ], "P_Value_Marginal", data$var_colnames, data$cat_interactions)
  })
  
  output$rb_marginal_manhattan_plot_table <- DT::renderDT({
    req(data$df[data$rb_marginal_nearest_points$index, ])
    variantTable(data$df[data$rb_marginal_nearest_points$index, ], "robust_P_Value_Marginal", data$var_colnames, data$cat_interactions)
  })
  
  output$mb_interaction_manhattan_plot_table <- DT::renderDT({
    req(data$df[data$mb_interaction_nearest_points$index, ])
    variantTable(data$df[data$mb_interaction_nearest_points$index, ], "P_Value_Interaction", data$var_colnames, data$cat_interactions)
  })
  
  output$rb_interaction_manhattan_plot_table <- DT::renderDT({
    req(data$df[data$rb_interaction_nearest_points$index, ])
    variantTable(data$df[data$rb_interaction_nearest_points$index, ], "robust_P_Value_Interaction", data$var_colnames, data$cat_interactions)
  })
  
  output$mb_joint_manhattan_plot_table <- DT::renderDT({
    req(data$df[data$mb_joint_nearest_points$index, ])
    variantTable(data$df[data$mb_joint_nearest_points$index, ], "P_Value_Joint", data$var_colnames, data$cat_interactions)
  })
  
  output$rb_joint_manhattan_plot_table <- DT::renderDT({
    req(data$df[data$rb_joint_nearest_points$index, ])
    variantTable(data$df[data$rb_joint_nearest_points$index, ], "robust_P_Value_Joint", data$var_colnames, data$cat_interactions)
  })
  
  
  
  # Summary Statistics Box------------------------------------------------------
  output$mb_marginal_ssTables <- renderUI({
    ssTable_box("mb_marginal_ssTable")
  })
  
  output$rb_marginal_ssTables <- renderUI({
    ssTable_box("rb_marginal_ssTable")
  })
  
  output$mb_interaction_ssTables <- renderUI({
    ssTable_box("mb_interaction_ssTable")
  })
  
  output$rb_interaction_ssTables <- renderUI({
    ssTable_box("rb_interaction_ssTable")
  })
  
  output$mb_joint_ssTables <- renderUI({
    ssTable_box("mb_joint_ssTable")
  })
  
  output$rb_joint_ssTables <- renderUI({
    ssTable_box( "rb_joint_ssTable")
  })
  
  
  
  # Variant Table Row Selected -------------------------------------------------
  observeEvent(ignoreInit = TRUE, list(data$mb_marginal_nearest_points, input$mb_marginal_manhattan_plot_table_rows_selected), {
    row <- input$mb_marginal_manhattan_plot_table_rows_selected
    ssTables(output, "mb", "marginal", data$df[data$mb_marginal_nearest_points$index[row], ], data$int_colnames, data$mb_beta, data$mb_se, data$mb_covs, data$mb_cov_rownames)
  })

  observeEvent(ignoreInit = TRUE, list(data$rb_marginal_nearest_points, input$rb_marginal_manhattan_plot_table_rows_selected), {
    row <- input$rb_marginal_manhattan_plot_table_rows_selected
    ssTables(output, "rb", "marginal", data$df[data$rb_marginal_nearest_points$index[row], ], data$int_colnames, data$rb_beta, data$rb_se, data$rb_covs, data$rb_cov_rownames)
  })

  observeEvent(ignoreInit = TRUE, list(data$mb_interaction_nearest_points, input$mb_interaction_manhattan_plot_table_rows_selected), {
    row <- input$mb_interaction_manhattan_plot_table_rows_selected
    ssTables(output, "mb", "interaction", data$df[data$mb_interaction_nearest_points$index[row], ], data$int_colnames, data$mb_beta, data$mb_se, data$mb_covs, data$mb_cov_rownames)
  })

  observeEvent(ignoreInit = TRUE, list(data$rb_interaction_nearest_points, input$rb_interaction_manhattan_plot_table_rows_selected), {
    row <- input$rb_interaction_manhattan_plot_table_rows_selected
    ssTables(output, "rb", "interaction", data$df[data$rb_interaction_nearest_points$index[row], ], data$int_colnames, data$rb_beta, data$rb_se, data$rb_covs, data$rb_cov_rownames)
  })

  observeEvent(ignoreInit = TRUE, list(data$mb_joint_nearest_points, input$mb_joint_manhattan_plot_table_rows_selected), {
    row <- input$mb_joint_manhattan_plot_table_rows_selected
    ssTables(output, "mb", "joint", data$df[data$mb_joint_nearest_points$index[row], ], data$int_colnames, data$mb_beta, data$mb_se, data$mb_covs, data$mb_cov_rownames)
  })

  observeEvent(ignoreInit = TRUE, list(data$rb_joint_nearest_points, input$rb_joint_manhattan_plot_table_rows_selected), {
    row <- input$rb_joint_manhattan_plot_table_rows_selected
    ssTables(output, "rb", "joint", data$df[data$rb_joint_nearest_points$index[row], ], data$int_colnames, data$rb_beta, data$rb_se, data$rb_covs, data$rb_cov_rownames)
  })
}