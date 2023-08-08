server <- function(input, output, session) {

  # UI - GENERAL ---------------------------------------------------------------
  data <- reactiveValues(df = NULL)
  
  session$onSessionEnded(function() {
    stopApp()
  })
  
  # READ INPUT FILE ------------------------------------------------------------
  observeEvent(input$inputFile, {
    req(input$inputFile)
    
    ## Read File
    df   <- fread(input$inputFile$datapath, sep = "\t", header = T)
    nvar <- nrow(df)
    coln <- colnames(df)
    
    setorderv(df, c("CHR", "POS"))
    chrom_lengths     <- chrom_lengths_hg38[df[, kit::funique(CHR)]]
    cum_chrom_lengths <- get_cumulative_length(chrom_lengths)
    df[, cumulative_pos := POS + cum_chrom_lengths[CHR]]
  
    ## P-values------------------------------------------------------------------
    pvalue_columns <- c("P_Value_Marginal",    "robust_P_Value_Marginal",
                        "P_Value_Interaction", "robust_P_Value_Interaction",
                        "P_Value_Joint",       "robust_P_Value_Joint")
    pvalue_columns <- pvalue_columns[pvalue_columns %in% coln]
    
    
    data$lambda <- df[, lapply(.SD, function(x) {gc(); formatC(median(qchisq(x, 1, lower.tail = FALSE)) / qchisq(0.5, 1), digits = 2, format = "f")}), .SDcols=pvalue_columns]
    gc()
    df[, (pvalue_columns) := lapply(.SD, function(x) -log10(x)), .SDcols = pvalue_columns]
    data$qq <- sapply(pvalue_columns,  function(x) {gc(); fastqq::drop_dense(df[[x]], -log10(stats::ppoints(nvar)))}, simplify = FALSE, USE.NAMES = TRUE)
    gc()
    
    mh_keep <- list()
    keep <- vector(mode = "logical", length = nvar)
    for (pcol in pvalue_columns) {
      round_df <- data.table(index = 1:nvar,
                             round_cpos = df[["cumulative_pos"]],
                             round_pval = round(df[[pcol]], digits = 1))
      if (nvar > 500000) {gc()}
      round_df[round_pval >   8, round_cpos := plyr::round_any(round_cpos, accuracy = 100000)]
      round_df[round_pval <=  8, ':='(round_pval = plyr::round_any(round_pval, accuracy = 0.2), round_cpos = plyr::round_any(round_cpos, accuracy = 500000))]
      if (nvar > 500000) {gc()}
      round_df <- round_df[!fduplicated(round_df[,c("round_cpos", "round_pval")]), ]
      if (nvar > 500000) {gc()}
      keep[df[[pcol]] > 8] <- TRUE
      keep[round_df[["index"]]] <- TRUE
      
      mh_keep[[pcol]] <- round_df[["index"]]
    }
    rm(round_df)

    df <- df[keep, ]
    gc()
    
    index   <- c(1:nvar)[keep]
    index2  <- c(1:length(index))
    mh_keep <- lapply(mh_keep, FUN = function(x) index2[match(x, index)])
    data$mh_nvar <- lapply(mh_keep, FUN = function(x) as.data.frame(df[x, .N, by = CHR]))
    data$mh_data <- mh_keep
    rm(index, index2, keep, mh_keep)
    gc(verbose = FALSE)
    
    
    # Get interactions including marginal columns
    interactions <- c("Beta_Marginal", "Beta_G", coln[grepl("^Beta_G-", coln)])
    interactions <- gsub("Beta_", "", interactions)
    int_colnames <- c("Marginal", "Main", gsub("G[-]", "G x ", interactions[-c(1:2)]))

    
    ## Betas--------------------------------------------------------------------
    beta_columns <- paste0("Beta_", interactions)
    data$mb_beta <- beta_columns
    data$rb_beta <- beta_columns
    data$mb_betaG_prefix   <- "Beta_G"
    data$rb_betaG_prefix   <- "Beta_G"
    data$mb_betaInt_prefix <- "Beta_G-"
    data$rb_betaInt_prefix <- "Beta_G-"
    if(any(grepl("^robust_Beta_", coln))){
      data$rb_beta <- paste0("robust_", beta_columns)
      data$rb_betaG_prefix   <- "robust_Beta_G"
      data$rb_betaInt_prefix <- "robust_Beta_G-"
    }
    
    ## SE-----------------------------------------------------------------------
    data$mb_se <- paste0("SE_Beta_", interactions)
    data$rb_se <- paste0("robust_SE_Beta_", interactions)
    data$int_colnames <- int_colnames
    
    
    ## Covariance----------------------------------------------------------------
    covariances  <- unlist(lapply(combn(interactions[-1], 2, simplify = FALSE), FUN = function(x) {paste0(x[1], "_", x[2])}))
    cov_rownames <- unlist(lapply(combn(int_colnames[-1], 2, simplify = FALSE), FUN = function(x) {paste0(x[1], ", ", x[2])}))
    data$mb_covs <- paste0("Cov_Beta_", covariances)
    data$rb_covs <- paste0("robust_Cov_Beta_", covariances)
    data$mb_cov_rownames <- paste0("Cov(", cov_rownames, ")")
    data$rb_cov_rownames <- paste0("Cov<sub>R</sub>(", cov_rownames, ")")
    
    
    ## Variant Info-------------------------------------------------------------
    var_columns  <- c("SNPID", "CHR", "POS", "Non_Effect_Allele", "Effect_Allele", "N_Samples", "AF")
    var_colnames <- c("-log<SUB>10</SUB>(p)", "SNP ID", "CHROM", "POS", "<center>NON-EFFECT<br>ALLELE</center>", "<center>EFFECT<br>ALLELE</center>", "<center>N<br>SAMPLES</center>", "AF")
    
    
    ## Categorical--------------------------------------------------------------
    interactions <- gsub("G[-]", "", interactions[-c(1,2)])
    cat_interactions  <- coln[grepl(paste0("^N[_]", interactions, collapse = "|"), coln)]
  
    mb_mxi_dfs <- vector(mode = "list", length = length(interactions))
    rb_mxi_dfs <- vector(mode = "list", length = length(interactions))
    names(mb_mxi_dfs) <- interactions
    names(rb_mxi_dfs) <- interactions
    categorical_ints  <- vector(mode = "character", length = 0)
    if (length(cat_interactions) != 0) {
      cat_interactions <- gsub("N[_]", "", cat_interactions)
      
      for (cat in cat_interactions) {
        cat_split <- strsplit(cat, "_")[[1]]
        stratum   <- as.numeric(cat_split[length(cat_split)])
        intname   <- paste0(cat_split[-length(cat_split)], collapse = "_")
        mb_mxi_dfs[[intname]] <- c(mb_mxi_dfs[[intname]], stratum)
        rb_mxi_dfs[[intname]] <- c(rb_mxi_dfs[[intname]], stratum)
      }
      
      for (x in names(mb_mxi_dfs)) {
        if (!is.null(mb_mxi_dfs[[x]])) {
          mb_mxi_dfs[[x]] <- data.frame(i = rep(x, length(mb_mxi_dfs[[x]])), e = as.factor(mb_mxi_dfs[[x]][order(mb_mxi_dfs[[x]])]), b = rep(paste0(data$mb_betaInt_prefix, x), length(mb_mxi_dfs[[x]])))
          rb_mxi_dfs[[x]] <- data.frame(i = rep(x, length(rb_mxi_dfs[[x]])), e = as.factor(rb_mxi_dfs[[x]][order(rb_mxi_dfs[[x]])]), b = rep(paste0(data$rb_betaInt_prefix, x), length(rb_mxi_dfs[[x]])))
        }
      }
      categorical_ints <- names(mb_mxi_dfs[!sapply(mb_mxi_dfs,is.null)])
      
      cat_n  <- paste0("<center>N<br>",  gsub("[_]", " - ", cat_interactions), "</center>")
      cat_af <- paste0("<center>AF<br>", gsub("[_]", " - ", cat_interactions), "</center>")
      cat_colnames <- unlist(lapply(1:length(cat_n), FUN = function(x) c(cat_n[x], cat_af[x])))
      var_colnames <- c(var_colnames, cat_colnames)
      cat_n  <- paste0("N_", cat_interactions)
      cat_af <- paste0("AF_", cat_interactions)
      cat_interactions <- unlist(lapply(1:length(cat_n), FUN = function(x) c(cat_n[x], cat_af[x])))
      var_columns <- c(var_columns, cat_interactions)
      data$cat_interactions <- cat_interactions
    }
    continuous_ints <- names(mb_mxi_dfs[sapply(mb_mxi_dfs, is.null)])
    
    data$minRange_names <- paste0("_minRange_", continuous_ints)
    data$maxRange_names <- paste0("_maxRange_", continuous_ints)
    data$var_colnames <- var_colnames
    data$var_columns  <- var_columns
    data$interactions <- interactions
    data$categorical_ints <- categorical_ints
    data$continuous_ints  <- continuous_ints
    data$mb_mxi_dfs <- mb_mxi_dfs
    data$rb_mxi_dfs <- rb_mxi_dfs
    
    
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
  
  
  # Button Updates--------------------------------------------------------------
  button <- reactiveValues(gwis_clicked = FALSE)
  
  observeEvent("", {
    show("gwas_panel")
    updateButton(session,'gwas', style = "warning", icon = icon("chart-column"))
  }, once = TRUE)
  
  observeEvent(input$gwas, {
    if (!button$gwis_clicked) {
      button$gwis_clicked = TRUE
      updateButton(session,'gwas', style = "warning", icon = icon("chart-column"))
      show("gwas_panel")
    } else {
      button$gwis_clicked = FALSE
      updateButton(session,'gwas', style = "secondary", icon = icon("chart-column"))
      hide("gwas_panel")
    }
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
    
    if (!is.null(data$mh_nvar$P_Value_Marginal)) {
      dfs$mb_marginal <- get_chr_colors(data$mh_nvar$P_Value_Marginal, colors)
    }
    
    if (!is.null(data$mh_nvar$robust_P_Value_Marginal)) {
      dfs$rb_marginal <- get_chr_colors(data$mh_nvar$robust_P_Value_Marginal, colors)
    }
    
    if (!is.null(data$mh_nvar$P_Value_Interaction)) {
      dfs$mb_interaction <- get_chr_colors(data$mh_nvar$P_Value_Interaction, colors)
    }
    
    if (!is.null(data$mh_nvar$robust_P_Value_Interaction)) {
      dfs$rb_interaction <- get_chr_colors(data$mh_nvar$robust_P_Value_Interaction, colors)
    }
    
    if (!is.null(data$mh_nvar$P_Value_Joint)) {
      dfs$mb_joint <- get_chr_colors(data$mh_nvar$P_Value_Joint, colors)
    }
    
    if (!is.null(data$mh_nvar$robust_P_Value_Joint)) {
      dfs$rb_joint <- get_chr_colors(data$mh_nvar$robust_P_Value_Joint, colors)
    }
    
    return(dfs)
  })
  

  
  # Manhattan Plot -------------------------------------------------------------
  output$mb_marginal_manhattan_plot <- renderPlot({
    req(data$df)
    if (!is.null(data$mh_data$P_Value_Marginal)) {
      manhattan_plot(data$df[data$mh_data$P_Value_Marginal, c("cumulative_pos", "P_Value_Marginal")], data$x_breaks, mh_sigThreshold(), mh_sigColor(), mh_chrColor()$mb_marginal)
    } else {
      manhattan_void()
    }
  })
  
  output$rb_marginal_manhattan_plot <- renderPlot({
    req(data$df)
    if (!is.null(data$mh_data$robust_P_Value_Marginal)) {
      manhattan_plot(data$df[data$mh_data$robust_P_Value_Marginal, c("cumulative_pos", "robust_P_Value_Marginal")], data$x_breaks, mh_sigThreshold(), mh_sigColor(), mh_chrColor()$rb_marginal)
    } else {
      manhattan_void()
    }
  })
  
  output$mb_interaction_manhattan_plot <- renderPlot({
    req(data$df)
    if (!is.null(data$mh_data$P_Value_Interaction)) {
      manhattan_plot(data$df[data$mh_data$P_Value_Interaction, c("cumulative_pos", "P_Value_Interaction")], data$x_breaks, mh_sigThreshold(), mh_sigColor(), mh_chrColor()$mb_interaction)
    } else {
      manhattan_void()
    }
  })
  
  output$rb_interaction_manhattan_plot <- renderPlot({
    req(data$df)
    if (!is.null(data$mh_data$robust_P_Value_Interaction)) {
      manhattan_plot(data$df[data$mh_data$robust_P_Value_Interaction, c("cumulative_pos", "robust_P_Value_Interaction")], data$x_breaks, mh_sigThreshold(), mh_sigColor(), mh_chrColor()$rb_interaction)
    } else {
      manhattan_void()
    }
  })
  
  output$mb_joint_manhattan_plot <- renderPlot({
    req(data$df)
    if (!is.null(data$mh_data$P_Value_Joint)) {
      manhattan_plot(data$df[data$mh_data$P_Value_Joint, c("cumulative_pos", "P_Value_Joint")], data$x_breaks, mh_sigThreshold(), mh_sigColor(), mh_chrColor()$mb_joint)
    } else {
      manhattan_void()
    }
  })

  output$rb_joint_manhattan_plot <- renderPlot({
    req(data$df)
    if (!is.null(data$mh_data$robust_P_Value_Joint)) {
      manhattan_plot(data$df[data$mh_data$robust_P_Value_Joint, c("cumulative_pos", "robust_P_Value_Joint")], data$x_breaks, mh_sigThreshold(), mh_sigColor(), mh_chrColor()$rb_joint)
    } else {
      manhattan_void()
    }
  })
  
  
  
  # Manhattan Tooltip-----------------------------------------------------------
  output$mb_marginal_manhattan_plot_hover_info <- renderUI({
    req(input[["mb_marginal_manhattan_plot_hover"]], data$mh_data$P_Value_Marginal)
    manhattan_tooltip(input[["mb_marginal_manhattan_plot_hover"]], data$df[data$mh_data$P_Value_Marginal, c("CHR", "POS", "cumulative_pos", "P_Value_Marginal")])
  })

  output$rb_marginal_manhattan_plot_hover_info <- renderUI({
    req(input[["rb_marginal_manhattan_plot_hover"]], data$mh_data$robust_P_Value_Marginal)
    manhattan_tooltip(input[["rb_marginal_manhattan_plot_hover"]], data$df[data$mh_data$robust_P_Value_Marginal, c("CHR", "POS", "cumulative_pos", "robust_P_Value_Marginal")])
  })

  output$mb_interaction_manhattan_plot_hover_info <- renderUI({
    req(input[["mb_interaction_manhattan_plot_hover"]], data$mh_data$P_Value_Interaction)
    manhattan_tooltip(input[["mb_interaction_manhattan_plot_hover"]], data$df[data$mh_data$P_Value_Interaction, c("CHR", "POS", "cumulative_pos", "P_Value_Interaction")])
  })

  output$rb_interaction_manhattan_plot_hover_info <- renderUI({
    req(input[["rb_interaction_manhattan_plot_hover"]], data$mh_data$robust_P_Value_Interaction)
    manhattan_tooltip(input[["rb_interaction_manhattan_plot_hover"]], data$df[data$mh_data$robust_P_Value_Interaction, c("CHR", "POS", "cumulative_pos", "robust_P_Value_Interaction")])
  })

  output$mb_joint_manhattan_plot_hover_info <- renderUI({
    req(input[["mb_joint_manhattan_plot_hover"]], data$mh_data$P_Value_Joint)
    manhattan_tooltip(input[["mb_joint_manhattan_plot_hover"]], data$df[data$mh_data$P_Value_Joint, c("CHR", "POS", "cumulative_pos", "P_Value_Joint")])
  })
  
  output$rb_joint_manhattan_plot_hover_info <- renderUI({
    req(input[["rb_joint_manhattan_plot_hover"]], data$mh_data$robust_P_Value_Joint)
    manhattan_tooltip(input[["rb_joint_manhattan_plot_hover"]], data$df[data$mh_data$robust_P_Value_Joint, c("CHR", "POS", "cumulative_pos", "robust_P_Value_Joint")])
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
    req(data$qq$P_Value_Marginal, data$lambda$P_Value_Marginal)
    qq_plot(data$qq$P_Value_Marginal, data$lambda$P_Value_Marginal)
  })
  
  output$rb_marginal_qq_plot <- renderPlot({
    req(data$qq$robust_P_Value_Marginal, data$lambda$robust_P_Value_Marginal)
    qq_plot(data$qq$robust_P_Value_Marginal, data$lambda$robust_P_Value_Marginal)
  })
  
  output$mb_interaction_qq_plot <- renderPlot({
    req(data$qq$P_Value_Interaction, data$lambda$P_Value_Interaction)
    qq_plot(data$qq$P_Value_Interaction, data$lambda$P_Value_Interaction)
  })
  
  output$rb_interaction_qq_plot <- renderPlot({
    req(data$qq$robust_P_Value_Interaction, data$lambda$robust_P_Value_Interaction)
    qq_plot(data$qq$robust_P_Value_Interaction, data$lambda$robust_P_Value_Interaction)
  })
  
  output$mb_joint_qq_plot <- renderPlot({
    req(data$qq$P_Value_Joint, data$lambda$P_Value_Joint)
    qq_plot(data$qq$P_Value_Joint, data$lambda$P_Value_Joint)
  })
  
  output$rb_joint_qq_plot <- renderPlot({
    req(data$qq$robust_P_Value_Joint, data$lambda$robust_P_Value_Joint)
    qq_plot(data$qq$robust_P_Value_Joint, data$lambda$robust_P_Value_Joint)
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
    data$mb_marginal_nearest_points <- nearPoints(data$df, input$mb_marginal_manhattan_plot_click, xvar = "cumulative_pos", yvar = "P_Value_Marginal")
  })
  
  observeEvent(input$rb_marginal_manhattan_plot_click, {
    data$rb_marginal_nearest_points <- nearPoints(data$df, input$rb_marginal_manhattan_plot_click, xvar = "cumulative_pos", yvar = "robust_P_Value_Marginal")
  })
  
  observeEvent(input$mb_interaction_manhattan_plot_click, {
    data$mb_interaction_nearest_points <- nearPoints(data$df, input$mb_interaction_manhattan_plot_click, xvar = "cumulative_pos", yvar = "P_Value_Interaction")
  })
  
  observeEvent(input$rb_interaction_manhattan_plot_click, {
    data$rb_interaction_nearest_points <- nearPoints(data$df, input$rb_interaction_manhattan_plot_click, xvar = "cumulative_pos", yvar = "robust_P_Value_Interaction")
  })
  
  observeEvent(input$mb_joint_manhattan_plot_click, {
    data$mb_joint_nearest_points <- nearPoints(data$df, input$mb_joint_manhattan_plot_click, xvar = "cumulative_pos", yvar = "P_Value_Joint")
  })
  
  observeEvent(input$rb_joint_manhattan_plot_click, {
    data$rb_joint_nearest_points <- nearPoints(data$df, input$rb_joint_manhattan_plot_click, xvar = "cumulative_pos", yvar = "robust_P_Value_Joint")
  })
  
  
  
  # Variant Table---------------------------------------------------------------
  output$mb_marginal_manhattan_plot_table <- DT::renderDT({
    req(data$mb_marginal_nearest_points)
    variantTable(data$mb_marginal_nearest_points, "P_Value_Marginal", data$var_columns, data$var_colnames)
  })
  
  output$rb_marginal_manhattan_plot_table <- DT::renderDT({
    req(data$rb_marginal_nearest_points)
    variantTable(data$rb_marginal_nearest_points, "robust_P_Value_Marginal", data$var_columns, data$var_colnames)
  })
  
  output$mb_interaction_manhattan_plot_table <- DT::renderDT({
    req(data$mb_interaction_nearest_points)
    variantTable(data$mb_interaction_nearest_points, "P_Value_Interaction", data$var_columns, data$var_colnames)
  })
  
  output$rb_interaction_manhattan_plot_table <- DT::renderDT({
    req(data$rb_interaction_nearest_points)
    variantTable(data$rb_interaction_nearest_points, "robust_P_Value_Interaction", data$var_columns, data$var_colnames)
  })
  
  output$mb_joint_manhattan_plot_table <- DT::renderDT({
    req(data$mb_joint_nearest_points)
    variantTable(data$mb_joint_nearest_points, "P_Value_Joint", data$var_columns, data$var_colnames)
  })
  
  output$rb_joint_manhattan_plot_table <- DT::renderDT({
    req(data$rb_joint_nearest_points)
    variantTable(data$rb_joint_nearest_points, "robust_P_Value_Joint", data$var_columns, data$var_colnames)
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
    ssTable_box("rb_joint_ssTable")
  })
  
  
  
  # Variant Table Row Selected -------------------------------------------------
  observeEvent(ignoreInit = TRUE, list(data$mb_marginal_nearest_points, input$mb_marginal_manhattan_plot_table_rows_selected), {
    row <- input$mb_marginal_manhattan_plot_table_rows_selected
    ssTables(output, "mb", "marginal", data$mb_marginal_nearest_points[row, ,drop = FALSE], data$int_colnames, data$mb_beta, data$mb_se, data$mb_covs, data$mb_cov_rownames)
  })

  observeEvent(ignoreInit = TRUE, list(data$rb_marginal_nearest_points, input$rb_marginal_manhattan_plot_table_rows_selected), {
    row <- input$rb_marginal_manhattan_plot_table_rows_selected
    ssTables(output, "rb", "marginal", data$rb_marginal_nearest_points[row, ,drop = FALSE], data$int_colnames, data$rb_beta, data$rb_se, data$rb_covs, data$rb_cov_rownames)
  })

  observeEvent(ignoreInit = TRUE, list(data$mb_interaction_nearest_points, input$mb_interaction_manhattan_plot_table_rows_selected), {
    row <- input$mb_interaction_manhattan_plot_table_rows_selected
    ssTables(output, "mb", "interaction", data$mb_interaction_nearest_points[row, ,drop = FALSE], data$int_colnames, data$mb_beta, data$mb_se, data$mb_covs, data$mb_cov_rownames)
  })

  observeEvent(ignoreInit = TRUE, list(data$rb_interaction_nearest_points, input$rb_interaction_manhattan_plot_table_rows_selected), {
    row <- input$rb_interaction_manhattan_plot_table_rows_selected
    ssTables(output, "rb", "interaction", data$rb_interaction_nearest_points[row, ,drop = FALSE], data$int_colnames, data$rb_beta, data$rb_se, data$rb_covs, data$rb_cov_rownames)
  })

  observeEvent(ignoreInit = TRUE, list(data$mb_joint_nearest_points, input$mb_joint_manhattan_plot_table_rows_selected), {
    row <- input$mb_joint_manhattan_plot_table_rows_selected
    ssTables(output, "mb", "joint", data$mb_joint_nearest_points[row, ,drop = FALSE], data$int_colnames, data$mb_beta, data$mb_se, data$mb_covs, data$mb_cov_rownames)
  })

  observeEvent(ignoreInit = TRUE, list(data$rb_joint_nearest_points, input$rb_joint_manhattan_plot_table_rows_selected), {
    row <- input$rb_joint_manhattan_plot_table_rows_selected
    ssTables(output, "rb", "joint", data$rb_joint_nearest_points[row, ,drop = FALSE], data$int_colnames, data$rb_beta, data$rb_se, data$rb_covs, data$rb_cov_rownames)
  })


  
  # G Effect vs Interaction Panel-----------------------------------------------
  observeEvent(data$interactions, {
    output$mb_marginal_ssTable_mxi_ui <- renderUI({
      mxi_panel("mb_marginal", data$interactions, lapply(data$continuous_ints, rangePanel, prefix="mb_marginal"))
    })
  })
  
  observeEvent(data$interactions, {
    output$rb_marginal_ssTable_mxi_ui <- renderUI({
      mxi_panel("rb_marginal", data$interactions, lapply(data$continuous_ints, rangePanel, prefix="rb_marginal"))
    })
  })
    
  observeEvent(data$interactions, {
    output$mb_interaction_ssTable_mxi_ui <- renderUI({
      mxi_panel("mb_interaction", data$interactions, lapply(data$continuous_ints, rangePanel, prefix="mb_interaction"))
    })
  })
  
  observeEvent(data$interactions, {
    output$rb_interaction_ssTable_mxi_ui <- renderUI({
      mxi_panel("rb_interaction", data$interactions, lapply(data$continuous_ints, rangePanel, prefix="rb_interaction"))
    })
  })
  
  observeEvent(data$interactions, {
    output$mb_joint_ssTable_mxi_ui <- renderUI({
      mxi_panel("mb_joint", data$interactions, lapply(data$continuous_ints, rangePanel, prefix="mb_joint"))
    })
  })
  
  observeEvent(data$interactions, {
    output$rb_joint_ssTable_mxi_ui <- renderUI({
      mxi_panel("rb_joint", data$interactions, lapply(data$continuous_ints, rangePanel, prefix="rb_joint"))
    })
  })
  
  
  
  # G Effect vs Interaction Choice Panel----------------------------------------
  observeEvent(input$mb_marginal_ssTable_mxi_select, {
    choice <- input$mb_marginal_ssTable_mxi_select
    ret    <- toggle_rangePanel("mb_marginal", choice, data$continuous_ints)
    observeEvent(ignoreInit = TRUE, list(input[[ret[[1]]]], input[[ret[[2]]]]), {
      req(input[[ret[[1]]]], input[[ret[[2]]]])
      data$mb_mxi_dfs[[choice]] <- mxiDF(choice, input[[ret[[1]]]], input[[ret[[2]]]], data$mb_betaInt_prefix)
      data$rb_mxi_dfs[[choice]] <- mxiDF(choice, input[[ret[[1]]]], input[[ret[[2]]]], data$rb_betaInt_prefix)
      updateRangeInputs(names(input)[grepl(paste0("_minRange_", choice), names(input))], session, "Min.", input[[ret[[1]]]])
      updateRangeInputs(names(input)[grepl(paste0("_maxRange_", choice), names(input))], session, "Max.", input[[ret[[2]]]])
    })
  })
  
  observeEvent(input$rb_marginal_ssTable_mxi_select, {
    choice <- input$rb_marginal_ssTable_mxi_select
    ret    <- toggle_rangePanel("rb_marginal", choice, data$continuous_ints)
    observeEvent(ignoreInit = TRUE, list(input[[ret[[1]]]], input[[ret[[2]]]]), {
      req(input[[ret[[1]]]], input[[ret[[2]]]])
      data$mb_mxi_dfs[[choice]] <- mxiDF(choice, input[[ret[[1]]]], input[[ret[[2]]]], data$mb_betaInt_prefix)
      data$rb_mxi_dfs[[choice]] <- mxiDF(choice, input[[ret[[1]]]], input[[ret[[2]]]], data$rb_betaInt_prefix)
      updateRangeInputs(names(input)[grepl(paste0("_minRange_", choice), names(input))], session, "Min.", input[[ret[[1]]]])
      updateRangeInputs(names(input)[grepl(paste0("_maxRange_", choice), names(input))], session, "Max.", input[[ret[[2]]]])
    })
  })
  
  observeEvent(input$mb_interaction_ssTable_mxi_select, {
    choice <- input$mb_interaction_ssTable_mxi_select
    ret    <- toggle_rangePanel("mb_interaction", choice, data$continuous_ints)
    observeEvent(ignoreInit = TRUE, list(input[[ret[[1]]]], input[[ret[[2]]]]), {
      req(input[[ret[[1]]]], input[[ret[[2]]]])
      data$mb_mxi_dfs[[choice]] <- mxiDF(choice, input[[ret[[1]]]], input[[ret[[2]]]], data$mb_betaInt_prefix)
      data$rb_mxi_dfs[[choice]] <- mxiDF(choice, input[[ret[[1]]]], input[[ret[[2]]]], data$rb_betaInt_prefix)
      updateRangeInputs(names(input)[grepl(paste0("_minRange_", choice), names(input))], session, "Min.", input[[ret[[1]]]])
      updateRangeInputs(names(input)[grepl(paste0("_maxRange_", choice), names(input))], session, "Max.", input[[ret[[2]]]])
    })
  })
  
  observeEvent(input$rb_interaction_ssTable_mxi_select, {
    choice <- input$rb_interaction_ssTable_mxi_select
    ret    <- toggle_rangePanel("rb_interaction", choice, data$continuous_ints)
    observeEvent(ignoreInit = TRUE, list(input[[ret[[1]]]], input[[ret[[2]]]]), {
      req(input[[ret[[1]]]], input[[ret[[2]]]])
      data$mb_mxi_dfs[[choice]] <- mxiDF(choice, input[[ret[[1]]]], input[[ret[[2]]]], data$mb_betaInt_prefix)
      data$rb_mxi_dfs[[choice]] <- mxiDF(choice, input[[ret[[1]]]], input[[ret[[2]]]], data$rb_betaInt_prefix)
      updateRangeInputs(names(input)[grepl(paste0("_minRange_", choice), names(input))], session, "Min.", input[[ret[[1]]]])
      updateRangeInputs(names(input)[grepl(paste0("_maxRange_", choice), names(input))], session, "Max.", input[[ret[[2]]]])
    })
  })
  
  observeEvent(input$mb_joint_ssTable_mxi_select, {
    choice <- input$mb_joint_ssTable_mxi_select
    ret    <- toggle_rangePanel("mb_joint", choice, data$continuous_ints)
    observeEvent(ignoreInit = TRUE, list(input[[ret[[1]]]], input[[ret[[2]]]]), {
      req(input[[ret[[1]]]], input[[ret[[2]]]])
      data$mb_mxi_dfs[[choice]] <- mxiDF(choice, input[[ret[[1]]]], input[[ret[[2]]]], data$mb_betaInt_prefix)
      data$rb_mxi_dfs[[choice]] <- mxiDF(choice, input[[ret[[1]]]], input[[ret[[2]]]], data$rb_betaInt_prefix)
      updateRangeInputs(names(input)[grepl(paste0("_minRange_", choice), names(input))], session, "Min.", input[[ret[[1]]]])
      updateRangeInputs(names(input)[grepl(paste0("_maxRange_", choice), names(input))], session, "Max.", input[[ret[[2]]]])
    })
  })

  observeEvent(input$rb_joint_ssTable_mxi_select, {
    choice <- input$rb_joint_ssTable_mxi_select
    ret    <- toggle_rangePanel("rb_joint", choice, data$continuous_ints)
    observeEvent(ignoreInit = TRUE, list(input[[ret[[1]]]], input[[ret[[2]]]]), {
      req(input[[ret[[1]]]], input[[ret[[2]]]])
      data$mb_mxi_dfs[[choice]] <- mxiDF(choice, input[[ret[[1]]]], input[[ret[[2]]]], data$mb_betaInt_prefix)
      data$rb_mxi_dfs[[choice]] <- mxiDF(choice, input[[ret[[1]]]], input[[ret[[2]]]], data$rb_betaInt_prefix)
      updateRangeInputs(names(input)[grepl(paste0("_minRange_", choice), names(input))], session, "Min.", input[[ret[[1]]]])
      updateRangeInputs(names(input)[grepl(paste0("_maxRange_", choice), names(input))], session, "Max.", input[[ret[[2]]]])
    })
  })
  
  

  # G Effect vs Interaction Plot -----------------------------------------------
  observeEvent(ignoreInit = TRUE, list(data$mb_marginal_nearest_points, input$mb_marginal_manhattan_plot_table_rows_selected, input$mb_marginal_ssTable_mxi_select), {
    row    <- input$mb_marginal_manhattan_plot_table_rows_selected
    choice <- input$mb_marginal_ssTable_mxi_select
    output$mb_marginal_ssTable_mxi <- renderPlot({mxi_plot(data$mb_marginal_nearest_points[row, ], data$mb_mxi_dfs[[choice]], choice, data$mb_betaG_prefix)})
  })

  observeEvent(ignoreInit = TRUE, list(data$rb_marginal_nearest_points, input$rb_marginal_manhattan_plot_table_rows_selected, input$rb_marginal_ssTable_mxi_select), {
    row    <- input$rb_marginal_manhattan_plot_table_rows_selected
    choice <- input$rb_marginal_ssTable_mxi_select
    output$rb_marginal_ssTable_mxi <- renderPlot({mxi_plot(data$rb_marginal_nearest_points[row, ], data$rb_mxi_dfs[[choice]], choice, data$rb_betaG_prefix)})
  })

  observeEvent(ignoreInit = TRUE, list(data$mb_interaction_nearest_points, input$mb_interaction_manhattan_plot_table_rows_selected, input$mb_interaction_ssTable_mxi_select), {
    row    <- input$mb_interaction_manhattan_plot_table_rows_selected
    choice <- input$mb_interaction_ssTable_mxi_select
    output$mb_interaction_ssTable_mxi <- renderPlot({mxi_plot(data$mb_interaction_nearest_points[row, ], data$mb_mxi_dfs[[choice]], choice, data$mb_betaG_prefix)})
  })

  observeEvent(ignoreInit = TRUE, list(data$rb_interaction_nearest_points, input$rb_interaction_manhattan_plot_table_rows_selected, input$rb_interaction_ssTable_mxi_select), {
    row    <- input$rb_interaction_manhattan_plot_table_rows_selected
    choice <-  input$rb_interaction_ssTable_mxi_select
    output$rb_interaction_ssTable_mxi <- renderPlot({mxi_plot(data$rb_interaction_nearest_points[row, ], data$rb_mxi_dfs[[choice]], choice, data$rb_betaG_prefix)})
  })

  observeEvent(ignoreInit = TRUE, list(data$mb_joint_nearest_points, input$mb_joint_manhattan_plot_table_rows_selected, input$mb_joint_ssTable_mxi_select), {
    row    <- input$mb_joint_manhattan_plot_table_rows_selected
    choice <- input$mb_joint_ssTable_mxi_select
    output$mb_joint_ssTable_mxi <- renderPlot({mxi_plot(data$mb_joint_nearest_points[row, ], data$mb_mxi_dfs[[choice]], choice, data$mb_betaG_prefix)})
  })

  observeEvent(ignoreInit = TRUE, list(data$rb_joint_nearest_points, input$rb_joint_manhattan_plot_table_rows_selected, input$rb_joint_ssTable_mxi_select), {
    row    <- input$rb_joint_manhattan_plot_table_rows_selected
    choice <- input$rb_joint_ssTable_mxi_select
    output$rb_joint_ssTable_mxi <- renderPlot({mxi_plot(data$rb_joint_nearest_points[row, ], data$rb_mxi_dfs[[choice]], choice, data$rb_betaG_prefix)})
  })
}