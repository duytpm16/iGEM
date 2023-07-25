server <- function(input, output, session) {

  # UI - GENERAL ---------------------------------------------------------------
  data <- reactiveValues(df = NULL)
  
  mh_data <- reactiveValues(mb_marginal    = NULL,
                            rb_marginal    = NULL,
                            mb_interaction = NULL,
                            rb_interaction = NULL,
                            mb_joint       = NULL,
                            rb_joint       = NULL)
  
  mh_nvar <- reactiveValues(mb_marginal    = NULL,
                            rb_marginal    = NULL,
                            mb_interaction = NULL,
                            rb_interaction = NULL,
                            mb_joint       = NULL,
                            rb_joint       = NULL)
  
  # READ INPUT FILE ------------------------------------------------------------
  observeEvent(input$inputFile, {
    req(input$inputFile)
    
    ## Read File
    df   <- fread(input$inputFile$datapath, sep = "\t", header = T)
    nvar <- nrow(df)
    coln <- colnames(df)
    gc(verbose = FALSE)
    
    setorderv(df, c("CHR", "POS"))
    chrom_lengths     <- chrom_lengths_hg38[df[, kit::funique(CHR)]]
    cum_chrom_lengths <- get_cumulative_length(chrom_lengths)
    df[, cumulative_pos := POS + cum_chrom_lengths[CHR]]
    gc(verbose = FALSE)
  
    ## P-values------------------------------------------------------------------
    pvalue_columns <- c("P_Value_Marginal",    "robust_P_Value_Marginal",
                        "P_Value_Interaction", "robust_P_Value_Interaction",
                        "P_Value_Joint",       "robust_P_Value_Joint")
    pvalue_columns <- pvalue_columns[pvalue_columns %in% coln]
    
    
    #genomic_inflation <- df[, lapply(.SD, function(x) median(qchisq(x, 1, lower.tail = FALSE)) / qchisq(0.5, 1)), .SDcols=pvalue_columns]
    #genomic_inflation2 <- df[, lapply(.SD, function(x) median(qnorm(x, lower.tail = FALSE)^2) / qchisq(0.5, 1)), .SDcols=pvalue_columns]
    df[, (pvalue_columns) := lapply(.SD, function(x) -log10(x)), .SDcols = pvalue_columns]
    data$qq <- sapply(pvalue_columns,  function(x) fastqq::drop_dense(df[[x]], -log10(stats::ppoints(nvar))), simplify = FALSE, USE.NAMES = TRUE)
    gc(verbose = FALSE)
    
    
    index <- vector(mode = "logical", length = nvar)
    round_cumpos <- plyr::round_any(df[["cumulative_pos"]], accuracy = 100000)
    
    for (pcol in pvalue_columns) {
      pcol_split <- strsplit(pcol, "_")[[1]]
      se   <- ifelse(pcol_split[1] == "P", "mb", "rb")
      test <- tolower(pcol_split[length(pcol_split)])
      
      if (pcol %in% coln) {
        sig_df <- data.table(pval = round(df[[pcol]], digits = 3), cpos = round_cumpos)
        
        index[sig_df[["pval"]] > 8] <- TRUE
        keep <- !fduplicated(sig_df) & (sig_df[["pval"]] > 5)
        keep[!(fduplicated(sig_df[["pval"]]) & fduplicated(round_cumpos)) & sig_df[["pval"]] <= 5] <- TRUE
        
        tmp <- sum(keep)
        if ((nvar > 125000) & (tmp < 125000)) {
          keep[sample(which(!keep & sig_df[["pval"]] <= 5), 125000 - tmp, replace = FALSE)] <- TRUE
        } else if (tmp > 125000) {
          keep[sample(which(keep & sig_df[["pval"]] <= 5), tmp - 125000, replace = FALSE)] <- FALSE
        } else if (nvar <= 125000) {
          keep[1:nvar] <- TRUE 
        }
        index[keep] <- TRUE
        
        mh_data[[paste0(se, "_", test)]] <- c(1:nvar)[keep] 
        mh_nvar[[paste0(se, "_", test)]] <- as.data.frame(df[mh_data[[paste0(se, "_", test)]], .N, by = CHR])

        gc(verbose = FALSE)
      }
    }
    rm(keep)
    rm(sig_df)
    
    new_index  <- c(1:nvar)[keep]
    new_index2 <- c(1:length(new_index))
    df <- df[index, ]
    rm(index)
    gc(verbose = FALSE)
    
    
    # Get interactions including marginal columns
    interactions <- c("Beta_Marginal", "Beta_G", coln[grepl("^Beta_G-", coln)])
    interactions <- gsub("Beta_", "", interactions)
    int_colnames <- c("Marginal", "Main", gsub("G[-]", "G x ", interactions[grepl("^G[-]", interactions)]))

    
    ## Betas and SE-------------------------------------------------------------
    data$mb_beta <- paste0("Beta_", interactions)
    data$rb_beta <- data$mb_beta
    if(any(grepl("^robust_Beta", coln))){
      data$rb_beta <- paste0("robust_Beta_", interactions)
    }

    data$mb_se   <- paste0("SE_Beta_", interactions)
    data$rb_se   <- paste0("robust_SE_Beta_", interactions)
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
    data$interactions <- gsub("G[-]", "", interactions[-c(1,2)])
    cat_interactions  <- data$interactions
    cat_interactions  <- coln[grepl(paste0("^N[_]", cat_interactions, collapse = "|"), coln)]
    
    
    if (length(cat_interactions) == 0) {
      data$mxi_df <- data.frame(i = unlist(lapply(data$interactions, FUN = function(x) rep(x, 7))),
                                e = rep(-3:3, length(data$interactions)))
      data$mxi_df$b <- paste0("Beta_G-", data$mxi_df$i)
    } else {
      cat_interactions <- gsub("N[_]", "", cat_interactions)
      
      tmp <- vector(mode = "list")
      for (cat in cat_interactions){
        cat_split <- strsplit(cat, "_")[[1]]
        stratum <- as.numeric(cat_split[length(cat_split)])
        intname <- paste0(cat_split[-length(cat_split)], collapse = "_")
        tmp[[intname]] <- c(tmp[[intname]], stratum)
      }
      
      tmp_i <- c()
      tmp_e <- c()
      for (x in data$interactions) {
        if (x %in% names(tmp)) {
          tmp_e <- c(tmp_e, c(0, tmp[[x]]))
          tmp_i <- c(tmp_i, rep(x, length(tmp[[x]]) + 1))
        } else {
          tmp_e <- c(tmp_e, -3:3)
          tmp_i <- c(tmp_i, rep(x, 7))
        }
      }
      data$mxi_df <- data.frame(i = tmp_i, e = tmp_e, b = paste0("Beta_G-", tmp_i))
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
    data$var_colnames <- var_colnames
    data$var_columns <- var_columns
    
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
  button <- reactiveValues(gwis_clicked   = FALSE,
                           mainxe_clicked = FALSE)
  
  observeEvent("", {
    show("gwas_panel")
    updateButton(session,'gwas', style = "warning", icon = icon("chart-column"))
  }, once = TRUE)
  
  observeEvent(input$gwas, {
    if (!button$gwis_clicked) {
      button$gwis_clicked = TRUE
      updateButton(session,'gwas',  style = "warning",   icon = icon("chart-column"))
      show("gwas_panel")
    } else {
      button$gwis_clicked = FALSE
      updateButton(session,'gwas',  style = "secondary",   icon = icon("chart-column"))
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
    
    if (!is.null(mh_nvar$mb_marginal)) {
      dfs$mb_marginal <- get_chr_colors(mh_nvar$mb_marginal, colors)
    }
    
    if (!is.null(mh_nvar$rb_marginal)) {
      dfs$rb_marginal <- get_chr_colors(mh_nvar$rb_marginal, colors)
    }
    
    if (!is.null(mh_nvar$mb_interaction)) {
      dfs$mb_interaction <- get_chr_colors(mh_nvar$mb_interaction, colors)
    }
    
    if (!is.null(mh_nvar$rb_interaction)) {
      dfs$rb_interaction <- get_chr_colors(mh_nvar$rb_interaction, colors)
    }
    
    if (!is.null(mh_nvar$mb_joint)) {
      dfs$mb_joint <- get_chr_colors(mh_nvar$mb_joint, colors)
    }
    
    if (!is.null(mh_nvar$rb_joint)) {
      dfs$rb_joint <- get_chr_colors(mh_nvar$rb_joint, colors)
    }
    
    return(dfs)
  })
  

  
  # Manhattan Plot -------------------------------------------------------------
  output$mb_marginal_manhattan_plot <- renderPlot({
    manhattan_plot(data$df[mh_data$mb_marginal, c("cumulative_pos", "P_Value_Marginal")], data$x_breaks, mh_sigThreshold(), mh_sigColor(), mh_chrColor()$mb_marginal)
  })
  
  output$rb_marginal_manhattan_plot <- renderPlot({
    manhattan_plot(data$df[mh_data$rb_marginal, c("cumulative_pos", "robust_P_Value_Marginal")], data$x_breaks, mh_sigThreshold(), mh_sigColor(), mh_chrColor()$rb_marginal)
  })
  
  output$mb_interaction_manhattan_plot <- renderPlot({
    manhattan_plot(data$df[mh_data$mb_interaction, c("cumulative_pos", "P_Value_Interaction")], data$x_breaks, mh_sigThreshold(), mh_sigColor(), mh_chrColor()$mb_interaction)
  })
  
  output$rb_interaction_manhattan_plot <- renderPlot({
    manhattan_plot(data$df[mh_data$rb_interaction, c("cumulative_pos", "robust_P_Value_Interaction")], data$x_breaks, mh_sigThreshold(), mh_sigColor(), mh_chrColor()$rb_interaction)
  })
  
  output$mb_joint_manhattan_plot <- renderPlot({
    manhattan_plot(data$df[mh_data$mb_joint, c("cumulative_pos", "P_Value_Joint")], data$x_breaks, mh_sigThreshold(), mh_sigColor(), mh_chrColor()$mb_joint)
  })

  output$rb_joint_manhattan_plot <- renderPlot({
    manhattan_plot(data$df[mh_data$rb_joint, c("cumulative_pos", "robust_P_Value_Joint")], data$x_breaks, mh_sigThreshold(), mh_sigColor(), mh_chrColor()$rb_joint)
  })
  
  
  
  # Manhattan Tooltip-----------------------------------------------------------
  output$mb_marginal_manhattan_plot_hover_info <- renderUI({
    manhattan_tooltip(input[["mb_marginal_manhattan_plot_hover"]], data$df[mh_data$mb_marginal, c("CHR", "POS", "cumulative_pos", "P_Value_Marginal")])
  })

  output$rb_marginal_manhattan_plot_hover_info <- renderUI({
    manhattan_tooltip(input[["rb_marginal_manhattan_plot_hover"]], data$df[mh_data$rb_marginal, c("CHR", "POS", "cumulative_pos", "robust_P_Value_Marginal")])
  })

  output$mb_interaction_manhattan_plot_hover_info <- renderUI({
    manhattan_tooltip(input[["mb_interaction_manhattan_plot_hover"]], data$df[mh_data$mb_interaction, c("CHR", "POS", "cumulative_pos", "P_Value_Interaction")])
  })

  output$rb_interaction_manhattan_plot_hover_info <- renderUI({
    manhattan_tooltip(input[["rb_interaction_manhattan_plot_hover"]], data$df[mh_data$rb_interaction, c("CHR", "POS", "cumulative_pos", "robust_P_Value_Interaction")])
  })

  output$mb_joint_manhattan_plot_hover_info <- renderUI({
      manhattan_tooltip(input[["mb_joint_manhattan_plot_hover"]], data$df[mh_data$mb_joint, c("CHR", "POS", "cumulative_pos", "P_Value_Joint")])
  })
  
  output$rb_joint_manhattan_plot_hover_info <- renderUI({
    manhattan_tooltip(input[["rb_joint_manhattan_plot_hover"]], data$df[mh_data$rb_joint, c("CHR", "POS", "cumulative_pos", "robust_P_Value_Joint")])
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
    qq_plot(data$qq$P_Value_Marginal, "")
  })
  
  output$rb_marginal_qq_plot <- renderPlot({
    qq_plot(data$qq$robust_P_Value_Marginal, "")
  })
  
  output$mb_interaction_qq_plot <- renderPlot({
    qq_plot(data$qq$P_Value_Interaction, "")
  })
  
  output$rb_interaction_qq_plot <- renderPlot({
    qq_plot(data$qq$robust_P_Value_Interaction, "")
  })
  
  output$mb_joint_qq_plot <- renderPlot({
    qq_plot(data$qq$P_Value_Joint, "")
  })
  
  output$rb_joint_qq_plot <- renderPlot({
    qq_plot(data$qq$robust_P_Value_Joint, "")
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
    ssTable_box( "rb_joint_ssTable")
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

  # G Effect vs Interaction Select--------------------------------------------
  observeEvent(data$interactions, {
    updateSelectInput(session, "mb_marginal_ssTable_mxi_select",
                      label   = "Select Interaction(s):",
                      choices = data$interactions)
  })

  observeEvent(data$interactions, {
    updateSelectInput(session, "rb_marginal_ssTable_mxi_select",
                      label   = "Select Interaction(s):",
                      choices = data$interactions)
  })

  observeEvent(data$interactions, {
    updateSelectInput(session, "mb_interaction_ssTable_mxi_select",
                      label   = "Select Interaction(s):",
                      choices = data$interactions)
  })

  observeEvent(data$interactions, {
    updateSelectInput(session, "rb_interaction_ssTable_mxi_select",
                      label   = "Select Interaction(s):",
                      choices = data$interactions)
  })

  observeEvent(data$interactions, {
    updateSelectInput(session, "mb_joint_ssTable_mxi_select",
                      label   = "Select Interaction(s):",
                      choices = data$interactions)
  })

  observeEvent(data$interactions, {
    updateSelectInput(session, "rb_joint_ssTable_mxi_select",
                      label   = "Select Interaction(s):",
                      choices = data$interactions)
  })


  # G Effect vs Interaction Plot --------------------------------------------
  observeEvent(ignoreInit = TRUE, list(data$mb_marginal_nearest_points, input$mb_marginal_manhattan_plot_table_rows_selected, input$mb_marginal_ssTable_mxi_select), {
    row <- input$mb_marginal_manhattan_plot_table_rows_selected
    output$mb_marginal_ssTable_mxi <- renderPlot({mxi_plot(data$mb_marginal_nearest_points[row, ], data$mxi_df, input$mb_marginal_ssTable_mxi_select)})
  })

  observeEvent(ignoreInit = TRUE, list(data$rb_marginal_nearest_points, input$rb_marginal_manhattan_plot_table_rows_selected, input$rb_marginal_ssTable_mxi_select), {
    row <- input$rb_marginal_manhattan_plot_table_rows_selected
    output$rb_marginal_ssTable_mxi <- renderPlot({mxi_plot(data$rb_marginal_nearest_points[row, ], data$mxi_df, input$rb_marginal_ssTable_mxi_select)})
  })

  observeEvent(ignoreInit = TRUE, list(data$mb_interaction_nearest_points, input$mb_interaction_manhattan_plot_table_rows_selected, input$mb_interaction_ssTable_mxi_select), {
    row <- input$mb_interaction_manhattan_plot_table_rows_selected
    output$mb_interaction_ssTable_mxi <- renderPlot({mxi_plot(data$mb_interaction_nearest_points[row, ], data$mxi_df, input$mb_interaction_ssTable_mxi_select)})
  })

  observeEvent(ignoreInit = TRUE, list(data$rb_interaction_nearest_points, input$rb_interaction_manhattan_plot_table_rows_selected, input$rb_interaction_ssTable_mxi_select), {
    row <- input$rb_interaction_manhattan_plot_table_rows_selected
    output$rb_interaction_ssTable_mxi <- renderPlot({mxi_plot(data$rb_interaction_nearest_points[row, ], data$mxi_df, input$rb_interaction_ssTable_mxi_select)})
  })

  observeEvent(ignoreInit = TRUE, list(data$mb_joint_nearest_points, input$mb_joint_manhattan_plot_table_rows_selected, input$mb_joint_ssTable_mxi_select), {
    row <- input$mb_joint_manhattan_plot_table_rows_selected
    output$mb_joint_ssTable_mxi <- renderPlot({mxi_plot(data$mb_joint_nearest_points[row, ], data$mxi_df, input$mb_joint_ssTable_mxi_select)})
  })

  observeEvent(ignoreInit = TRUE, list(data$rb_joint_nearest_points, input$rb_joint_manhattan_plot_table_rows_selected, input$rb_joint_ssTable_mxi_select), {
    row <- input$rb_joint_manhattan_plot_table_rows_selected
    output$rb_joint_ssTable_mxi <- renderPlot({mxi_plot(data$rb_joint_nearest_points[row, ], data$mxi_df, input$rb_joint_ssTable_mxi_select)})
  })
}