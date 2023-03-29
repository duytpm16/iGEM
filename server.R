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
                      choices = c("", "GWAS", "Diagnostics", "GWAS Marginal", "GWAS Interaction", "GWAS Joint"),
                      label = "",
                      selected = x
    )
  }
  
  observeEvent(input$gwas, {
    update_all("GWAS")
  })
  observeEvent(input$gwas_marginal, {
    update_all("GWAS Marginal")
  })
  observeEvent(input$gwas_interaction, {
    update_all("GWAS Interaction")
  })
  observeEvent(input$gwas_joint, {
    update_all("GWAS Joint")
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
    hide("gwas_marginal_panel")
    hide("gwas_interaction_panel")
    hide("gwas_joint_panel")
  }, once = TRUE)
  
  observeEvent(input$gwas, {
    show("gwas_panel")
    hide("diagnostics_panel")
  })
  observeEvent(input$diagnostics, {
    show("diagnostics_panel")
    hide("patients_panel")
  })
  
  observeEvent(input$gwas_marginal, {
    show("gwas_panel")
    show("gwas_marginal_panel")
    hide("diagnostics_panel")
    hide("gwas_interaction_panel")
    hide("gwas_joint_panel")
  })
  
  observeEvent(input$gwas_interaction, {
    show("gwas_panel")
    show("gwas_interaction_panel")
    hide("gwas_marginal_panel")
    hide("gwas_joint_panel")
    hide("diagnostics_panel")
    
  })
  observeEvent(input$gwas_joint, {
    show("gwas_panel")
    show("gwas_joint_panel")
    hide("gwas_marginal_panel")
    hide("gwas_interaction_panel")
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
    var_colnames <- c("SNP ID", "CHROM", "POS", "<center>NON-EFFECT<br>Allele</center>", "<center>EFFECT<br>ALLELE</center>", "<center>N<br>SAMPLES</center>", "AF", cat_colnames)
    
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
  
  # UI - PATIENTS - 1 ----------------------------------------------------------
  
  output$box_pat <- renderUI({
      box(
        title = p("Manhattan Plot", style = 'font-size:21px;'),
        status = "primary",
        collapsible = FALSE,
        solidHeader = FALSE,
        width = 12,
        withSpinner(
          plotOutput("mb_marginal_manhattan_plot", height = 300,
                     click = "mb_marginal_manhattan_plot_click")
        )
      )
  })
  
  observeEvent(input$gwas_marginal, {

    y.max <- floor(max(data$df$P_Value_Marginal)) + 5

    mhplot <- ggplot(data$df, aes(x=cumulative_pos, y=P_Value_Marginal, color=color)) +
                geom_point() +
                ggtitle("") +
                xlab("Chromosome") +
                ylab(expression(-log[10](italic(p)))) +
                scale_x_continuous(expand = c(0.01,0),
                                   breaks = data$x_breaks,
                                   labels = names(data$x_breaks)) +
                scale_y_continuous(expand = c(0.01,0), limits = c(0, y.max)) +
                scale_color_manual(values = data$color_map,
                                   guide  = 'none') +
                theme(panel.background = element_blank(),
                      panel.grid       = element_line(color = "grey97"),
                      axis.line        = element_line(linewidth = 0.6),
                      axis.title       = element_text(size = 13),
                      axis.text        = element_text(size = 11),
                      legend.position = "none")
    output$mb_marginal_manhattan_plot <- renderPlot({
      mhplot
    })
  })
  
  
  # UI - PATIENTS - 2 -------------------------------------------------------

  output$box_pat2 <- renderUI({
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
  

  # UI - PATIENTS - 3 -------------------------------------------------------

  observeEvent(input$mb_marginal_manhattan_plot_click, {
    output$box_pat3 <- renderUI({
      box(
        title = p("Variants Table", style = 'font-size:21px;'),
        status = "primary",
        collapsible = FALSE,
        solidHeader = FALSE,
        width = 12,
        withSpinner(
          dataTableOutput("mb_marginal_manhattan_plot_table", height = 360)
        )
      )
    })
    
    data$mb_mrg_nearest_points <- nearPoints(data$df, input$mb_marginal_manhattan_plot_click,
                                             xvar = "cumulative_pos", yvar = "P_Value_Marginal")
    
    output$mb_marginal_manhattan_plot_table <- DT::renderDT({
      DT::datatable(
        data$mb_mrg_nearest_points[, c("SNPID", "CHR", "POS", "Non_Effect_Allele", "Effect_Allele", "N_Samples", "AF", data$cat_interactions)],
        colnames = data$var_colnames,
        rownames = FALSE,
        style = "bootstrap",
        selection = 'single',
        escape = FALSE,
        caption = "Select a state to examine it in more detail",
        options = list(
          dom = 'tp',
          search = list(regex = TRUE, caseInsensitive = TRUE),
          pageLength = 5,
          ordering = TRUE,
          stateSave = TRUE,
          columnDefs = list(list(targets = "_all", className = "dt-center"))
        )
      )
    })
  })
  
  
  # UI - PATIENTS - 4 -------------------------------------------------------
  
  observeEvent(input$mb_marginal_manhattan_plot_table_rows_selected, {
    row <- input$mb_marginal_manhattan_plot_table_rows_selected
    
    
    output$box_pat4 <- renderUI({
      box(
        title = p(data$mb_mrg_nearest_points$SNPID[row], style = 'font-size:21px;text-underline-position: under;text-decoration: underline;'),
        status = "primary",
        collapsible = FALSE,
        solidHeader = FALSE,
        width = 12,
        grid_page(
          layout = c(
            "       400px   300px",
            "175px  table2  table3",
            "150px  table4  table3"
          ),
          grid_card(
            "table2",
            card_body(
              dataTableOutput("mb_marginal_manhattan_plot_table2")
            )
          ),
          grid_card(
            "table3",
            card_body(
              dataTableOutput("mb_marginal_manhattan_plot_table3")
            )
          ),
          grid_card(
            "table4",
            card_body(
              dataTableOutput("mb_marginal_manhattan_plot_table4")
            )
          )
        )
      )
    })
    
    
    beta_se <- list()
    beta_se[["beta"]] <- data$mb_mrg_nearest_points[row, data$beta_columns, drop = F]
    beta_se[["se"]] <- data$mb_mrg_nearest_points[row, data$se_columns, drop = F]
    beta_se <- rbindlist(beta_se, use.names = FALSE)
    beta_se <- signif(beta_se, digits = 6)
    output$mb_marginal_manhattan_plot_table2 <- DT::renderDT({
      DT::datatable(
        beta_se,
        colnames = data$int_colnames,
        rownames = c("Coefficients", "Std. Errors"),
        style = "bootstrap",
        selection = 'none',
        caption = "Table 1. Coefficient Estimates and Standard Errors",
        options = list(
          dom = 't',
          pageLength = 2,
          scrollX = TRUE,
          columnDefs = list(list(width = '300px', targets = "_all", className = "dt-center"))
        )
      )
    })

    print(data$covariances)
    covs <- cbind(data$cov_rownames, c(data$mb_mrg_nearest_points[row, data$covariances, drop = T]))
    print(covs)
    output$mb_marginal_manhattan_plot_table3 <- DT::renderDT({
      DT::datatable(
        covs,
        colnames = c("", "Covariances"),
        rownames = FALSE,
        style = "bootstrap",
        selection = 'none',
        caption = "Table 2. Model-based Covariances",
        options = list(
          dom = 't',
          scrollX = TRUE,
          scrollY = "auto",
          columnDefs = list(list(width = '150px', targets = "_all", className = "dt-center"))
        )
      )
    })

    pvals <- data$mb_mrg_nearest_points[row, c("P_Value_Marginal", "P_Value_Interaction", "P_Value_Joint")]
    pvals <- signif(10^-pvals, digits = 6)
    output$mb_marginal_manhattan_plot_table4 <- DT::renderDT({
      DT::datatable(
        pvals,
        rownames = "P-Values",
        colnames = c("Marginal", "Interaction", "Joint"),
        style = "bootstrap",
        selection = 'none',
        caption = "Table 3. Model-based P-Values",
        options = list(
          dom = 't',
          pageLength = 1,
          scrollX = TRUE,
          columnDefs = list(list(width = '150px', targets = "_all", className = "dt-center"))
        )
      )
    })
    
    
    
  })
  
  # # UI - DIAGNOSTICS - 1 ------------------------------------------------------------------
  # 
  # output$box5 <- renderUI({
  #   div(
  #     style = "position: relative",
  #     tabBox(
  #       id = "box5",
  #       width = NULL,
  #       height = 400,
  #       tabPanel(
  #         title = "Diagnostics in selected patients",
  #         div(
  #           style = "position: absolute; left: 0.5em; bottom: 0.5em;",
  #           dropdown(
  #             radioGroupButtons(
  #               inputId = "box5.1",
  #               label = "Change time", 
  #               choiceNames = c("Year", "Quarter", "Month"), 
  #               choiceValues = c("year", "yearquarter_adm", "yearmonth_adm"), 
  #               selected = "year", 
  #               direction = "vertical"
  #             ),
  #             radioGroupButtons(
  #               inputId = "box5.2",
  #               label = "Change plot", 
  #               choiceNames = c("Count", "Proportion"), 
  #               choiceValues = c("dodge", "fill"), 
  #               selected = "dodge", 
  #               direction = "vertical"
  #             ),
  #             size = "xs",
  #             icon = icon("gear", class = "opt"), 
  #             up = TRUE
  #           )
  #         ),
  #         div(
  #           style = "position: absolute; left: 4em;bottom: 0.5em;",
  #           dropdown(
  #             downloadButton(outputId = "down_box_5", label = "Download plot"),
  #             size = "xs",
  #             icon = icon("download", class = "opt"), 
  #             up = TRUE
  #           )
  #         ),
  #         withSpinner(
  #           plotOutput("plot_dia_adm", height = 300),
  #           type = 4,
  #           color = "#d33724",
  #           size = 0.7
  #         )
  #       ),
  #       tabPanel(
  #         title = "Timing of selected diagnostics",
  #         div(
  #           style = "position:absolute;left:0.5em;bottom: 0.5em;",
  #           dropdown(
  #             downloadButton(outputId = "down_box_6", label = "Download plot"),
  #             size = "xs",
  #             icon = icon("download", class = "opt"), 
  #             up = TRUE
  #           )
  #         ),
  #         withSpinner(
  #           plotOutput("plot_dia_timing", height = 300),
  #           type = 4,
  #           color = "#d33724",
  #           size = 0.7
  #         )
  #       ),
  #       div(
  #         style = "position: absolute; right: 0.5em; bottom: 0.5em;",
  #         conditionalPanel(
  #           "input.box5 == 'Diagnostics in selected patients'",
  #           actionBttn(
  #             inputId = "dia_adm",
  #             icon = icon("search-plus", class = "opt"),
  #             style = "fill",
  #             color = "danger",
  #             size = "xs"
  #           )
  #         )
  #       ),
  #       div(
  #         style = "position: absolute; right: 0.5em; bottom: 0.5em;",
  #         conditionalPanel(
  #           "input.box5 == 'Timing of selected diagnostics'",
  #           actionBttn(
  #             inputId = "dia_timing",
  #             icon = icon("search-plus", class = "opt"),
  #             style = "fill",
  #             color = "danger",
  #             size = "xs"
  #           )
  #         )
  #       )
  #     )
  #   )
  # })
  # 
  # observeEvent((input$dia_adm), {
  #   showModal(modalDialog(
  #     renderPlot({
  #       dia_adm() + theme(
  #         axis.title = element_text(size = 20),
  #         text = element_text(size = 20),
  #         plot.title = element_text(size = 26)
  #       )
  #     }, height = 600),
  #     easyClose = TRUE,
  #     size = "l",
  #     footer = NULL
  #   ))
  # })
  # 
  # observeEvent((input$dia_timing), {
  #   showModal(modalDialog(
  #     renderPlot({
  #       plot_dia_timing() + theme(
  #         axis.title = element_text(size = 20),
  #         text = element_text(size = 20),
  #         plot.title = element_text(size = 26)
  #       )
  #     }, height = 600),
  #     easyClose = TRUE,
  #     size = "l",
  #     footer = NULL
  #   ))
  # })
  # 
  # # UI - DIAGNOSTICS - 2 ------------------------------------------------------------------
  # 
  # output$box6 <- renderUI({
  #   div(
  #     style = "position: relative",
  #     tabBox(
  #       id = "box6",
  #       width = NULL,
  #       height = 400,
  #       tabPanel(
  #         title = "Diagnostics in relation",
  #         div(
  #           style = "position: absolute; left: 0.5em; bottom: 0.5em;",
  #           dropdown(
  #             radioGroupButtons(
  #               inputId = "box6.1",
  #               label = "Select group", 
  #               choiceNames = c("Antimicrobial - Groups", "Antimicrobials", "Year", "Specialty", "Subspecialty", "Origin"),
  #               choiceValues = c("ab_group", "ab_type", "year", "specialty", "sub_specialty", "adm_route"), 
  #               selected = "year", 
  #               direction = "vertical"
  #             ),
  #             size = "xs",
  #             icon = icon("gear", class = "opt"), 
  #             up = TRUE
  #           )
  #         ),
  #         div(
  #           style = "position:absolute;left:4em;bottom: 0.5em;",
  #           dropdown(
  #             downloadButton(outputId = "down_box_7", label = "Download plot"),
  #             size = "xs",
  #             icon = icon("download", class = "opt"), 
  #             up = TRUE
  #           )
  #         ),
  #         withSpinner(
  #           plotOutput("plot_dia_perform", height = 300),
  #           type = 4,
  #           color = "#d33724",
  #           size = 0.7
  #         )
  #       ),
  #       tabPanel(
  #         title = "Table - Proportion performed",
  #         div(
  #           style = "position: absolute; left: 0.5em; bottom: 0.5em;",
  #           dropdown(
  #             radioGroupButtons(
  #               inputId = "box6.3",
  #               label = "Select group", 
  #               choiceNames = c("Antimicrobials", "Antimicrobial - Groups", "Year", "Specialty", "Subspecialty", "Origin"),
  #               choiceValues = c("ab_group", "ab_type", "year", "specialty", "sub_specialty", "adm_route"), 
  #               selected = "year", 
  #               direction = "vertical"
  #             ),
  #             size = "xs",
  #             icon = icon("gear", class = "opt"), 
  #             up = TRUE
  #           )
  #         ),
  #         withSpinner(
  #           DT::dataTableOutput("dia_table"),
  #           type = 4,
  #           color = "#d33724",
  #           size = 0.7
  #         )
  #       ),
  #       div(
  #         style = "position:absolute;right:0.5em;bottom: 0.5em;",
  #         conditionalPanel(
  #           "input.box6 == 'Diagnostics in relation'",
  #           actionBttn(
  #             inputId = "dia_perform",
  #             icon = icon("search-plus", class = "opt"),
  #             style = "fill",
  #             color = "danger",
  #             size = "xs"
  #           )
  #         )
  #       )
  #     )
  #   )
  # })
  # 
  # observeEvent((input$dia_perform), {
  #   showModal(modalDialog(
  #     renderPlot({ 
  #       plot_dia_perform() + theme(
  #         axis.title = element_text(size = 20),
  #         text = element_text(size = 20),
  #         plot.title = element_text(size = 26)
  #       )
  #     }),
  #     easyClose = TRUE,
  #     size = "l",
  #     footer = NULL
  #   ))
  # })
  # 
  # 
  # # UI - DIAGNOSTICS - 3 ------------------------------------------------------------------
  # 
  # output$box7 <- renderUI({
  #   div(
  #     style = "position: relative",
  #     tabBox(
  #       id = "box7",
  #       width = NULL,
  #       height = 400,
  #       tabPanel(
  #         title = "First isolates in selected diagnostics", 
  #         div(
  #           style = "position: absolute; left: 0.5em; bottom: 0.5em;",
  #           dropdown(
  #             radioGroupButtons(
  #               inputId = "box7.0",
  #               label = "Select group", 
  #               choiceNames = c("All", "Antimicrobial - Groups", "Antimicrobials", "Year", "Gender", "Specialty", "Subspecialty", "Origin"),
  #               choiceValues = c("fullname", "ab_group", "ab_type", "year", "gender", "specialty", "sub_specialty", "adm_route"), 
  #               selected = "fullname", 
  #               direction = "vertical"
  #             ),
  #             size = "xs",
  #             icon = icon("gear", class = "opt"), 
  #             up = TRUE
  #           )
  #         ),
  #         div(
  #           style = "position: absolute; left: 4em; bottom: 0.5em;",
  #           dropdown(
  #             sliderInput(
  #               inputId = "box7.1",
  #               label = "Show top ...", 
  #               min = 0, 
  #               max = 50, 
  #               value = c(25), 
  #               step = 5
  #             ),
  #             size = "xs",
  #             icon = icon("search-plus", class = "opt"), 
  #             up = TRUE
  #           )
  #         ),
  #         div(
  #           style = "position: absolute; left: 7.5em; bottom: 0.5em;",
  #           dropdown(
  #             downloadButton(outputId = "down_box_micro", label = "Download plot"),
  #             size = "xs",
  #             icon = icon("download", class = "opt"), 
  #             up = TRUE
  #           )
  #         ),
  #         withSpinner(
  #           plotOutput("micro_plot", height = 300),
  #           type = 4,
  #           color = "#d33724",
  #           size = 0.7
  #         )
  #       ),
  #       tabPanel(
  #         title = "First isolates - table",
  #         div(
  #           style = "position: absolute; left: 0.5em; bottom: 0.5em;",
  #           dropdown(
  #             radioGroupButtons(
  #               inputId = "box7.2",
  #               label = "Select group", 
  #               choiceNames = c("All", "Year", "Gender", "Specialty", "Subspecialty", "Origin"),
  #               choiceValues = c("fullname", "year", "gender", "specialty", "sub_specialty", "adm_route"), 
  #               selected = "fullname", 
  #               direction = "vertical"
  #             ),
  #             size = "xs",
  #             icon = icon("gear", class = "opt"), 
  #             up = TRUE
  #           )
  #         ),
  #         withSpinner(
  #           DT::dataTableOutput("micro_table", height = 300),
  #           type = 4,
  #           color = "#d33724",
  #           size = 0.7
  #         )
  #       ),
  #       div(
  #         style = "position:absolute;right:0.5em;bottom: 0.5em;",
  #         conditionalPanel(
  #           "input.box7 == 'First isolates in selected diagnostics'",
  #           actionBttn(
  #             inputId = "micro_plus",
  #             icon = icon("search-plus", class = "opt"),
  #             style = "fill",
  #             color = "danger",
  #             size = "xs"
  #           )
  #         )
  #       )
  #     )
  #   )
  # })
  # 
  # observeEvent((input$micro_plus), {
  #   showModal(modalDialog(
  #     renderPlot({
  #       micro_plot() + theme(
  #         axis.title = element_text(size = 20),
  #         text = element_text(size = 20),
  #         plot.title = element_text(size = 26)
  #       )
  #     }, height = 600),
  #     easyClose = TRUE,
  #     size = "l",
  #     footer = NULL
  #   ))
  # })
  # 
  # # UI - DIAGNOSTICS - 4 ------------------------------------------------------------------
  # 
  # output$box8 <- renderUI({
  #   div(
  #     style = "position: relative",
  #     div(
  #       style = "position: absolute; left: 4em; bottom: 0.5em;",
  #       dropdown(
  #         selectizeInput(
  #           inputId = "box8.1",
  #           label = "Select isolates",
  #           choices = sort(unique(microbiology$fullname)[!is.na(unique(microbiology$fullname))]),
  #           multiple = TRUE
  #         ),
  #         size = "xs",
  #         label = "Isolates",
  #         up = TRUE
  #       )
  #     ),
  #     div(
  #       style = "position: absolute; left: 9.5em; bottom: 0.5em;",
  #       dropdown(
  #         checkboxGroupButtons(
  #           width = "300px", 
  #           inputId = "box8.2",
  #           label = "Select antimicrobials",
  #           size = "xs", checkIcon = list("yes" = icon("check")),
  #           individual = TRUE,
  #           choiceValues = 
  #             sort(colnames(microbiology %>% select_if(is.rsi))),
  #           choiceNames = paste(
  #             ab_name(
  #               sort(colnames(microbiology %>% select_if(is.rsi)))),
  #             sep = ", "
  #           )
  #         ),
  #         size = "xs",
  #         label = "Antimicrobials",
  #         up = TRUE
  #       )
  #     ),
  #     tabBox(
  #       id = "box8",
  #       width = NULL,
  #       height = 400,
  #       tabPanel(
  #         title = "Resistance",
  #         div(
  #           style = "position: absolute; left: 0.5em; bottom: 0.5em;",
  #           dropdown(
  #             radioGroupButtons(
  #               inputId = "box8.0",
  #               label = "Select group",
  #               choiceNames = c("All", "Year", "Gender", "Specialty", "Subspecialty", "Origin"),
  #               choiceValues = c("fullname", "year", "gender", "specialty", "sub_specialty", "adm_route"),
  #               selected = "fullname",
  #               direction = "vertical"
  #             ),
  #             radioGroupButtons(
  #               inputId = "box8.0.1",
  #               label = "Select calculation",
  #               choices = c("Count", "Proportion"),
  #               selected = "Count",
  #               direction = "vertical"
  #             ),
  #             size = "xs",
  #             icon = icon("gear", class = "opt"),
  #             up = TRUE
  #           )
  #         ),
  #         div(
  #           style = "position: absolute; left: 17.5em; bottom: 0.5em;",
  #           dropdown(
  #             downloadButton(outputId = "down_box_res", label = "Download plot"),
  #             size = "xs",
  #             icon = icon("download", class = "opt"), 
  #             up = TRUE
  #           )
  #         ),
  #         withSpinner(
  #           plotOutput("isolate_plot", height = 300),
  #           type = 4,
  #           color = "#d33724",
  #           size = 0.7
  #         )
  #       ),
  #       tabPanel(
  #         title = "Resistance - over time",
  #         div(
  #           style = "position: absolute; left: 0.5em; bottom: 0.5em;",
  #           dropdown(
  #             radioGroupButtons(
  #               inputId = "box8.3",
  #               label = "Select time", 
  #               choiceNames = c("per year", "per month", "per quarter"),
  #               choiceValues = c("year", "yearmonth_test", "yearquarter_test"), 
  #               selected = "yearquarter_test", 
  #               direction = "vertical"
  #             ),
  #             size = "xs",
  #             icon = icon("gear", class = "opt"), 
  #             up = TRUE
  #           )
  #         ),
  #         div(
  #           style = "position: absolute; left: 17.5em; bottom: 0.5em;",
  #           dropdown(
  #             downloadButton(outputId = "down_box_res_ts", label = "Download plot"),
  #             size = "xs",
  #             icon = icon("download", class = "opt"), 
  #             up = TRUE
  #           )
  #         ),
  #         withSpinner(
  #           plotOutput("isolate_ts", height = 300),
  #           type = 4,
  #           color = "#d33724",
  #           size = 0.7
  #         )
  #       ),
  #       tabPanel(
  #         title = "Resistance - table",
  #         div(
  #           style = "position: absolute; left: 0.5em; bottom: 0.5em;",
  #           dropdown(
  #             radioGroupButtons(
  #               inputId = "box8.7",
  #               label = "Select group", 
  #               choiceNames = c("None", "Year", "Month", "Quarter", "Gender", "Specialty", "Subspecialty", "Origin"),
  #               choiceValues = c("fullname", "year", "yearmonth_test", "yearquarter_test", "gender", "specialty", "sub_specialty", "adm_route"), 
  #               selected = "fullname", 
  #               direction = "vertical"
  #             ),
  #             radioGroupButtons(
  #               inputId = "box8.8",
  #               label = "Select calculation", 
  #               choices = c("Count", "Proportion"),
  #               selected = "Count", 
  #               direction = "vertical"
  #             ),
  #             size = "xs",
  #             icon = icon("gear", class = "opt"), 
  #             up = TRUE
  #           )
  #         ),
  #         withSpinner(
  #           DT::dataTableOutput("isolate_table", height = 300),
  #           type = 4,
  #           color = "#d33724",
  #           size = 0.7
  #         )
  #       ),
  #       div(
  #         style = "position:absolute;right:0.5em;bottom: 0.5em;",
  #         conditionalPanel(
  #           "input.box8 == 'Resistance'",
  #           actionBttn(
  #             inputId = "res_plus",
  #             icon = icon("search-plus", class = "opt"),
  #             style = "fill",
  #             color = "danger",
  #             size = "xs"
  #           )
  #         )
  #       ),
  #       div(
  #         style = "position:absolute;right:0.5em;bottom: 0.5em;",
  #         conditionalPanel(
  #           "input.box8 == 'Resistance - over time'",
  #           actionBttn(
  #             inputId = "res_ts",
  #             icon = icon("search-plus", class = "opt"),
  #             style = "fill",
  #             color = "danger",
  #             size = "xs"
  #           )
  #         )
  #       )
  #     )
  #   )
  # })
  # 
  # observeEvent((input$res_plus), {
  #   showModal(modalDialog(
  #     renderPlot({
  #       isolate_plot() + theme(
  #         axis.title = element_text(size = 20),
  #         text = element_text(size = 20),
  #         plot.title = element_text(size = 26)
  #       )
  #     }, height = 600),
  #     easyClose = TRUE,
  #     size = "l",
  #     footer = NULL
  #   ))
  # })
  # 
  # observeEvent((input$res_ts), {
  #   showModal(modalDialog(
  #     renderPlot({
  #       isolate_ts() + theme(
  #         axis.title = element_text(size = 20),
  #         text = element_text(size = 20),
  #         plot.title = element_text(size = 26)
  #       )
  #     }, height = 600),
  #     easyClose = TRUE,
  #     size = "l",
  #     footer = NULL
  #   ))
  # })
}
  

  
