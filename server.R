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
    df <- as.data.frame(fread(input$inputFile$datapath, sep = "\t", header = T))
    data$df <- df
  })
  
  # UI - PATIENTS - 1 ----------------------------------------------------------
  
  output$box_pat <- renderUI({
      box(
        title = "Manhattan Plot",
        status = "primary",
        collapsible = FALSE,
        solidHeader = FALSE,
        width = 12,
        withSpinner(
          plotlyOutput("mb_marginal_manhattan_plot", height = 400)
        )
      )
  })
  
  observeEvent(input$gwas_marginal, {
    
    chrom_lengths <- chrom_lengths_hg38[extract_which_chr(data$df)]
    data$df <- add_cumulative_pos(data$df, chrom_lengths)
    data$df <- add_color(data$df, color1 = "black", color2 = "grey")

    y.max <- floor(max(data$df$P_Value_Marginal)) + 2

    x_breaks=get_x_breaks(chrom_lengths_hg38)
    names(x_breaks)[20]="20"
    names(x_breaks)[22]="22"
    color_map=unique(data$df$color)
    names(color_map)=unique(data$df$color)


    mhplot <- ggplot(data$df, aes(x=cumulative_pos, y=-log10(P_Value_Marginal), color=color)) +
                geom_point() +
                ggtitle("") +
                xlab("Chromosome") +
                ylab("-log<sub>10</sub>(<i>p</i>)") +
                scale_x_continuous(expand = c(0.01,0),
                                   breaks = x_breaks,
                                   labels = names(x_breaks)) +
                scale_y_continuous(expand = c(0.01,0)) +
                scale_color_manual(values = color_map,
                                   guide  = 'none') +
                theme(panel.background = element_blank(),
                      panel.grid       = element_line(color = "grey97"),
                      axis.line        = element_line(linewidth = 0.6),
                      axis.title       = element_text(size = 17, face = "bold"),
                      axis.text        = element_text(size = 12, face = "bold"),
                      legend.position = "none")
    output$mb_marginal_manhattan_plot <- renderPlotly({
      mhplot
    })
  })
  
  
  # # UI - PATIENTS - 2 -------------------------------------------------------
  # 
  # output$box_pat2 <- renderUI({
  #   div(
  #     style = "position: relative",
  #     tabBox(
  #       id = "box_pat2",
  #       width = NULL,
  #       height = 400,
  #       tabPanel(
  #         title = "Subspecialties - table",
  #         htmlOutput("patients_total"),
  #         withSpinner(
  #           DT::dataTableOutput("table_pat_all"),
  #           type = 4,
  #           color = "#d33724",
  #           size = 0.7
  #         )
  #       ),
  #       tabPanel(
  #         title = "Patient age",
  #         div(
  #           style = "position: absolute; left: 0.5em; bottom: 0.5em;",
  #           dropdown(
  #             radioGroupButtons(
  #               inputId = "box_pat1.1",
  #               label = "Select group", 
  #               choiceNames = c("All", "Gender"),
  #               choiceValues = c("all", "gender"), 
  #               selected = "all", 
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
  #             downloadButton(outputId = "down_age_select", label = "Download plot"),
  #             size = "xs",
  #             icon = icon("download", class = "opt"), 
  #             up = TRUE
  #           )
  #         ),
  #         withSpinner(
  #           plotOutput("plot_age_select", height = 300),
  #           type = 4,
  #           color = "#d33724",
  #           size = 0.7
  #         )
  #       )
  #     )
  #   )
  # })
  # 
  # output$patients_total <- renderText({
  #   HTML(
  #     paste("Total number of admissions:", 
  #           strong(
  #             unique(
  #               paste(set_reac_1()$id, set_reac_1()$adm_id) %>% 
  #                 length()
  #             )
  #           )
  #     )
  #   )
  # })
  # 
  # 
  # 
  # # UI - PATIENTS - 3 -------------------------------------------------------
  # 
  # output$box_year <- renderUI({
  #   div(
  #     style = "position: relative",
  #     tabBox(
  #       id = "box_year",
  #       width = NULL,
  #       height = 400,
  #       tabPanel(
  #         title = "Number of patients per year",
  #         div(
  #           style = "position: absolute; left:0.5em; bottom: 0.5em;",
  #           dropdown(
  #             radioGroupButtons(
  #               inputId = "box_year1",
  #               label = "Select time period", 
  #               choiceNames = c("Years", "Quarter", "Months"),
  #               choiceValues = c("years", "yearquarter_adm", "yearmonth_adm"), 
  #               selected = "years", 
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
  #             downloadButton(outputId = "down_year_select", label = "Download plot"),
  #             size = "xs",
  #             icon = icon("download", class = "opt"), 
  #             up = TRUE
  #           )
  #         ),
  #         withSpinner(
  #           plotOutput("plot_year_select", height = 300),
  #           type = 4,
  #           color = "#d33724",
  #           size = 0.7
  #         )
  #       )
  #     )
  #   )
  # })
  # 
  # 
  # 
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
  

  
