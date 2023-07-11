options(shiny.maxRequestSize=6000*1024^2) 

ui <- fluidPage(
  theme = bs_theme(bootswatch = "simplex"),
  
  navbarPage(
    
    title = p("iGEM", style = "font-size: 30px; font-weight: bold; padding-top: 50px"),
    theme = bs_theme(bootswatch = "simplex"),
    tabPanel(
      title = h4("GWIS", style = "padding-left: 20px; font-weight:bold;"),
      # SIDEBAR -----------------------------------------------------------------
      sidebarLayout(
        sidebarPanel(
          width = 2,
          tags$style(".progress-bar{color: white;}"),
          fileInput("inputFile", h5("Input File:", style = "font-weight: bold;"), 
                    accept = c("text/plain", ".txt", ".out")),
          br(),
          h5("Manhattan Plot:", style = "font-weight: bold;"),
          fluidRow(
            column(
              width = 9,
              numericInput("mh_sigThreshold", label = h6("Significance Threshold:", style = "font-weight: bold;"), value = 8, min = 0),
            )
          ),
          fluidRow(
            column(
              width = 9,
              textInput("mh_sigColor", label = h6("Significance Color:", style = "font-weight: bold;"), value = "red"),
            )
          ),
          fluidRow(
            column(
              width = 9,
              textInput("mh_chrColor", label = h6("Chromosome Colors:", style = "font-weight: bold;"),  value = "black;darkgray"),
            )
          )
        ),
        
        
        # BODY --------------------------------------------------------------------
        
        mainPanel(
          width = 10,
          
          useShinyjs(),
          
          tags$head(
            tags$link(
              rel = "stylesheet", 
              type = "text/css", 
              href = "style.css")
          ),
          tags$style(HTML('table.dataTable tr.active td, table.dataTable tr.active {background-color: red !important;}')),
          # MAIN BODY ---------------------------------------------------------------
          
          fluidRow(
            column(
              width = 12,
              style = 'padding-left:0px; padding-right:0px; padding-top:10px; padding-bottom:0px',
              bsButton("gwas", 
                       label = "GWIS", 
                       icon  = icon("chart-column"),
                       style = "secondary")
            )
          ),
          
          
          fluidRow(
            id = "gwas_panel",
            column(
              width = 2,
              style = 'padding-left:0px; padding-right:0px; padding-top:11px; padding-bottom:0px',
              selectInput(inputId = "gwis_choice",
                          label   = h5("Select a test:"),
                          choices = c("Joint"       = "joint",
                                      "Interaction" = "interaction",
                                      "Marginal"    = "marginal"))
            ),
            column(
              width = 6,
              style = 'padding-left:10px; padding-right:0px; padding-top:11px; padding-bottom:0px',
              selectInput(inputId = "se_choice",
                          label   = h5("Select standard errors:"),
                          choices = c("Model-based" = "modelbased",
                                      "Robust"      = "robust"))
            ),
            br(),
          ),
          
          hidden(fluid_design("mb", "marginal")),
          hidden(fluid_design("rb", "marginal")),
          hidden(fluid_design("mb", "interaction")),
          hidden(fluid_design("rb", "interaction")),
          hidden(fluid_design("mb", "joint")),
          hidden(fluid_design("rb", "joint")),
          br(),
          br(),
          br(),
        )
      )
    )
  )
)