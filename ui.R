options(shiny.maxRequestSize=6000*1024^2) 

ui <- dashboardPage(
  skin = "black",
  
  # HEADER ------------------------------------------------------------------
  
  dashboardHeader(
    # App title visible in browser tab
    # App title visible
    tags$li(class = "dropdown title", tags$h1("iGEM", style = "50px")),
    # App current version
    tags$li(class = "dropdown version", tags$p("1.0.0")),
    # App time range
    tags$li(class = "dropdown time-range", tags$p(""))
  ),
  
  # SIDEBAR -----------------------------------------------------------------
  
  dashboardSidebar(
    width = 300,
    sidebarMenu(
      fileInput("inputFile", h6("Input File:", style = "font-size: 20px; color: #000000"), 
                accept = c("text/plain", ".txt", ".out")),
      br(),
      h5("Plot Configurations:",
         style = "font-family: Source Sans Pro; padding-left:10px; font-size: 20px ;color: #000000;"),
      menuItem(
        h6("Manhattan", style = "font-family: Source Sans Pro; font-size: 20px; color: #000000"),
      
        fluidRow(
          column(
            width = 9,
            numericInput("mh_sigthreshold", label = "Significance Threshold", value = 1e-8, min = 0, max = 1),
          )
        ),
        fluidRow(
          column(
            width = 9,
            textInput("mh_sigcolor", label = "Significance Color", value = "red"),
          )
        ),
        fluidRow(
          column(
            width = 9,
            textInput("mh_chrcolor", label = "Chromosome Colors",  value = "black;darkgray"),
          )
        )
      )
    )
  ),
  
  
  # BODY --------------------------------------------------------------------
  
  dashboardBody(
    tags$head(
      tags$link(
        rel = "stylesheet", 
        type = "text/css", 
        href = "radar_style.css")
    ),
    
    useShinyjs(),
    
    # MAIN BODY ---------------------------------------------------------------
    
    fluidRow(
      column(
          width = 12,
          style = 'padding-left:0px; padding-right:0px; padding-top:10px; padding-bottom:0px',
          bsButton("gwas", 
                   label = "GWIS", 
                   icon  = icon("chart-column"),
                   style = "success")
      )
    ),

    tags$head(tags$style('.btn-default{ margin-left: -3px;border-radius: 5px; border: 2px solid black;}')),  # add the spacing,
    fluidRow(
        id = "gwas_panel",
        column(
            width = 12,
            style = 'padding-left:0px; padding-right:0px; padding-top:11px; padding-bottom:0px',
            selectInput(inputId = "gwis_choice",
                        label   = "Select test:",
                        choices = c("Joint"       = "joint",
                                    "Interaction" = "interaction",
                                    "Marginal"    = "marginal"))
        ),
        column(
            width = 12,
            style = 'padding-left:0px; padding-right:0px; padding-top:11px; padding-bottom:0px',
            selectInput(inputId = "se_choice",
                        label   = "Select standard errors:",
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
