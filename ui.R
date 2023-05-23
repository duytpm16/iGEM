options(shiny.maxRequestSize=6000*1024^2) 

ui <- dashboardPage(
  # HEADER ------------------------------------------------------------------
  
  dashboardHeader(
    # App title visible in browser tab
    title = "iGEM",
    # App title visible
    tags$li(class = "dropdown title", tags$h1("iGEM")),
    # App current version
    tags$li(class = "dropdown version", tags$p("1.0.0")),
    # App time range
    tags$li(class = "dropdown time-range", tags$p(""))
  ),
  
  # SIDEBAR -----------------------------------------------------------------
  
  dashboardSidebar(
    width = 300,
    fileInput("inputFile", "INPUT FILE:", 
              accept = c("text/plain", ".txt", ".out"))
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
                   icon  = icon("chart-column"))
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
