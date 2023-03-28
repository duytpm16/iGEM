library(modules)
CONSTS <- use("constants.R")

ui <- dashboardPage(
  # HEADER ------------------------------------------------------------------
  
  dashboardHeader(
    # App title visible in browser tab
    title = CONSTS$APP_TITLE,
    # App title visible
    tags$li(class = "dropdown title", tags$h1(CONSTS$APP_TITLE)),
    # App current version
    tags$li(class = "dropdown version", tags$p(CONSTS$APP_VERSION)),
    # App time range
    tags$li(class = "dropdown time-range", tags$p(CONSTS$APP_TIME_RANGE))
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
          bsButton("gwas", 
                   label = "GWAS RESULTS", 
                   icon = icon("chart-bar"), 
                   style = "success"),
          bsButton("diagnostics", 
                   label = "DIAGNOSTICS", 
                   icon = icon("flask", class = "flask-box"), 
                   style = "success")
      )
    ),
    
    fluid_design("diagnostics_panel", "box5", "box6", "box7", "box8"),
    
    fluidRow(
      div(
        id = "gwas_panel",
        fluidRow(
          column(
            width = 12,
            bsButton("gwas_marginal",
                     label = "Marginal",
                     style = "default"),
            bsButton("gwas_interaction",
                     label = "Interaction",
                     style = "default"),
            bsButton("gwas_joint",
                     label = "Joint",
                     style = "default"),
          ),
        ),
        br()
      )
    )
    ,
    fluidRow(
      div(
        id = "gwas_marginal_panel",
        column(
          width = 8,
          uiOutput("box_pat")
        ),
        column(
          width = 6,
          uiOutput("box_pat2")
        ),
        column(
          width = 6,
          uiOutput("box_year")
        )
      )
    )
  )
)
