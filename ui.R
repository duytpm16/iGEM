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
    
    # MAIN BODY ---------------------------------------------------------------
    
    fluidRow(
      column(
          width = 12,
          style = 'padding-left:28px; padding-right:0px; padding-top:10px; padding-bottom:0px',
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
    tags$head(tags$style('.btn-default{ margin-left: 15px;}')),  # add the spacing,
    fluidRow(
        id = "gwas_panel",
        fluidRow(
          column(
            width = 12,
            style = 'padding-left:28px; padding-right:0px; padding-top:10px; padding-bottom:0px',
            bsButton("gwas_marginal",
                     label = "MARGINAL",
                     style = "default"),
            bsButton("gwas_interaction",
                     label = "INTERACTION",
                     style = "default"),
            bsButton("gwas_joint",
                     label = "JOINT",
                     style = "default"),
          ),
        ),
        br()
    )
    ,
    fluidRow(
        id = "gwas_marginal_panel",
        column(
          width = 8,
          uiOutput("box_pat")
        ),
        column(
          width  = 4,
          offset = 0, 
          style  = 'padding-left:0px; padding-right:0px; padding-top:0px; padding-bottom:0px',
          uiOutput("box_pat2")
        )
    ),
    br(),
    fluidRow(
        id = "gwas_marginal_ss_panel",
        fluidRow(
          id = "gwas_marginal_tables",
          column(width  = 6,
                 offset = 0, 
                 style  = 'padding-left:30px; padding-right:0px; padding-top:0px; padding-bottom:0px',
                 uiOutput("box_pat3")),
          column(width = 6,
                 uiOutput("box_pat4")),          
        )
    ),
    br(),
    br(),
  )
)
