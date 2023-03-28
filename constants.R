import("shiny")

APP_TITLE <- "iGEM"
APP_TIME_RANGE <- ""
APP_VERSION <- "1.0.0"

lsgxe_website <- "https://github.com/large-scale-gxe-methods"
marketplace_website <- "https://appsilon.com/"
marketplace_name <- "Appsilon's Shiny Marketplace"

COLORS <- list(
  white = "#FFF",
  black = "#0a1e2b",
  primary = "#0099F9",
  secondary = "#15354A",
  ash = "#B3B8BA",
  ash_light = "#e3e7e9"
)


appsilon_footer <- tags$h3(
  class = "footer-heading",
  tags$span("This application is powered by a template from"),
  tags$a(
    class = "footer-link",
    href = marketplace_website,
    target = "_blank",
    rel = "nofollow noreferrer",
    marketplace_name
  )
)
