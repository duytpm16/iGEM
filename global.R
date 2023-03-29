# INSTALL DEPENDENCIES ----------------------------------------------------

source('dependencies.R')
# load all packages
lapply(required_packages, require, character.only = TRUE)

# https://genome.ucsc.edu/goldenpath/help/hg38.chrom.sizes
chrom_lengths_hg38=c("1"  = 248956422, "2"  = 242193529, "3"  = 198295559,
                     "4"  = 190214555, "5"  = 181538259, "6"  = 170805979,
                     "7"  = 159345973, "8"  = 145138636, "9"  = 138394717,
                     "10" = 133797422, "11" = 135086622, "12" = 133275309,
                     "13" = 114364328, "14" = 107043718, "15" = 101991189,
                     "16" = 90338345,  "17" = 83257441,  "18" = 80373285,
                     "19" = 58617616,  "20" = 64444167,  "21" = 46709983,
                     "22" = 50818468,   "X" = 156040895,  "Y" = 57227415)


# FLUID DESIGN FUNCTION ---------------------------------------------------

fluid_design <- function(id, test, model) {
  fluidRow(
    id = id,
    div(
      column(
        width = 8,
        uiOutput(paste0(test, "_", model, "_", "manhattan_box"))
      ),
      column(
        width  = 4,
        offset = 0, 
        style  = 'padding-left:0px; padding-right:0px; padding-top:0px; padding-bottom:0px',
        uiOutput(paste0(test, "_", model, "_", "qq_box"))
      ),
      column(width  = 6,
             offset = 0, 
             style  = 'padding-left:30px; padding-right:0px; padding-top:25px; padding-bottom:0px',
             uiOutput(paste0(test, "_", model, "_", "variants_table"))),
      column(width = 6,
             style  = 'padding-left:0px; padding-right:0px; padding-top:25px; padding-bottom:0px',
             uiOutput(paste0(test, "_", model, "_", "ss_tables")))
    )
  )
}


manhattan_box <- function(plotOutputId) {
  box(
    title = p("Manhattan Plot", style = 'font-size:21px;'),
    status = "primary",
    collapsible = FALSE,
    solidHeader = FALSE,
    width = 12,
    withSpinner(
      plotOutput(plotOutputId, height = 300,
                 click = paste0(plotOutputId, "_click"))
    )
  )
}


variant_box <- function(tableOutputId) {
  box(
    title = p("Variants Table", style = 'font-size:21px;'),
    status = "primary",
    collapsible = FALSE,
    solidHeader = FALSE,
    width = 12,
    withSpinner(
      dataTableOutput(tableOutputId, height = 350)
    )
  )
}

plot_manhattan <- function(df, x_breaks, color_map, test) {
  pcol <- paste0("P_Value_", test)
  y.max <- floor(max(df[,pcol])) + 5
  colnames(df)[which(colnames(df) == pcol)] <- "PV"
  
  ggplot(df, aes(x=cumulative_pos, y=PV, color=color)) +
    geom_point() +
    ggtitle("") +
    xlab("Chromosome") +
    ylab(expression(-log[10](italic(p)))) +
    scale_x_continuous(expand = c(0.01,0),
                       breaks = x_breaks,
                       labels = names(x_breaks)) +
    scale_y_continuous(expand = c(0.01,0), limits = c(0, y.max)) +
    scale_color_manual(values = color_map,
                       guide  = 'none') +
    theme(panel.background = element_blank(),
          panel.grid       = element_line(color = "grey97"),
          axis.line        = element_line(linewidth = 0.6),
          axis.title       = element_text(size = 13),
          axis.text        = element_text(size = 11),
          legend.position = "none")
}


variant_table <- function(df, variant_colnames, cat_interactions) {
  DT::datatable(
    df[, c("SNPID", "CHR", "POS", "Non_Effect_Allele", "Effect_Allele", "N_Samples", "AF", cat_interactions)],
    colnames = variant_colnames,
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
}

ss_box <- function(boxTitle, tableOutputPrefix) {
  tableOutputIds <- paste0(tableOutputPrefix, 2:4)
    
  box(
      title = p(boxTitle, style = 'font-size:21px;text-underline-position: under;text-decoration: underline;'),
      status = "primary",
      collapsible = FALSE,
      solidHeader = FALSE,
      width = 12,
      tags$div(
        class = "main-content-grid advanced-grid",
        dataTableOutput(tableOutputIds[1]),
        dataTableOutput(tableOutputIds[2]),
        dataTableOutput(tableOutputIds[3]),
      )
  )
}


mb_ss_tables <- function(output, test, df, row, beta_columns, se_columns, int_colnames, covariances, cov_rownames) {
  beta_se <- list()
  beta_se[["beta"]] <- df[row, beta_columns, drop = F]
  beta_se[["se"]] <- df[row, se_columns, drop = F]
  beta_se <- rbindlist(beta_se, use.names = FALSE)
  beta_se <- signif(beta_se, digits = 6)
  output[[paste0("mb_", test, "_manhattan_plot_table2")]] <- DT::renderDT({
    DT::datatable(
      beta_se,
      colnames = int_colnames,
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
  
  covs <- cbind(cov_rownames, c(df[row, covariances, drop = T]))
  output[[paste0("mb_", test, "_manhattan_plot_table3")]] <- DT::renderDT({
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
        scrollY = "250px",
        pageLength = 1000,
        columnDefs = list(list(width = '150px', targets = "_all", className = "dt-center"))
      )
    )
  })
  
  pvals <- df[row, c("P_Value_Marginal", "P_Value_Interaction", "P_Value_Joint")]
  pvals <- signif(10^-pvals, digits = 6)
  output[[paste0("mb_", test, "_manhattan_plot_table4")]] <- DT::renderDT({
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
}
