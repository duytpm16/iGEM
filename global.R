# INSTALL DEPENDENCIES ---------------------------------------------------------

source('dependencies.R')


# FLUID DESIGN FUNCTION --------------------------------------------------------

fluid_design <- function(test, model) {
  fluidRow(
    id = paste0("gwis_", test, "_", model, "_panel"),
    div(
      fluidRow(
        column(
          width  = 8,
          offset = 0,
          style  = 'padding-left:0px; padding-right:0px; padding-top:0px; padding-bottom:0px',
          uiOutput(paste0(test, "_", model, "_", "manhattan_box"))
        ),
        column(
          width  = 4,
          offset = 0, 
          style  = 'padding-left:15px; padding-right:0px; padding-top:0px; padding-bottom:0px',
          uiOutput(paste0(test, "_", model, "_", "qq_box"))
        )
      ),
      fluidRow(
        column(width  = 6,
               offset = 0, 
               style  = 'padding-left:0px; padding-right:0px; padding-top:0px; padding-bottom:0px;',
               uiOutput(paste0(test, "_", model, "_", "variants_table"))
        ),
        column(width = 6,
               style  = 'padding-left:15px; padding-right:0px; padding-top:0px; padding-bottom:0px;',
               uiOutput(paste0(test, "_", model, "_", "ssTables"))
        )
      )
    )
  )
}


manhattan_box <- function(plotOutputId) {
  bslib::card(
    class = "card bg-secondary mb-3",
    style = "box-shadow: 5px 10px #D3D3D3; font-weight: bold;",
    card_header(
      style = "font-size: 20px;",
      "Manhattan Plot"
    ),
    card_body(
      div(
        style = "position: relative; overflow:hidden;",
        plotOutput(plotOutputId, 
                   height = 315,
                   click  = paste0(plotOutputId, "_click"),
                   hover  = hoverOpts(paste0(plotOutputId,"_hover"), delay = 100, delayType = "debounce")),
        uiOutput(paste0(plotOutputId, "_hover_info"))
      )
    )
  )
}

manhattan_plot <- function(df, x_breaks, sig_threshold, sig_color, chr_color) {
  if (is.null(df)) {
    return(NULL)
  }

  
  y.max <- ceiling(max(df$LOGP)) + 5
  
  ggplot(df, aes(x=POS, y=LOGP)) +
    geom_point(color = chr_color, size = 2.5, alpha = 0.5) +
    geom_hline(yintercept = sig_threshold, color = sig_color, linetype = "dashed") +
    ggtitle("") +
    xlab("Chromosome") +
    ylab(expression(-log[10](italic(p)))) +
    scale_x_continuous(expand = c(0.01,0), breaks = x_breaks, labels = names(x_breaks)) +
    scale_y_continuous(expand = c(0.01,0), limits = c(0, y.max), breaks = round(seq(0, y.max, length.out = 6))) +
    theme(panel.background = element_blank(),
          panel.border     = element_rect(colour = "black", fill=NA, linewidth =1.5),
          panel.grid       = element_line(color = "grey97"),
          axis.line        = element_blank(),
          axis.title       = element_text(size = 16),
          axis.title.x     = element_text(margin =  margin(t = 10, r = , b = 0, l = 0)),
          axis.title.y     = element_text(margin =  margin(t = 0, r = 10, b = 0, l = 0)),
          axis.text        = element_text(size = 14, face = "bold"),
          legend.position  = "none")
}

manhattan_tooltip <- function (hover, df) {
  if (is.null(df)) {
    return(NULL)
  }
  
  point <- nearPoints(df, hover, maxpoints =  1,
                      xvar = "POS",  yvar = "LOGP")
  
  if (is.null(point) || nrow(point) == 0) return(NULL)
  
  style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                  "left:", hover$coords_css$x + 2, "px; top:", hover$coords_css$y + 2, "px;")
  wellPanel(
    style = style,
    p(HTML(paste0("<b> CHR: </b>", point$CHR, "<br/>",
                  "<b> POS: </b>", point$POS, "<br/>",
                  "<b> -log10(p): </b>", format(round(point$LOGP, 2), nsmall = 2), "<br/>"))
    )
  )
}


qq_box <- function(plotOutputId) {
  bslib::card(
    class = "card bg-secondary mb-3",
    style = "box-shadow: 5px 10px #D3D3D3; font-weight: bold;",
    card_header(
      style = "font-size: 20px;",
      "Quantile-Quantile Plot"
    ),
    card_body(
      height = "348px",
      plotOutput(plotOutputId, height = 315)
    )
  )
}

qq_plot <- function(df, h) {
  if (is.null(df)) {
    return(NULL)
  }

  ggplot(df, aes(x = y, y = x)) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    geom_point() +
    ggtitle(bquote(lambda==.(h))) +
    theme(panel.background = element_blank(),
          panel.border     = element_rect(colour = "black", fill=NA, linewidth =1.5),
          panel.grid       = element_line(color = "grey97"),
          plot.title       = element_text(size = 18, face = "bold", hjust = 0.1, vjust = -15),
          axis.line        = element_line(linewidth = 0.6),
          axis.title       = element_text(size = 16, face = "bold"),
          axis.text        = element_text(size = 15, face = "bold")) +
    ylab(expression(paste('Observed ', -log[10](italic(p))))) +
    xlab(expression(paste('Expected ', -log[10](italic(p)))))
}


variantTable_box <- function(tableOutputId) {
  bslib::card(
    class = "card bg-secondary mb-3",
    style = "box-shadow: 5px 10px #D3D3D3; font-weight: bold;",
    card_header(
      style = "font-size: 20px",
      "Variants in Manhattan Plot Region"
    ),
    card_body(
      style = "height: 435px",
      dataTableOutput(tableOutputId, height = 400)
    )
  )
}

variantTable <- function(df, pcol, variant_colnames, cat_interactions) {
  DT::datatable(
    df[, c(pcol, "SNPID", "CHR", "POS", "Non_Effect_Allele", "Effect_Allele", "N_Samples", "AF", cat_interactions)],
    colnames  = variant_colnames,
    rownames  = FALSE,
    escape    = FALSE,
    selection = 'single',
    caption   = htmltools::tags$caption(
                  style = 'caption-side: top; text-align: left;',
                  "Select a row to view the summary statistics."
                ),
    options   = list(
      dom = 'tp',
      search = list(regex = TRUE, caseInsensitive = TRUE),
      pageLength = 5,
      ordering  = TRUE,
      stateSave = TRUE,
      autoWidth = TRUE,
      scrollX   = TRUE,
      columnDefs = list(list(targets = "_all", className = "dt-center", width = "75px"))
    )
  ) %>% formatRound(columns=c(1), digits=2)
}

ssTable_box <- function(tableOutputPrefix) {
  tableOutputIds <- paste0(tableOutputPrefix, 1:3)
  bslib::card(
    class = "card bg-secondary mb-3",
    style = "box-shadow: 5px 10px #D3D3D3; font-weight: bold;",
    card_header(
      style = "font-size: 20px;",
      "Summary Statistics"
    ),
    card_body(
      style = "height: 435px; overflow: hidden",
      uiOutput(paste0(tableOutputPrefix, "_title")),
      div(
        class = "main-content-grid advanced-grid",
        dataTableOutput(tableOutputIds[1]),
        dataTableOutput(tableOutputIds[2]),
        dataTableOutput(tableOutputIds[3]),
      )
    )
  )
}

ssTables <- function(output, se, test, df, int_colnames, beta_columns, se_columns, covariances, cov_rownames) {
  if (se == "mb") {
    ss_caption1 <- "Table 1: Coefficient Estimates and Model-based Standard Errors."
    ss_caption2 <- "Table 2: Model-based Covariances."
    ss_caption3 <- "Table 3: Model-based P-Values."
    ss_colname2 <- c("", "Covariances")
    ss_rowname1 <- c("Coefficients", "SE")
    ss_rowname3 <- c("P-Value")
    pcols <- c("P_Value_Marginal", "P_Value_Interaction", "P_Value_Joint")
  } else {
    ss_caption1 <- "Table 1: Coefficient Estimates and Robust Standard Errors."
    ss_caption2 <- "Table 2: Robust Covariances."
    ss_caption3 <- "Table 3: Robust P-Values."
    ss_colname2 <- c("", "Covariances<sub>R</sub>")
    ss_rowname1 <- c("Coefficients", "SE<sub>R</sub>")
    ss_rowname3 <- c("P-Value<sub>R</sub>")
    pcols <- c("robust_P_Value_Marginal", "robust_P_Value_Interaction", "robust_P_Value_Joint")
  }
  
  
  if (!is.null(df) || nrow(df) != 0) {
    beta_se <- list(beta = df[, beta_columns], se = df[, se_columns])
    beta_se <- rbindlist(beta_se, use.names = FALSE)
    beta_se <- apply(beta_se, 2, FUN = function(x) {formatC(x, format = "e", digits = 2)})
    
    covs <- df[, covariances, drop = FALSE]
    covs <- formatC(unlist(covs), format = "e", digits = 2)

    pvals <- df[, pcols, drop = FALSE]
    pvals <- apply(pvals, 1, FUN = function(x) {formatC(10^-x, format = "e", digits = 2)})
  }
  
  output[[paste0(se, "_", test, "_ssTable_title")]] <- renderText({
    req(df)
    HTML(paste0("<h3>", df$SNPID, "</h3>"))
  })
  
  output[[paste0(se, "_", test, "_ssTable1")]] <- DT::renderDT({
    req(beta_se)
    DT::datatable(
      beta_se,
      colnames  = int_colnames,
      rownames  = ss_rowname1,
      escape    = FALSE, 
      selection = 'none',
      caption   = htmltools::tags$caption(
                    style = 'caption-side: top; text-align: left;',
                    ss_caption1
                  ),
      options   = list(
                    autoWidth  = TRUE,
                    dom        = 't',
                    scrollX    = TRUE,
                    pageLength = 2,
                    columnDefs = list(list(targets = "_all", className = "dt-center", width = "75px"))
                  )
      )
  })
  
  output[[paste0(se, "_", test, "_ssTable2")]] <- DT::renderDT({
    req(covs)
    DT::datatable(
      as.data.frame(covs),
      colnames  = ss_colname2,
      rownames  = cov_rownames,
      escape    = FALSE,
      selection = 'none',
      caption   = htmltools::tags$caption(
                    style = 'caption-side: top; text-align: left;',
                    ss_caption2
                  ),
      options   = list(
                    dom        = 't',
                    pageLength = 1000,
                    ordering   = TRUE,
                    stateSave  = TRUE,
                    autoWidth  = TRUE,
                    scrollX    = TRUE,
                    scrollY    = "250px",
                    columnDefs = list(list(targets = "_all", className = "dt-center"))
                  )
      )
  })


  output[[paste0(se, "_", test, "_ssTable3")]] <- DT::renderDT({
    req(pvals)
    DT::datatable(
      t(pvals),
      rownames  = ss_rowname3,
      colnames  = c("Marginal", "Interaction", "Joint"),
      escape    = FALSE,
      selection = 'none',
      caption   = htmltools::tags$caption(
                    style = 'caption-side: top; text-align: left;',
                    ss_caption3
                  ),
      options   = list(
                    autoWidth  = TRUE,
                    dom        = 't',
                    scrollX    = TRUE,
                    pageLength = 1,
                    columnDefs = list(list(targets = "_all", className = "dt-center"))
                  )
      )
  })
}




# https://genome.ucsc.edu/goldenpath/help/hg38.chrom.sizes
chrom_lengths_hg38=c("1"  = 248956422, "2"  = 242193529, "3"  = 198295559,
                     "4"  = 190214555, "5"  = 181538259, "6"  = 170805979,
                     "7"  = 159345973, "8"  = 145138636, "9"  = 138394717,
                     "10" = 133797422, "11" = 135086622, "12" = 133275309,
                     "13" = 114364328, "14" = 107043718, "15" = 101991189,
                     "16" = 90338345,  "17" = 83257441,  "18" = 80373285,
                     "19" = 58617616,  "20" = 64444167,  "21" = 46709983,
                     "22" = 50818468,   "X" = 156040895,  "Y" = 57227415)



get_cumulative_length <- function(chrom_lengths) {
  cumulative_length <- 0
  
  if(length(chrom_lengths) > 1){
    cumulative_length <- head(c(0,cumsum(x = unname(chrom_lengths))), -1)
    names(cumulative_length) <- names(chrom_lengths)
  }
  return(cumulative_length)
}

get_x_breaks <- function(chrom_lengths) {
  cumulative_length <- get_cumulative_length(chrom_lengths)
  x_breaks <-cumulative_length+round(chrom_lengths/2)
  names(x_breaks)=gsub('chr', '', names(x_breaks))
  if(length(chrom_lengths) == 21){
    names(x_breaks)[20]='20'
  }
  if(length(chrom_lengths) > 21){
    names(x_breaks)[20]='20'
    names(x_breaks)[22]='22'
  }
  return(x_breaks)
}


get_chr_colors <- function(df, colors) {
  nchr <- nrow(df)
  cols <- rep(colors, ceiling(nchr / length(colors)))
  
  return(unlist(lapply(1:nchr, function(x) data.frame(color = rep(cols[x], df$n[x])))))
}


subset_data <- function(subDF, pcol, nvar) {
  colnames(subDF)[colnames(subDF) == pcol]  <- "LOGP"
  colnames(subDF)[colnames(subDF) == "cumulative_pos"] <- "POS"
  
  subDF[, index := 1:nvar]
  if (nvar > 100) {
    subDF[, round_pcol := round(LOGP, digits = 3)]
    subDF[, round_pos  := plyr::round_any(POS, 100000)]
    subDF[, duplicated := (fduplicated(subDF$round_pos) & fduplicated(subDF$round_pcol))]
    subDF <- subDF[!subDF$duplicated | subDF$LOGP > 8, ]
  }
  
  return(as.data.frame(subDF[,c("index", "CHR", "POS", "LOGP", "duplicated")]))
}
