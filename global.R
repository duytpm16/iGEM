# INSTALL DEPENDENCIES ---------------------------------------------------------

source('dependencies.R')

# https://genome.ucsc.edu/goldenpath/help/hg38.chrom.sizes
chrom_lengths_hg38=c("1"  = 248956422, "2"  = 242193529, "3"  = 198295559,
                     "4"  = 190214555, "5"  = 181538259, "6"  = 170805979,
                     "7"  = 159345973, "8"  = 145138636, "9"  = 138394717,
                     "10" = 133797422, "11" = 135086622, "12" = 133275309,
                     "13" = 114364328, "14" = 107043718, "15" = 101991189,
                     "16" = 90338345,  "17" = 83257441,  "18" = 80373285,
                     "19" = 58617616,  "20" = 64444167,  "21" = 46709983,
                     "22" = 50818468,   "X" = 156040895,  "Y" = 57227415)


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
          style  = 'padding-left:0px; padding-right:0px; padding-top:0px; padding-bottom:0px',
          uiOutput(paste0(test, "_", model, "_", "qq_box"))
        )
      ),
      fluidRow(
        column(width  = 6,
               offset = 0, 
               style  = 'padding-left:0px; padding-right:0px; padding-top:25px; padding-bottom:0px',
               uiOutput(paste0(test, "_", model, "_", "variants_table"))),
        
        column(width = 6,
               style  = 'padding-left:0px; padding-right:0px; padding-top:25px; padding-bottom:0px',
               uiOutput(paste0(test, "_", model, "_", "ss_tables")))
      )
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
      plotOutput(plotOutputId, 
                   height = 300,
                   click  = paste0(plotOutputId, "_click"))
    )
  )
}


qq_box <- function(plotOutputId) {
  box(
    title = p("Quantile-Quantile Plot", style = 'font-size:21px;'),
    status = "primary",
    collapsible = FALSE,
    solidHeader = FALSE,
    width = 12,
    withSpinner(
      plotOutput(plotOutputId, 
                 height = 300)
    )
  )
}


variantTable_box <- function(tableOutputId) {
  box(
    title = p("Variants Table", style = 'font-size:21px;'),
    status = "primary",
    collapsible = FALSE,
    solidHeader = FALSE,
    width = 12,
    withSpinner(
      dataTableOutput(tableOutputId, 
                      height = 400)
    )
  )
}

ssTable_box <- function(boxTitle, tableOutputPrefix) {
  tableOutputIds <- paste0(tableOutputPrefix, 1:3)
  
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


plot_manhattan <- function(df, x_breaks, sig_threshold, sig_color, chr_color, test, robust) {
  if (is.null(df)) {
    return(NULL)
  }
  
  pcol <- ifelse(robust, paste0("robust_", paste0("P_Value_", test)), paste0("P_Value_", test))
  y.max <- floor(max(df[,pcol])) + 5
  colnames(df)[which(colnames(df) == pcol)] <- "PV"
  
  ggplot(df, aes(x=cumulative_pos, y=PV)) +
    geom_hline(yintercept = -log10(sig_threshold), color = sig_color, alpha = 0.5, linetype = "dashed") +
    geom_point(color = chr_color$color) +
    ggtitle("") +
    xlab("Chromosome") +
    ylab(expression(-log[10](italic(p)))) +
    scale_x_continuous(expand = c(0.01,0),
                       breaks = x_breaks,
                       labels = names(x_breaks)) +
    scale_y_continuous(expand = c(0.01,0), limits = c(0, y.max)) +
    theme(panel.background = element_blank(),
          panel.grid       = element_line(color = "grey97"),
          axis.line        = element_line(linewidth = 0.6),
          axis.title       = element_text(size = 13),
          axis.text        = element_text(size = 11),
          legend.position = "none")
}


plot_qq <- function(df, pcol) {
  if (is.null(df)) {
    return(NULL)
  }
  
  qq_df <- drop_dense(sort(df[[pcol]], decreasing = T), -log10(stats::ppoints(length(df))))
  ggplot(qq_df, aes(x = y, y = x)) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    geom_point() +
    theme(panel.background = element_blank(),
          panel.grid       = element_line(color = "grey97"),
          axis.line        = element_line(linewidth = 0.6)) +
    ylab(expression(paste('Observed ', -log[10](italic(p))))) +
    xlab(expression(paste('Expected ', -log[10](italic(p)))))
}


variant_table <- function(df, pcol, variant_colnames, cat_interactions) {
  DT::datatable(
    df[, c(pcol, "SNPID", "CHR", "POS", "Non_Effect_Allele", "Effect_Allele", "N_Samples", "AF", cat_interactions)],
    colnames  = variant_colnames,
    rownames  = FALSE,
    escape    = FALSE,
    style     = "bootstrap",
    selection = 'single',
    caption   = "Select a row to view the summary statistics",
    options   = list(
                  dom = 'tp',
                  search = list(regex = TRUE, caseInsensitive = TRUE),
                  pageLength = 5,
                  ordering = TRUE,
                  stateSave = TRUE,
                  columnDefs = list(list(targets = "_all", className = "dt-center"))
                )
    ) %>% formatRound(columns=c(1), digits=2)
}


ss_tables <- function(output, se, test, df, row, int_colnames, beta_columns, se_columns, covariances, cov_rownames) {
  if (se == "mb") {
    ss_caption1 <- "Table 1: Coefficient Estimates and Model-based Standard Errors."
    ss_caption2 <- "Table 2: Model-based Covariances."
    ss_caption3 <- "Table 3: Model-based P-Values."
    ss_colname2 <- c("", "Covariances")
    ss_rowname1 <- c("Coefficients", "SE")
    ss_rowname3 <- c("P-Values")
    pcols <- c("P_Value_Marginal", "P_Value_Interaction", "P_Value_Joint")
  } else {
    ss_caption1 <- "Table 1: Coefficient Estimates and Robust Standard Errors."
    ss_caption2 <- "Table 2: Robust Covariances."
    ss_caption3 <- "Table 3: Robust P-Values."
    ss_colname2 <- c("", "Covariances<sub>R</sub>")
    ss_rowname1 <- c("Coefficients", "SE<sub>R</sub>")
    ss_rowname3 <- c("P-Values<sub>R</sub>")
    pcols <- c("robust_P_Value_Marginal", "robust_P_Value_Interaction", "robust_P_Value_Joint")
  }
  
  
  beta_se <- list()
  beta_se[["beta"]] <- df[row, beta_columns, drop = F]
  beta_se[["se"]]   <- df[row, se_columns,   drop = F]
  beta_se <- rbindlist(beta_se, use.names = FALSE)
  beta_se <- signif(beta_se, digits = 6)
  output[[paste0(se, "_", test, "_ss_table1")]] <- DT::renderDT({
    DT::datatable(
      beta_se,
      colnames  = int_colnames,
      rownames  = ss_rowname1,
      escape    = FALSE, 
      style     = "bootstrap",
      selection = 'none',
      caption   = ss_caption1,
      options   = list(
                    dom        = 't',
                    scrollX    = TRUE,
                    pageLength = 2,
                    columnDefs = list(list(width = '300px', targets = "_all", className = "dt-center"))
                  ))
  })
  
  
  covs <- cbind(cov_rownames, c(df[row, covariances, drop = T]))
  output[[paste0(se, "_", test, "_ss_table2")]] <- DT::renderDT({
    DT::datatable(
      covs,
      colnames  = ss_colname2,
      rownames  = FALSE,
      escape    = FALSE,
      style     = "bootstrap",
      selection = 'none',
      caption   = ss_caption2,
      options   = list(
                    dom        = 't',
                    scrollX    = TRUE,
                    scrollY    = "250px",
                    pageLength = 1000,
                    columnDefs = list(list(className = "dt-right",  targets = 0, width = '125px'),
                                      list(className = "dt-center", targets = 1, width = '100px'))
                  ))
  })

  
  pvals <- df[row, pcols]
  pvals <- signif(10^-pvals, digits = 6)
  output[[paste0(se, "_", test, "_ss_table3")]] <- DT::renderDT({
    DT::datatable(
      pvals,
      rownames  = ss_rowname3,
      colnames  = c("Marginal", "Interaction", "Joint"),
      escape    = FALSE,
      style     = "bootstrap",
      selection = 'none',
      caption   = ss_caption3,
      options   = list(
                    dom        = 't',
                    scrollX    = TRUE,
                    pageLength = 1,
                    columnDefs = list(list(width = '150px', targets = "_all", className = "dt-center"))
                  ))
  })
}




reduce_data <- function(ms, pcol) {
  colnames(ms) <- c("CHR", "LOGP")
  numc <- length(sort(unique(ms$CHR)));
  
  # part 3: reduce size for fast plotting --------------------------------------------------------------------------------------------------------------------------------------------
  quants <- c(0,0.002,0.5,0.998,1); 
  quants <- quantile(ms$LOGP, quants); # minimum, 0.2th percentile, median, 99.8th percentile and maximum are being calculated
  right  <- (quants[5] - quants[4])/(quants[4] - quants[3]); # measure of significance of right tail
  left   <- (quants[1] - quants[2])/(quants[2] - quants[3]); # measure of significance of left tail
  if (nrow(ms) < 1E5) { #if there are lest than 100k rows then full data is rounded to 3 digits
    digs <- 3; 
    ms$LOGP <- round(ms$LOGP, digits=digs); 
    ms      <-ms[!duplicated(ms$LOGP),]
  }
  else { # round lower and upper parts separately
    if (right>0.1) {  # significant right tail
      if (left>0.1) { # significant left tail
        f1 <- ms$LOGP <= quants[4] & ms$LOGP >= quants[2];
        digs1 <- 2; 
        digs2 <- 3; 
        digs  <- 0*f1 + digs2; 
        digs[f1] <- digs1;
        ms$LOGP  <- round(ms$LOGP, digits=digs); 
        rm(digs);
        f=NULL; #store vector of duplicated rows
        for (i in 1:numc) {
          m1 <- ms[ms$C==i, c("LOGP"), drop = FALSE]; # subset by chromosome for smaller memory
          f1 <- duplicated(m1$LOGP);
          f  <- c(f,f1);
          rm(f1, m1);
        }
        ms <- ms[!f,]; 
        rm(f);
      }
      else { # insignificant left tail
        f1 <- ms$LOGP <= quants[4];
        digs1 <- 2; 
        digs2 <- 3; 
        digs  <- 0*f1 + digs2; 
        digs[f1] <- digs1;
        ms$LOGP  <- round(ms$LOGP, digits=digs); 
        rm(digs);
        f=NULL; #store vector of duplicated rows
        for (i in 1:numc) {
          m1 <- ms[ms$C==i, c("LOGP"), drop = FALSE]; # subset by chromosome for smaller memory
          f1 <- duplicated(m1$LOGP);
          f  <- c(f,f1);
          rm(f1, m1);
        }
        ms=ms[!f, ]; 
        rm(f);
      }
    }
    else { # insignificant right tail
      if (left>0.1) { # significant left tail
        f1 <- ms$LOGP >= quants[2];
        digs1 <- 2; 
        digs2 <- 3; 
        digs  <- 0*f1 + digs2; 
        digs[f1] <- digs1;
        ms$LOGP  <- round(ms$LOGP, digits=digs); 
        rm(digs);
        f=NULL; #store vector of duplicated rows
        for (i in 1:numc) {
          m1 <- ms[ms$C==i,c("LOGP"), drop = FALSE]; # subset by chromosome for smaller memory
          f1 <- duplicated(m1$LOGP);
          f  <- c(f,f1);
          rm(f1, m1);
        }
        ms=ms[!f,]; 
        rm(f);
      }
      else { # insignificant left tail
        digs <- 3;
        ms$LOGP <- round(ms$LOGP, digits=digs);
        ms <- ms[!duplicated(ms$LOGP),]; 
      }
    }
  }
  
  return(ms)
}


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


