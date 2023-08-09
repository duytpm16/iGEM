# iGEM

An R Shiny application to summarize and visualize gene-environment interaction (GEI) results from the [GEM](https://github.com/large-scale-gxe-methods/GEM), [REGEM](https://github.com/large-scale-gxe-methods/REGEM), and [METAGEM](https://github.com/large-scale-gxe-methods/METAGEM) software programs.

<br />

## 1. Requirements  

- bslib (>= <b>0.5.0</b>)
- data.table
- DT
- ggplot2
- kit
- plyr
- shiny
- shinyBS
- shinycssloaders
- shinyjs
<br />
- R (>= <b>3.5.0</b>)

<br />

## 2. Usage

1. **To run iGEM locally**

It is recommended to run iGEM locally for large file sizes (> 500 MB). 

```
shiny::runApp("iGEM")
```

2. **Using the web**

Alternatively, iGEM can be accessed at [shinyapps.io](https://posit.co/download/rstudio-desktop/).
> **Note** 
> 
> iGEM is deployed for free at [shinyapps.io](https://posit.co/download/rstudio-desktop/). It allows up to 1 GB of memory usage. Therefore, a large file will not be uploaded and/or summarized successfully. In this case, please run iGEM locally.

## 3. Overview  

**Figure 1.** Overview of the iGEM application

<img src="https://github.com/duytpm16/iGEM/tree/main/figures/igem_fig1.png" width = "75%"/>

