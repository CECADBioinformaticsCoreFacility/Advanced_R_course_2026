#' ---
#' title: "R Advanced Course 2026"
#' subtitle: "Differential Gene Expression"
#' date: "2026-11-24"
#' author: 
#'   - name: "Ulrike Göbel"
#'     affiliations: "Bioinformatics Core Facility CECAD"
#' title-slide-attributes:
#'   data-background-color: "#035AA6"
#' format:
#'   revealjs: #pptx beamer revealjs
#'     width: 1600
#'     height: 900
#'     preview-links: auto
#'     theme: [default, custom.scss]
#'     logo: images/CECAD-Logo-24-white.svg
#'     slide-number: true
#'     embed-resources: false #true
#'     highlight: zenburn
#'     smaller: false
#'     # chalkboard: 
#'     #   buttons: false
#' include-in-header: 
#'   text: |
#'     <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
#' ---
#' 
#' ## Slides & Code {.smaller}
#' 
#' :::: {.columns}
#' 
#' ::: {.column width=30%}
## ----QRcode, echo=FALSE, fig.width=2, fig.height=2, fig.align='left'----------
#fig.width=2, fig.height=2,
plot(qrcode::qr_code("https://cecadbioinformaticscorefacility.github.io/Intermediate_R_Course_2026/"))

#' 
#' - *[f]* Full screen
#' - *[o]* Slide Overview
#' - *[c]* Notes
#' - *[h]* help
#' 
#' :::
#' 
#' ::: {.column width=70%}
#' 
#' ### `r fontawesome::fa("github")` git repo 
#' 
#' [R-Basic](https://github.com/CECADBioinformaticsCoreFacility/Intermediate_R_Course_2026)
#' 
#' <br>
#' 
#' ### Clone repo
#' 
#' `git clone https://github.com/CECADBioinformaticsCoreFacility/Intermediate_R_Course_2026.git`
#' 
#' <br>
#' 
#' ### Slides Directly
#' 
#' [https://cecadbioinformaticscorefacility.github.io/Advanced_R_Course_2026/](https://cecadbioinformaticscorefacility.github.io/Advanced_R_Course_2026/)
#' 
#' :::
#' ::::
#' 
#' ::: {.notes}
#' Speaker notes go here.
#' :::
#' 
#' # {.center background-color="#03A696"}
#' 
#' ::: title-head
#' **Session 4:: Differential Gene Expression -- Models**
#' :::
#' 
## ----setup--------------------------------------------------------------------
#| echo: false
#| eval: true

## Silent Preliminaries

.libPaths("/home/ugoebel/R-4.4.3/lib/R/library")

library(dplyr)
library(ggplot2)
library(pheatmap)
library(DESeq2)
library(org.Mm.eg.db)
source("R/functions.R")
source("R/plotPCA2.R")

library(Cairo)
options(bitmapType="cairo")


#' ##
#' 
## ----load_from_Session2-------------------------------------------------------
#| echo: false
#| eval: true

colData <- readRDS(file="colData.rds")


#' 
#' ![](images/intro_4.png)
#' 
#' ## Understanding Model Formulae 
#' \         
#' 
#' Model formulae like  __`~ 1 + genotype`__ provide  a **high level language**
#' 
#' - for the **experimenter** to __`conceptually express a hypothesis`__ about his or her data
#' - for **modeling software** to know the __`structure of the tests to compute`__
#' \        
#' 
#' - **Independent of the actual mathematics of the tests !**
#' 
#' ## Understanding Model Formulae 
#' 
#' ![](images/model_page1.png)
#' 
#' ## Understanding Model Formulae 
#' 
#' ![](images/model_page2.png)
#' 
#' ## Understanding Model Formulae
#' 
#' ![](images/model_page3.png)
#' 
#' ## Understanding Model Formulae
#' 
#' ![](images/model_page4.png)
#' 
#' ## Understanding Model Formulae
#' 
#' ![](images/model_page5.png)
#' 
#' ## Understanding Model Formulae
#' 
#' ![](images/model_page6.png)
#' 
#' ## Understanding Model Formulae
#' 
#' ![](images/model_page7.png)
#' 
#' ## Understanding Model Formulae
#' 
#' ![](images/model_page8.png)
#' 
#' ## Understanding Model Formulae 
#' 
#' ![](images/model_page9.png)
#' 
#' ## Understanding Model Formulae -- Multiple Factors
#' 
#' ![](images/model_page10.png)
#' 
#' ## Understanding Model Formulae -- Multiple Factors
#' 
#' ![](images/model_page11_NEW3.png)
#' 
#' ## Understanding Model Formulae -- Multiple Factors
#' 
#' 
#' 
#' ![](images/model_page12.png)
#' 
#' 
#' ## Understanding Model Formulae -- Multiple Factors
#' 
#' ![](images/model_page13.png)
#' 
#' ## Understanding Model Formulae -- Multiple Factors
#' 
#' ![](images/model_page14_NEW4.png)
#' 
#' ## 
#' 
#' ::: {}
#' 
#' ![](images/intro_5.png)
#' :::
#' 
#' 
#' ## An Additive Model of Gene Expression Including All Factors
#' 
#' ::::: {}
#' 
#' :::: {}
#' 
#' We can use a model ["~ age + sex + tissue + genotype"]{style="background-color: lightyellow"} to describe the full dataset. At least in theory [the model should partition the total variance]{style="background-color: lightyellow"} in the data into the contributions of the experimental variables, such that we can [study the genotype effect unconfounded]{style="background-color: lightyellow"}. 
#' 
## ----reload_dds_additive------------------------------------------------------
#| echo: false
#| eval: true

library(DESeq2)
dds_additive <- readRDS("dds_additive.rds")

#'          
#' ::::
#' 
#' :::: {.fragment}
#' 
## ----dds_additive-------------------------------------------------------------
#| echo: true
#| eval: false

# library(DESeq2)
# 
# dds_additive <-
#   # High level step 1: set up the object structure
#     DESeqDataSetFromMatrix(countData = M,
#                            colData = colData,
#                            design= ~ age + sex + tissue + genotype
#                            ) |>
#     DESeq()
# 

#' ::::
#' 
#' :::: {.fragment}
#' These are the available test results:
#' ::::
#' 
#' :::: {.fragment}
## ----dds_additive_2-----------------------------------------------------------
#| echo: true
#| eval: true
resultsNames(dds_additive)

#' ::::
#' 
## ----save_dds_additive--------------------------------------------------------
#| echo: false
#| eval: false

# library(DESeq2)
# saveRDS(dds_additive, file="dds_additive.rds")

#' 
#' :::::
#' 
#' ## Let's Have a Look at the Genotype Effect in this Model
#' :::: {}
## ----res_additive-------------------------------------------------------------
#| echo: true
#| eval: true

res_additive <- 
    results(dds_additive, name="genotype_nmrhas2_vs_creer")

#' ::::
#' 
#' :::: {.fragment}
## ----res_additive_significant_dummy-------------------------------------------
#| eval: true
#| echo: true

dt <- DT::datatable(
           res_additive |> as.data.frame() |>
           dplyr::filter(!is.na(padj), padj <= 0.05) |>
           dplyr::arrange(padj),

           options=list(scrollY="50000px",
                        lengthMenu = list(c(10, -1), 
                                          c('10','All'))
                        )
     ) 

# convert numeric columns to scientific notation
dt <- dt |>
      DT::formatSignif(sapply(dt$x$data, is.numeric),
                       digits=2
      )

#' 
#' ::::
#' 
#' ## Let's Have a Look at the Genotype Effect in this Model
#' 
## ----res_additive_significant-------------------------------------------------
#| eval: true
#| echo: false

dt

#' 
#' We do get quite some more significant genes than with a genotype-only model. However the majority of them has adjusted p-values very close to the 0.05 threshold. 
#' 
#' ## Let's Have a Look at the Genotype Effect in this Model -- Your Turn!
#' 
#' \         
#' 
#' ::: {style="text-align: center;"}
#' 
#' __`Exercise 5:`__ 
#' 
#' - Use function [quick_GOcheck()]{style="background-color: lightyellow"} [see Exercises] to 
#'   - extract a list of significant genes from res_additive
#'   - feed the list to [clusterProfiler::enrichGO()]{style="background-color: lightyellow"}
#' 
#'   - One [biological topic]{style="background-color: lightyellow"} is highly significant, although only at relatively small absolute log2FoldChange values. What is it?
#' 
#' 
#' :::
#' ##
#' 
#' ::: {}
#' 
#' ![](images/intro_6.png)
#' :::
#' 
#' ## Yet Another Alternative: Group Contrasts
#' 
#' ::: {.incremental}
#' There is another way to include all variables when testing for a genotype effect:
#' 
#' - use the ["group" factor]{style="background-color: lightyellow"} (= the levels of all variables concatenated)
#' - __`compare two level of "group" at a time,`__ \          
#' __`differing only in genotype`__
#' - [this  captures the effect of genotype in the context of the shared levels]{style="background-color: lightyellow"}
#' - notably __`all these contrasts are computed in a single model`__, \       
#' using information from the full dataset
#' 
#' :::
#' 
#' ## Yet Another Alternative: Group Contrasts
#' 
#' ::: {}
#' The idea is embedded in a [framework similar to Anova. Like Anova, it proceeds in two steps:]{style="background-color: lightyellow"}
#' 
#' - First we run __`DESeq() with a 0 + group`__ model. This model tests whether the mean read count at a given "group" level is significantly different from zero. It is only meaningful as raw material for the later tests. 
#' 
#' :::
#' 
#' ::: {.fragment}
## ----load_dds_group-----------------------------------------------------------
#| echo: false
#| eval: true
dds_group <- readRDS(file="dds_group.rds")


#' 
## ----dds_group----------------------------------------------------------------
#| echo: true
#| eval: false
# ## See https://notchbio.app/blog/deseq2-multiple-conditions-contrasts-tutorial/
# 
# dds_group <-
#   DESeqDataSetFromMatrix(countData = M,
#                          colData = colData,
#                          design= ~ 0 + group) |>
#   DESeq()
# 
# 

#' 
## ----save_dds_group-----------------------------------------------------------
#| echo: false
#| eval: false
# saveRDS(dds_group, file="dds_group.rds")
# 

#' 
#' :::
#' 
#' ## Yet Another Alternative: Group Contrasts
#' 
#' ::: {}
#' For the second step we will need a table of the contrasts we want to test:
#' 
## ----load_contrasts-----------------------------------------------------------
#| echo: false
#| eval: false

# contrasts <- readRDS("contrasts.rds")

#' 
## ----contrast_table-----------------------------------------------------------
#| echo: true
#| eval: true

contrasts <-
  data.frame(group= "group",
             A  = colData |> 
                    dplyr::filter(genotype == "nmrhas2") |>
                    dplyr::pull(group),
             B=   colData |> 
                    dplyr::filter(genotype == "creer") |>
                    dplyr::pull(group)
  ) |> dplyr::distinct(.keep_all = TRUE) # there is initially one copy per replicate!

head(contrasts,n=6)
nrow(contrasts)

#' 
## -----------------------------------------------------------------------------
#| echo: false
#| eval: false
# saveRDS(contrasts,file="contrasts.rds")

#' :::
#' 
#' ## Yet Another Alternative: Group Contrasts
#' 
#' ::: {}
#' 
#' - Then in the __`second "Anova" step`__, we call  __`results() on the group means model with a contrast argument`__. This evaluates for every pair of "group" levels in the contrast table, [whether their mean counts are significantly different]{style="background-color: lightyellow"}. 
#' 
#' :::
#' 
#' ::: {.fragment}
#' 
## ----load_res_contrasts-------------------------------------------------------
#| echo: false
#| eval: true

res_contrasts <- readRDS("res_contrasts.rds") 


#' 
## ----res_contrasts------------------------------------------------------------
#| echo: true
#| eval: false

# # a matrix allows to directly extract rows as vectors
# # (unlike a data.frame):
# m <- as.matrix(contrasts)
# 
# res_contrasts <- list()
# for(i in 1:nrow(m)) {
# 
#     res_contrasts[[i]] <- results(dds_group,
#                                   contrast=m[i,]
#                           )
# }
# # the second name in each contrast pair without the trailing ".creer"
# # is the background -- use as name:
# names(res_contrasts)  <- sub("\\.creer$","",contrasts[,"B"])
# 

#' 
#' :::
#' 
#' 
#' ::: {.fragment}
## ----save_res_contrasts-------------------------------------------------------
#| echo: true
#| eval: false

# saveRDS(res_contrasts, file = "res_contrasts.rds")

#' :::
#' 
#' ## Yet Another Alternative: Group Contrasts -- Your Turn!
#' \     
#' 
#' ::: {style="text-align: center;"}
#' 
#' __`Exercise 6:`__ 
#' 
#' - Use function quick_GOcheck() [see Exercises] to 
#'   - extract a list of significant genes from liver-specific res_contrasts entries 
#'   - feed each list to clusterProfiler::enrichGO() and describe your observations
#' 
#' 
#' :::
#' ##
#' 
#' ::: {}
#' ![](images/intro_7.png)
#' :::
#' 
