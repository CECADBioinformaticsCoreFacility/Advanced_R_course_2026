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
#' **Session 2:: Differential Gene Expression:: Basics **
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


#' 
#' ## 
#' 
#' ::: {}
#' ![](images/intro_1.png)
#' :::
#' 
#' ## The Practice Dataset
#' 
#' ::::: {.columns}
#' 
#' :::: {.column width="50%"}
#' 
#' ![](images/gorbunova.png)
#' ::::
#' :::: {.column width="50%"}
#' :::{.fragment}
#' ![](images/zhang2023_abstract.png)
#' :::
#' ::::
#' :::::
#' 
#' ## The Practice Dataset
#' ::::: {.columns}
#' :::: {.column width="50%"}
#' ![](images/zhang2023_GEO.png)
#' ::::
#' 
#' :::: {.column width="50%"}
#' ::: {.fragment}
#' \ 
## ----download_matrix----------------------------------------------------------
#| echo: TRUE
data_dir <- "." ## the current folder
matrix_file <- "GSE234563_Matrix.txt.gz"

# download it ..
#download.file(
#  url="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE234563&format=file&file=GSE234563%5FMatrix%2Etxt%2Egz",
#  destfile=matrix_file
#)

#' :::
#' ::: {.fragment}
## ----read_matrix--------------------------------------------------------------
#| echo: TRUE
# .. read it:
M <- read.table(fs::path(data_dir,
                         matrix_file),
                header = TRUE,
                sep = "\t",
                quote = ""
                )
# .. this is a huge matrix!
dim(M)

#' :::
#' 
#' ::: {.fragment}
## ----inspect_matrix_dummy-----------------------------------------------------
#| echo: true
#| eval: false
#| 
# # .. display the first 6 rows (output next page):
# dt <- DT::datatable(M |> head())
# dt
# 

#' :::
#' ::::
#' :::::
#' 
#' ## The Practice Dataset
#' ::::: {}
## ----inspect_matrix-----------------------------------------------------------
#| echo: false
#| eval: true

dt <- DT::datatable(M |> head())
dt

#' :::::
#' 
#' ## Count Matrix Preprocessing
#' ::::: {}
#' \         
#' 
#' :::: {} 
#' First of all, we **convert the column names to all-lowercase**, to avoid case issues:
## ----tolower------------------------------------------------------------------
#| echo: true

colnames(M) <- tolower(colnames(M))

#' ::::
#' \         
#' 
#' :::: {.fragment}
#' Then we **discard genes with a total count < 10**, as was done in the paper: 
#' 
## ----count_filter1------------------------------------------------------------
#| echo: true
#| eval: true

min_count <- 10

M <- M[rowSums(M) >= min_count,] 
#M <- M[(rowSums(M) >= min_count) &  ## for testing
#      apply(M,1,function(x) !any(x >= 10^3)) ,]

nrow(M)



#' ::::
#' \         
#' 
#' :::: {.fragment}
#' Finally we **map matrix entries with fractional digits to the next higher integer value**\        
#' (the differential expression analysis expects integer counts):
#' 
## ----count_filter2------------------------------------------------------------
#| echo: true
#| eval: true
M <- ceiling(M)

#' ::::
#' :::::
#' 
#' ## Introducing the DESeq2 Package
#' ::::: {}
#' :::: {.incremental}
#' - **A Bioconductor package for Differential Gene Expression (DGE) analysis**         
#' - Based on sophisticated mathematical theory __`(more on that later by Ali)`__      
#' - **But its high-level functions allow to run an entire DGE analysis in three steps**!             
#' - See `browseVignettes("DESeq2")`           
#' - Read [the original paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)
#' 
#' ::::
#' 
#' :::: {.fragment}
## ----DESEq2_1-----------------------------------------------------------------
#| echo: true
#| eval: false

# library(DESeq2)
# 
# browseVignettes("DESeq2")

#' ::::
#' :::::
#' 
#' ## Setting Up the Input for DESeq2
#' 
#' ::::: {}
#' ![](images/DESeq2_input1.png)
#' :::::
#' 
#' ## Setting Up the Input for DESeq2
#' ::::: {}
#' ![](images/DESeq2_input2.png)
#' :::::
#' 
#' ## Setting Up the Input for DESeq2
#' ::::: {}
#' 
#' We do not have the colData -- but it is actually easy to create this matrix starting from the count matrix column names:
#' 
#' :::: {.fragment}
## ----build_colData1-----------------------------------------------------------
#| echo: true
#| eval: true

# Step 1: split the count matrix column names at "_"
#         (result: a list of vectors, where each vector holds the fragments extracted from one name)
colData <- strsplit(colnames(M),"_") 

head(colData, n=2) 

#' ::::
#' \            
#' 
#' :::: {.fragment}
## ----build_colData2-----------------------------------------------------------
#| echo: true
#| eval: true
#| 
# Step 2: make a matrix, in which each vector from Step 1 is a row
colData <- Reduce(rbind,colData)

head(colData, n=2) 

#' ::::
#' 
#' :::::
#' 
#' ## Setting Up the Input for DESeq2
#' 
## ----build_colData3-----------------------------------------------------------
#| echo: true
#| eval: true
#| 
# Step 3: We know what the row and column names should be, so let's set them:
dimnames(colData) <- list(colnames(M),
c("age","sex","genotype","tissue","replicate")
                          )
head(colData, n=2) 

#' 
#' ## Setting Up the Input for DESeq2
#'             
## ----build_colData4-----------------------------------------------------------
#| echo: true
#| eval: true
#| 

# Step 4: Add two grouping factors, which combine the levels of > 1 variable:  

colData <- 
  colData |> as.data.frame() |> 
  tidyr::unite("background",  c(age,sex,tissue),
               sep=".",remove=FALSE) |> 
  tidyr::unite("group",c(age,sex,tissue,genotype), 
               sep=".",remove=FALSE)


# Have look at the result:
head(colData,n=2) 

#' 
## ----save_for_Session4--------------------------------------------------------
#| echo: false
#| eval: true

saveRDS(colData, file="colData.rds")


#' 
#' ## This is a HUGE Matrix !
#' 
#' ![](images/dataset_structure1.png)
#' 
#' ## This is a HUGE Matrix !
#' 
#' The simplest way to deal with it: __`subset!`__
#' 
#' ![](images/dataset_structure2.png)
#' 
#' ## This is a HUGE Matrix !
#' 
#' The second-simplest way: __`group!`__
#' 
#' ![](images/dataset_structure3.png)
#' 
#' 
#' 
#' ## A Bold Try: Model the Full Data With Genotype Only
#' 
#' :::: {}
#' The most simple question we can ask to our data matrix is: \           
#' **Is the nmrHas transgene effect strong enough to manifest irrespective of age, sex and tissue?**  
#' Starting from this question, we will dig deeper both into the dataset and into modeling.
#' \                    
#' 
## ----reload_dds_genotype------------------------------------------------------
#| echo: false
#| eval: true

library(DESeq2)
dds_genotype <- readRDS("dds_genotype.rds")

#' \         
#' 
#' ::: {.fragment}
#' 
## ----dds_genotype_1-----------------------------------------------------------
#| echo: true
#| eval: false

# library(DESeq2)
# 
# dds_genotype <-
#   # High level step 1: set up the object structure
#     DESeqDataSetFromMatrix(countData = M,
#                            colData = colData,
#                            design= ~ 1 + genotype
#                            ) |>
#   # High level step 2: "do the math" (--> Ali!)
#   #                    size factors factors,
#   #                    dispersion estimates,
#   #                    statistical tests ...
#     DESeq()
# 
# 

#' :::
#' 
#' ::: {}
## -----------------------------------------------------------------------------
#| eval: false
#| echo: false

# saveRDS(dds_genotype, file = "dds_genotype.rds")

#' :::
#' 
#' ::: {.fragment}
#' At this point we can have the available test results listed:
#' :::
#' 
#' ::: {.fragment}
## ----dds_genotype_2-----------------------------------------------------------
#| echo: true
#| eval: true
resultsNames(dds_genotype)

#' :::
#' 
#' ::::
#' 
#' ## A Bold Try: Model the Full Data With Genotype Only
#' 
#' :::: {.fragment}
#' We want to see the genotype effect ...
## ----dds_genotype_4-----------------------------------------------------------
#| echo: true
#| eval: true

res_genotype <- 
    results(dds_genotype, 
            name="genotype_nmrhas2_vs_creer"
            )
saveRDS(res_genotype, file="res_genotype.rds")

#' 
#' 
#' ::::
#' 
#' :::: {.fragment}
#' ... in a nice table (see next slide):
## ----res_genotype_significant_dummy-------------------------------------------
#| eval: true
#| echo: true

dt <- DT::datatable(
           res_genotype |> as.data.frame() |>
           dplyr::filter(!is.na(padj), padj <= 0.05) |>
           dplyr::arrange(padj),

           options=list(scrollY="500px")
     ) 

# convert numeric columns to scientific notation
dt <- dt |>
      DT::formatSignif(
           sapply(dt$x$data, is.numeric),
           digits=2
      )

#' ::::
#' 
#' 
#' ## A Bold Try: Model the Full Data With Genotype Only
#' 
## ----res_genotype_significant-------------------------------------------------
#| eval: true
#| echo: false

dt

#' 
#' ## MA Plot and Shrinkage
#' 
#' The __`plotMA()`__ function of DESeq2 draws a [scatterplot of the mean normalized counts against the log2FoldChange values]{style="background-color: lightyellow"}. It helps to see the significant genes in the context of the entire dataset, and notably also to pinpoint unusual features of the entire count distribution.
#' 
## ----res_genotype_MAplot------------------------------------------------------
#| echo: true
#| output-location: column-fragment

plotMA(res_genotype) 


#' 
#' ::: {.fragment .incremental}
#' Observe:
#' 
#' - the number of significant genes is REALLY small compared to the entire distribution  
#' - there is a slight bimodality in the point cloud -- this is unusual
#' 
#' :::
#' 
#' ## MA Plot and Shrinkage
#' 
#' The __`lfcShrink()`__ function of DESeq2 
#' 
#' - is an alternative to results()
#' - [shifts noisy, unreliable log2FoldChange values towards zero]{style="background-color: lightyellow"}
#' - see __`?lfcShrink`__ and the DESeq2 vignette for available shrinkage methods
#' 
## ----res_genotype_shrunken----------------------------------------------------
#| echo: true
#| output-location: column-fragment

res_genotype_shrunken <-
  lfcShrink(dds_genotype, 
            ## MIND: "coef", 
            ##       not "name" as in results() !
            coef="genotype_nmrhas2_vs_creer", 
            type="ashr")

plotMA(res_genotype_shrunken)


#' 
#' 
#' ::: {.fragment .incremental}
#' The bimodality is much more clear after shrinkage -- this is really unusual
#' :::
#' 
#' 
#' 
#' ## 
#' ::: {}
#' 
#' ![](images/intro_2.png)
#' :::
#' 
#' ## Why So Few Significant Genes?
#' 
#' 
#' ## Why So Few Significant Genes?
#' 
#' ::: {}
#' I had on purpose omitted a step which should really be done in the beginning: __`Get to know the sources of variability in your count matrix!`__ 
#' \          
#' 
#' All known such sources should be either eliminated or stated in the design given to DESeq2.  Otherwise the **neglected factors will be treated as random variance** by the model and **this may mask a real effect**.  
#' :::
#' 
#' ## Why So Few Significant Genes?
#' 
#' ::: {}
#' __`Principal Component Analysis (PCA)`__ projects the transcriptome of each sample onto a point on a two-dimensional canvas, such that the among-sample variance is maximal. 
#' 
#' It helps to answer questions like:
#' 
#' :::
#' 
#' ::: {.incremental}
#' 
#' - Do the samples cluster as your hypotheses suggest? 
#' - If not, does the cluster structure suggest __`hidden factors`__ which can be identified and included in the model? 
#' - Or does __`poor clustering of expected groups`__ suggest experimental issues?
#' :::
#' 
#' ## Why So Few Significant Genes?
#' 
#' ::: {}
#' __`Heatmaps`__ display the transcriptomes directly as colored matrix columns. Clustering the columns can identify groups of samples with similar global gene expression, like PCA.
#' :::
#' 
#' ::: {.fragment}
#' __`Boxplots`__ of the matrix columns visualize the distributional properties of the per-sample counts.
#' :::
#' 
#' ::: {.fragment}
#' ... and of course the __`diagnostic tools of the DESeq2 package`__:
#' 
#' - __`MA plot`__: graphs the mean normalized count agains the log2FoldChange           
#' - __`meanSdPlot`__: [tries to transform the count distribution such that the standard deviation is ~identical for ~all of the genes]{style="background-color:lightyellow"}     
#'   - DESeq2's differential gene expression analysis assumes that this is possible
#'   - [if it is not possible, then the DGE results will likely be misleading!]{style="background-color:lightyellow"}
#' :::
#' 
#' 
#' ## Preparation for PCA
#' 
#' DESeq2's  __`vst()`__ function  implements this __`variance stabilizing transformation`__.\   
#' [It is meant to eliminate the increase in variance with count which is natural to count data.]{style="background-color:lightyellow"} 
#' 
#' We will use vst() to transform the input to PCA and heatmap functions.
#' 
#' 
## ----vsd----------------------------------------------------------------------
#| echo: true
#| eval: true

# The "blind" parameter tells the function NOT to use the input object's
# model formula for identifying "valid" (expected) variation. 
# We do this here because we have not yet found an optimal formula.
# NOTE that if we have a trusted formula and want to use vst() to 
# prepare an object for downstream analysis, then blind=FALSE should 
# be used.
vsd_blind <- vst(dds_genotype, 
                 blind=TRUE
                )
saveRDS(vsd_blind,file="vsd_blind.rds")

#' 
#' ## A Barebones PCA Plotting Function
#' 
#' :::{}
#' The __`DESeq2 package`__ also provides a function __`plotPCA()`__, which executes a Principal Components Analysis and plots the result. It is an implementation of a __`generic function`__ `BiocGenerics::plotPCA`. To see its code, type
#' 
## -----------------------------------------------------------------------------
#| echo: true
#| eval: true
#| output-location: fragment

# The triple colon operator shows non-exported code of a package
DESeq2:::plotPCA.DESeqTransform


#' 
#' :::
#' 
#' ## A Barebones PCA Plotting Function
#' 
#' :::{}
#' I provide a simplified version of it, which removes much of the pre-processing code.\         
#' 
#' Instead, the function relies on the __`dot-dot-dot`__ argument (__`...`__),\       
#' which is a list of __`parameter=argument pairs`__. \     
#' 
#' It is kind of __`a self-service store for functions called in the`__\     
#' __`body of a function`__. \           
#' 
#' ::: 
#' \                
#' 
#' ::: {.fragment}
#' 
## -----------------------------------------------------------------------------
#| output-location: column-fragment
#| error: true
#| echo: true

try({
f2 <- function(a) cat(a,"\n")
f1 <- function(i,j) f2() 

## Simply passing a value for f2's "a" 
## through f1's interface fails
f1(i=1,j=2,a="hello world!")

})

#' :::
#' 
#' ::: {.fragment}
## -----------------------------------------------------------------------------
#| output-location: column-fragment
#| echo: true

## but if f1() has dot-dot-dot in its interface 
## and passes it on to f2() ---
f1 <- function(i,j, ...) f2(...)

#' :::
#' 
#' ::: {.fragment}
## -----------------------------------------------------------------------------
#| output-location: column-fragment
#| echo: true

## --- it succeeds! 
f1(i=1,j=2,a="hello world!")

#' :::
#' 
#' 
#' 
#' ## A Barebones PCA Plotting Function
## ----plotpca------------------------------------------------------------------
#| echo: true

plotPCA2 <-
    
function (object, ntop = 500, returnData = FALSE,
          pcsToUse=c(x=1,y=2),
          pointSize=3,
          hjust=-0.75,vjust=0.5,
          just_print_select=FALSE,
          print_labels=FALSE, 
          use_coord_fixed=FALSE,
          ...
          )  {

    ## Get the ntop most highly variable genes:  
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                       length(rv)))]
    if(just_print_select) {
        return(select)
    }
    
    ## Compute the PCA:
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)

    d <- data.frame(pca$x[, pcsToUse[1]], pca$x[, pcsToUse[2]], 
                    as.data.frame(colData(object)),
                    name = colnames(object))
    colnames(d)[1:2] <- paste("PC",pcsToUse,sep="")

    if (returnData) {
        attr(d, "percentVar") <- percentVar[pcsToUse]
        return(d)
    }
    
    ## Draw the scatterplot:
    p <- 
        ggplot(data = d,
               aes_string(x = paste("PC",pcsToUse["x"],sep=""),
                          y = paste("PC",pcsToUse["y"],sep="")
                          )
        ) + 
        geom_point(mapping = aes(
                      ...
                   ),
                   size = pointSize  ## set last, to override any "size"
                                     ## argument set inside aes()
                   
        ) +
        xlab(paste0("PC",pcsToUse["x"],": ", round(percentVar[pcsToUse["x"]] * 100), "% variance",sep="")) + 
        ylab(paste0("PC",pcsToUse["y"],": ", round(percentVar[pcsToUse["y"]] * 100), "% variance",sep=""))


    ## Add point labels if desired 
    if(print_labels) {
        p <- p +  geom_text(aes(label=rownames(d),hjust=hjust, vjust=vjust))
    }
    
    ## Set the ratio of the units on x and y axis to 1:
    if(use_coord_fixed) { 
        p <- p + coord_fixed() 
    }

    p
}

#' 
#' ## What Drives the Structure?
#' 
## ----PCA1---------------------------------------------------------------------
#| echo: true
#| eval: true

plotPCA2(vsd_blind,
         ntop=500, 
         pcsToUse=c(x=1,y=2))

#' There is very clear structure in the  unlabelled dataset ...
#' 
#' ## What Drives the Structure?
#' 
#' ::: {}
## ----PCA2---------------------------------------------------------------------
#| echo: true
#| eval: true
#| output-location: column

plotPCA2(vsd_blind,
         color=genotype,
         ntop=500, 
         pcsToUse=c(x=1,y=2) 
         )  
  

#' :::
#' 
#' ::: {.fragment}
#' ... but genotype is more or less unrelated to this global structure!
#' :::
#' 
#' ::: {.fragment}
#' __`What about the other variables?`__
#' :::
#' 
#' ## What Drives the Structure?
#' ::: {}
## ----PCA3---------------------------------------------------------------------
#| echo: true
#| eval: true
#| output-location: column-fragment

plotPCA2(vsd_blind,
         ntop=1000, 
         color = sex,
         shape = age, 
         alpha = genotype,
         pcsToUse=c(x=1,y=2) 
         
         ) + 
  # mark empty-vector samples as transparent,
  # NMR transgene samples as opaque:
  scale_alpha_manual(values = c("creer" = 0.3, 
                                "nmrhas2" = 1)) 

#' :::
#' 
#' ## What Drives the Structure? -- Your Turn!
#' 
#' \         
#' 
#' 
#' ::: {style="text-align: center;"}
#' 
#' __`Exercises 1+2:`__ 
#'       
#' - After variance stabilization, the standard deviation of the counts should be roughly independent of their mean. Can this be achieved with the full dataset?
#' - Visualize the distribution of the experimental variables on the PCA plot  
#'   - [Do any of the variables correlate strongly with the PCA structure?]{style="background-color: lightyellow"}
#' 
#' 
#' :::
#' 
#' ## What Drives the Structure?
#' ::: {}
## ----PCA4---------------------------------------------------------------------
#| echo: true
#| eval: true
#| output-location: column-fragment

plotPCA2(vsd_blind,
         ntop=500, 
         color = sex,
         alpha = genotype,
         pcsToUse=c(x=1,y=2) 
         
         ) + 
  # mark empty-vector samples as transparent,
  # NMR transgene samples as opaque:
  scale_alpha_manual(values = c("creer" = 0.3, 
                                "nmrhas2" = 1))

#' :::
#' 
#' ## What Drives the Structure?
#' ::: {}
## ----PCA5---------------------------------------------------------------------
#| echo: true
#| eval: true
#| output-location: column-fragment

plotPCA2(vsd_blind,
         ntop=500, 
         color = paste(sex,age,sep="."),
         alpha = genotype,
         pcsToUse=c(x=1,y=2) 
         
         ) + 
  # mark empty-vector samples as transparent,
  # NMR transgene samples as opaque:
  scale_alpha_manual(values = c("creer" = 0.3, 
                                "nmrhas2" = 1))

#' :::
#' 
#' 
#' ## What Drives the Structure?
#' 
#' ::: {}
#' 
## ----PCA6---------------------------------------------------------------------
#| echo: true
#| eval: true
#| output-location: column-fragment

plotPCA2(vsd_blind,
         ntop=500, 
         color = tissue,
         shape = paste(sex, age, sep="."),
         alpha = genotype,
         pcsToUse=c(x=1,y=2) 
         
         ) + 
  # mark empty-vector samples as transparent,
  # NMR transgene samples as opaque:
  scale_alpha_manual(values = c("creer" = 0.3, 
                                "nmrhas2" = 1))

#' 
#' :::
#' 
#' ## What Drives the Structure?
#' 
#' ::: {}
#' 
## ----PCA6_ntop1000------------------------------------------------------------
#| echo: true
#| eval: true
#| output-location: column-fragment

plotPCA2(vsd_blind,
         ntop=1000, 
         color = tissue,
         shape = paste(sex, age, sep="."),
         alpha = genotype,
         pcsToUse=c(x=1,y=2) 
         
         ) + 
  # mark empty-vector samples as transparent,
  # NMR transgene samples as opaque:
  scale_alpha_manual(values = c("creer" = 0.3, 
                                "nmrhas2" = 1))

#' 
#' :::
#' 
#' ## Heatmap of the 500 Most Highly Expressed Genes
#' \               
#' 
#' ::: { .fragment}
#' 
## ----heatmap_vsd500-----------------------------------------------------------
#| echo: true
#| eval: true
#| output-location: column

selected <- 
  order(rowMeans(counts(dds_genotype,
                        normalized=TRUE)
                ),
        decreasing=TRUE)[1:500]

df <- (colData(dds_genotype) |>
       as.data.frame())[,c("age",
                           "sex",
                           "genotype",
                           "tissue")]

pheatmap(assay(vsd_blind)[selected,], 
         cluster_rows=FALSE, 
         show_rownames=FALSE,
         show_colnames=FALSE,
         cluster_cols=TRUE, 
         annotation_col=df)

#' 
#' :::
#' 
#' ::: { .fragment .incremental}
#'           
#' - With the exception of spleen and WAT, all tissues have strongly distinct expression patterns
#' - Both levels of the genotype factor are compatible with any of these global patterns
#' :::
#' 
#' ## Boxplots of the Per-Sample Count Distributions 
#' 
#' We want to `plot the raw and normalized count distributions of each sample side by side, using ggplot()`. For this we have to "melt" the log10 transformed matrices:
#' 
#' ::: {.fragment}
#' 
## ----count_boxplots1----------------------------------------------------------
#| echo: true
#| eval: true

df <-
    rbind(reshape2::melt(log10(counts(dds_genotype,normalized=FALSE)+1)) |>
          mutate(mode="raw"),

          reshape2::melt(log10(counts(dds_genotype,normalized=TRUE)+1)) |>
          mutate(mode="normalized")
          ) |>
  # rename to make variable names more meaningful:
    dplyr::rename(gene=Var1, sample=Var2, count=value) |>
  
  # remove the trailing sample ID to create a common ID for all replicates
  # of a given factor level combination:
    mutate(group=sub("_\\d+$","",sample)) |>
  
  # ggplot needs to know the order of sub-plots ("raw" should come first)
    mutate(mode = factor(mode,
                         levels=c("raw","normalized")
                         )
         )
           


#' :::
#' 
#' ## Boxplots of the Per-Sample Count Distributions 
#' 
#' Then plot it:
#' 
## ----count_boxplots2_dummy----------------------------------------------------
#| echo: true
#| eval: false

# ggplot(df,
#        # aes(y=count,x=sample,fill=group) ## too dense
#        aes(x=count,y=sample,fill=group)
#        ) +
#     geom_boxplot(show.legend=FALSE) +
#     facet_wrap(vars(mode),ncol=2) +
#     theme_classic() +
#     theme(axis.text=element_text(size=2)) #,
#           #axis.text.x = element_text(angle = 90,
#           #                           vjust = 0.5,
#           #                           hjust=1))

#' 
#' ## Boxplots of the Per-Sample Count Distributions 
#' 
#' ::: {}
#' ![](images/count_boxplots.png)
#' :::
#' 
#' ## Boxplots of the Per-Sample Count Distributions 
#' 
#' ::: {.incremental}
#' - While raw counts may differ in shape and location between samples, these **differences should go away after a successful [DESeq2 "median of ratios" normalization](https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html)**.
#' 
#' -  __`Failure to normalize indicates an issue!`__. 
#' \           
#' 
#' - Here, the reason is very likely the fact that **the samples come from different tissues, with very different gene expression signatures**. DESeq2's approach to use the __`geometric expression mean over all samples of each gene as a pseudo-reference count may not be justified in this case`__. 
#' 
#' :::
#' 
#' ##
#' 
#' ::: {}
#' 
#' ![](images/intro_3.png)
#' :::
#' 
#' ## Subset the Dataset to One Tissue At a Time
#' 
#' ::: {style="font-size: 50%;"}
#' The below function accepts a single metadata variable `use_var` and a vector `use_levels` of its levels. It returns a new dds object, which
#' 
#' - contains only samples from the selected levels of the selected variable
#' - has appropriate sizeFactors
#' - no longer contains genes whose total count is less than min_count after sample restriction     
#' - no longer contains/needs `use_var` -- all samples are now from `use_var` by definition
#' 
#' :::
#' \          
#' 
## ----subset_dds---------------------------------------------------------------
#| echo: true
#| eval: true

subset_dds <- 
  function(dds, use_var, use_levels, 
           design = NULL, min_count=10) {
    
    M <- counts(dds,normalize=FALSE)
    colData <- colData(dds)
    
    if(is.null(design)) { # use the design of the input object
      use_design <- design(dds) 
    } else { # set a custom design
      use_design = design
    }
    
    j <- colData[,use_var] %in% use_levels ## samples to keep
    i <- rowSums(M[,j]) >= min_count       ## rows to keep
    
    DESeqDataSetFromMatrix(countData = M[i,j],
                           colData = colData[j,] |> ## natively a DFrame ! 
                             as.data.frame() |>
                             dplyr::select(-!!use_var),
                           design= as.formula(use_design)
    ) |> DESeq()
  }


#' 
#' 
#' 
#' ## Subset the Dataset to One Tissue At a Time -- Your Turn!
#' 
#' 
#' ::: {style="text-align: center; font-size: 80%"}
#' 
#' __`Exercise 3:`__ 
#'       
#' - For each level of the tissue variable,\          
#'   [construct a dds object with samples exclusively from this level]{style="background-color: lightyellow"}
#' - Collect the per-tissue dds objects in a list, and save it in an RDS file for later use
#' - For each such dds object
#'   - [run vst(object) |> plotPCA2()]{style="background-color: lightyellow"}
#'   - mark sex on the plot by color, age by shape, and genotype by alpha
#' - [Can we now see more fine structure?]{style="background-color: lightyellow"} 
#' 
#' 
#' __`Exercise 4:`__ 
#'       
#' - For each per-tissue dds object, 
#'   - draw an MA plot with and without shrinkage
#'   - tabulate the number of significant genes (padj <= 0.05) with and without shrinkage
#'   - [Is there evidence of tissue-specific gene expression programs?]{style="background-color: lightyellow"}
#' 
#' 
#' :::
#' 
#' 
#' ## What Did We Learn in Part 1?
#' 
#' ::: {.incremental}
#' - The __`different tissues have very different transcriptomes`__ (as can be expected!)
#' - It is therefore __`sensible to analyse them one at a time`__, 
#' - unless (-> Part 2): 
#'   - [we include all experimental variables in the model formula!]{style="background-color: lightyellow"}
#' :::
#' 
#' ## Before We Go On ... 
#' 
#' ## ... Ali On the Mathematics Behind DESeq2!
#' 
#' ##
#' 
