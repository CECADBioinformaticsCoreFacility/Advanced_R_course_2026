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
    rv <- matrixStats::rowVars(assay(object))
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


    ## Add point labels if desired (NOT with a long list of genes ..)
    if(print_labels) {
        p <- p +  geom_text(aes(label=rownames(d),hjust=hjust, vjust=vjust))
    }
    
    ## Set the ratio of the units on x and y axis to 1:
    if(use_coord_fixed) { 
        p <- p + coord_fixed() 
    }

    p
}
