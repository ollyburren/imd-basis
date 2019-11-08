
library(ggplot2)
library(parallel)
library(pheatmap)
library(cluster)
library(dendextend)
library(factoextra)
library(cowplot)
library(gridGraphics) # recordPlot


dplotter <- function(b,z,p,fixed=NULL,patt=NULL,what=c("signfdr","b"),
                     dist.method="euclidean",hclust.method="ward.D2",k=4,
                     pal=redblue,show.legend=TRUE,horizontal=TRUE) {
    what <- match.arg(what)
    if(!is.null(patt))
        fixed <- c(fixed,grep(paste(patt,collapse="|"),rownames(b),value=TRUE))
    if(!is.null(fixed)) {
        fixed <- intersect(fixed,rownames(b))
        b <- b[fixed,]
        z <- z[fixed,]
        p <- p[fixed,]
    }
    db <- get_dist(b,method=dist.method)
    d <- hclust(db,method=hclust.method)
    if(what=="signfdr") {
        zp <- cut(sign(z) * (1-p),breaks=c(-2,-0.99,-0.95,-0.9,-0.8,0,0.8,0.9,0.95,0.99,2),include.lowest=TRUE)  %>%
          as.numeric()  
        cols <- matrix(pal(10)[zp], 
                       nrow(z),ncol(z),dimnames=dimnames(z))
    } else {
        tb <- t(b)
        tb <- (tb/apply(abs(tb),1,max))  %>% t()
        zp <- cut(tb ,breaks=seq(-1,1,by=0.1),include.lowest=TRUE)  %>%  t()  %>% as.numeric()  
        cols <- matrix(pal(21)[zp], nrow(b),ncol(b),dimnames=dimnames(b))
    }
    
    par(mfrow=c(1,1))
    
    dd <- as.dendrogram(d)
    cuts <- cutree(dd,k=k)
    leafcols <- tol7qualitative[cuts[labels(dd)]]
    ## labels_colors(dd) <-  trait2col[ labels(dd) ]
    dd %>%
      dendextend::set("leaves_pch", 19)  %>%
      dendextend::set("leaves_col", pcols[labels(dd)])  %>% 
      dendextend::set("labels_col", pcols[labels(dd)])  %>% 
      ## dendextend::set("leaves_col", ifelse(grepl("UKBB",labels(dd)), "grey40", "steelblue"))  %>% 
      ## dendextend::set("labels_col", ifelse(grepl("UKBB",labels(dd)), "grey40", "steelblue"))   %>%  
      ## dendextend::set("leaves_col", leafcols) %>%
      ## dendextend::set("branches_k_color", unique(leafcols), k = k)
    plot(.,axes=FALSE,horiz=horizontal)

    ## add coloured bars
    M <- ncol(b)
    ord <- structure(d$order,names=d$order.lab)
    colored_boxes(colors = cols[labels(dd),M:1],dend=dd,sort_by_labels_order=FALSE,horiz=horizontal)
    
    ## tapply(p[fixed,],zp,summary)
    if(show.legend) {
        if(what=="signfdr") {
        legend("topleft",fill=pal(10)[c(1:4,10:7)],
               ncol=2,bty="n",
               text.width=0,
               legend=c(rep("",4),c("<0.01","<0.05","<0.1","<0.2")),
               title="FDR component association")
        } else if(horizontal==TRUE) {
            par(srt = 0)
            legend("topleft",fill=pal(21)[c(1,11,21)],
               bty="n",
               ## legend=c("1","0.5","0",""),
               legend=c("-","0","+"),
               title="component value")
        }
        else {
        legend("topright",fill=pal(21)[c(1,11,21)],
               bty="n",
               ## legend=c("1","0.5","0",""),
               legend=c("-","0","+"),
               title="component value")
    }
        }
    invisible(d)
} 
## dplotter(b,z,p,fixed=tr)

