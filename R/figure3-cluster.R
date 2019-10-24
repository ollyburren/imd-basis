#!/usr/bin/env Rscript
##' MAKE FIGURE 2
##' clusters of significant projected traits

## basis p values
## install_github('ollyburren/cupcake')
library(pheatmap)
library(cluster)
library(dendextend)
library(factoextra)
library(cowplot)
library(gridGraphics) # recordPlot

dt2mat <- function(dt,...) {
    tmp <- dcast(dt,...)
    rn <- tmp[[1]]
    m <- as.matrix(tmp[,-1])
    rownames(m) <- rn
    m
}

library(randomFunctions)
library(magrittr)
library(data.table)
library(ggplot2)

################################################################################

source("R/cw-reader.R")
source("R/cw-colours2.R")
source("R/cw-renamer.R")
sparse=reader()

proj <- copy(sparse)

#proj[,FDR:=p.adjust(p.value,method="BH"),by="PC"][,stars:=ifelse(FDR<0.01,"*","")]
#proj[,FDR.overall:=p.adjust(p.overall,method="BH"),by="PC"]
proj[grepl("addison",trait),trait:="addisons disease"]# too many characters
proj[grepl("addison",trait),trait.label:="addisons disease"]# too many characters
proj[,trait.label:=sub("UKBB ","",trait.label)]
proj[,lab:=paste(category.label,trait.label)]

sigs <- proj[fdrcat=="general"][fdr.overall <= ifelse(category.label %in% c("UKBB","Psychiatric"),0.01,1) |
               trait %in% c("UKBB_NEALE:SRD:ankylosing.spondylitis")]$lab  %>% unique()
proj <- proj[lab %in% sigs]
proj[,stars:=ifelse(newfdr<0.01,"*","")]

stars <- dt2mat(proj, trait ~ PC, value.var="stars")[,paste0("PC",1:13)]
obsData <- dt2mat(proj, trait ~ PC, value.var="delta")[,paste0("PC",1:13)]
obsVars <- dt2mat(proj, trait ~ PC, value.var="var.proj")[,paste0("PC",1:13)]^2
fdr <- dt2mat(proj, trait ~ PC, value.var="newfdr")[,1,drop=FALSE]
b <- obsData

renamer <- unique(proj[,.(category,trait,trait.label,category.label,lab)])
renamer[is.na(lab)| is.na(trait)]


m <- match(rownames(b), renamer$trait)
summary(m)
rownames(stars) <- rownames(b) <- sub("_combined","combined",renamer$lab[m])

## clust plot
    par(mar = c(30,2,1,1),mfrow=c(1,1))
## z <- b/sqrt(obsVars)
## z <- matrix(ifelse(abs(z)>1.96,sign(z),0),nrow(b),ncol(b))
## dimnames(z) <- dimnames(b)
## db <- get_dist(z,method="euclidean")

tb <- t(b)
zp <- (tb/apply(abs(tb),1,max))  %>% t()  %>% 
    cut(. ,breaks=seq(-1,1,by=0.1),include.lowest=TRUE)  %>%  t()  %>% as.numeric()  
cols <- matrix(grnvi(21)[zp], nrow(b),ncol(b),dimnames=dimnames(b))



## source("~/Projects/auto-basis/scripts/files13.R")
## pc.emp <- readRDS(BASIS_FILE)

## norm <- function(x) {
##     ## ((t(x) - colMeans(x))/apply(x,2,sd))  %>% t()
##     ( t(x)/pc.emp$sdev[1:ncol(x)] )  %>% t()
## }

D <-  get_dist(b,method="euclidean")  
cl <- hclust(D,method="ward.D2")
dd <- as.dendrogram(cl)
## sapply(2:10, function(i) dunn(D, cutree(cl,i)))

## plot(dd)

k <- 4
col4 <- tol5qualitative[c(1,2,3,5)]
## col4 <- c("#015501","#010155","#015501","#010155")
col4 <- rep(c("royalblue","grey30"),length.out=k)
pch4 <- rep(c(19,18),length.out=k)
cuts <- cutree(dd,k=k)[labels(dd)]
colbb <- ifelse(grepl("UKBB",labels(dd)), "grey30",tol5qualitative[5])
pchbb <- ifelse(grepl("UKBB",labels(dd)), 18,19)

## if(J==1)
## if(J==2)
## pdf("~/basis-11-big-cluster.pdf",height=6,width=12,pointsize=10)

cuts[!duplicated(cuts)]
## cuts[cuts==2] <- 0
## cuts[cuts==1] <- 2
## cuts[cuts==0] <- 1


#split.screen(c(1,2)

## mat <- matrix(c(1,2), 1,2)
## mat
## layout(mat, c(3,1), c(1))
pdf("~/share/as_basis/figures/figure3-big-cluster.pdf",height=12,width=8,pointsize=10)
par( mar = c(3,2,1,29))
## par(mar = c(28,2,1,2))
dd  %>% 
    ## color_branches(., k=k, col=col4)  %>%  #[c(1,3,5,7)])  %>%
    ## color_labels(., k=k, col=col4)  %>%  #[c(1,3,5,7)])  %>%
    ## dendextend::set("leaves_pch",pch4[1])  %>%
    dendextend::set("leaves_pch",pchbb)  %>%
    dendextend::set("leaves_col",colbb)  %>% 
    dendextend::set("labels_col",colbb)  %>% 
    ## dendextend::set("leaves_col",col4[cuts])  %>% 
    ## dendextend::set("labels_col",col4[cuts])  %>% 
  plot(.,axes=FALSE,horiz=TRUE)
M <- ncol(b)
    ## colored_bars(colors = cols[ord[fixed],M:1],dend=dd,sort_by_labels_order=FALSE)
    ## colored_bars(colors = cols[d$order,M:1],dend=dd,sort_by_labels_order=FALSE)
    ## colored_boxes(colors = cols[labels(dd),M:1],dend=dd,sort_by_labels_order=FALSE)
colored_labelled_boxes(colors = cols[labels(dd),M:1],
                       labels=stars[labels(dd),M:1],
                       dend=dd,sort_by_labels_order=FALSE,horiz=TRUE)
rectcol <- tol5qualitative[1]
rectx <- 0.15 #0.45 #0.25
recth <- 1 #0.85
recthdiff=-0.3
r <- rect.dendrogram(dd, k=4,border=rectcol,horiz=TRUE,prop_k_height=0.8,lty=2,
                     lower_rect=-0.61) # -1.46 for 6in wide, -0.5 for 8in wide
text(rectx,length(r[[1]])*recth + recthdiff,labels="I",col=rectcol)
text(rectx,length(r[[1]]) + length(r[[2]])*recth + recthdiff,labels="II",col=rectcol)
text(rectx,length(unlist(r[1:2])) + length(r[[3]])*recth + recthdiff,labels="III",col=rectcol)
text(rectx,length(unlist(r[1:3])) + length(r[[4]])*recth + recthdiff,labels="IV",col=rectcol)
dev.off()    

if(interactive())
    system("evince ~/share/as_basis/figures/figure3-big-cluster.pdf &")
## ## fdr0
## par(mar=c(3,25,1,2),mgp=c(2,1,0),mfg=c(1,2))
## fdr[fdr<1e-4] <- 1e-4
## rownames(fdr) <- NULL
## barplot(-t(log10(fdr[order.dendrogram(dd),,drop=FALSE])),horiz=TRUE,xlab="-log10 FDR",col="grey70",axes=FALSE)
## axis(1,at=c(0,2,4),labels=c("1","0.01","<0.0001"))



