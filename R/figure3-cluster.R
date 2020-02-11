#!/usr/bin/env Rscript
##' MAKE FIGURE 2
##' clusters of significant projected traits

## basis p values
## install_github('ollyburren/cupcake')
## library(pheatmap)
library(cluster)
library(dendextend)
library(factoextra) # get_dist
library(cowplot)
library(gridGraphics) # recordPlot
## library(randomFunctions)
library(magrittr)
library(data.table)
library(ggplot2)

################################################################################

source("R/cw-reader.R")
source("R/cw-colours2.R")
source("R/cw-palette.R")
source("R/cw-renamer.R")
sparse=reader()

proj <- copy(sparse)

#proj[,FDR:=p.adjust(p.value,method="BH"),by="PC"][,stars:=ifelse(FDR<0.01,"*","")]
#proj[,FDR.overall:=p.adjust(p.overall,method="BH"),by="PC"]
proj[grepl("addison",trait),trait:="addisons disease"]# too many characters
proj[grepl("addison",trait),trait.label:="addisons disease"]# too many characters
proj[,trait.label:=sub("UKBB ","",trait.label)]
proj[,lab:=paste(category.label,trait.label) %>% sub("^ ","",.)]

sigs <- proj[fdrcat=="general"][fdr.overall <= ifelse(category.label %in% c("UKBB","Psychiatric"),0.01,1) |
                                trait %in% c("UKBB_NEALE:SRD:ankylosing.spondylitis")]$lab  %>%
                                        # c(.,c("Vasculitis MPO+  1","Vasculitis PR3+  1"))  %>%
        unique()
proj <- proj[lab %in% sigs]
proj[,stars:=ifelse(newfdr<0.01,"*","")]

## number of traits in fig 3
length(unique(proj$trait))


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
tb <- t(b)
zp <- (tb/apply(abs(tb),1,max))  %>% t()  %>% 
    cut(. ,breaks=seq(-1,1,by=0.1),include.lowest=TRUE)  %>%  t()  %>% as.numeric()  
cols <- matrix(grnvi(21)[zp], nrow(b),ncol(b),dimnames=dimnames(b))

D <-  get_dist(b,method="euclidean")  
cl <- hclust(D,method="ward.D2")
dd <- as.dendrogram(cl)
## this group are largely null - push them towards one end of display rather than in middle (topologically this doesn't change dendrogram)
nulls <-  c("UKBB basal cell carcinoma", "UKBB malignant melanoma", 
"UKBB allergy hypersensitivity anaphylaxis", "UKBB eczema dermatitis", 
"UKBB psoriasis", "methotrexate", "UKBB high cholesterol", "UKBB hypertension", 
"UKBB diabetes", "UKBB emphysema chronic bronchitis")
dd %<>% rotate(., order=c(setdiff(labels(dd),nulls),intersect(labels(dd),nulls)))

k <- 4
col4 <- tol5qualitative[c(1,2,3,5)]
## col4 <- c("#015501","#010155","#015501","#010155")
## col4 <- rep(c("dodgerblue","grey30"),length.out=k)
col4 <- rep(c(mygreen,"grey10"),length.out=k)
pch4 <- rep(c(19,18),length.out=k)
cuts <- cutree(dd,k=k)[labels(dd)]
## colbb <- ifelse(grepl("UKBB",labels(dd)),"grey20", "steelblue") #tol5qualitative[1])
colbb <- ifelse(grepl("UKBB",labels(dd)),"grey10", mygreen) #darken(mygreen,1.2)) #tol5qualitative[1])
pchbb <- ifelse(grepl("UKBB",labels(dd)), 18,19)

pdf("~/share/as_basis/figures/figure3-big-cluster.pdf",height=10,width=8,pointsize=10)
par( mar = c(1,0,0,29))
## par(mar = c(28,2,1,2))
dd  %>% 
    dendextend::set("leaves_pch",pchbb)  %>%
    dendextend::set("labels_cex",1)  %>%
    dendextend::set("leaves_col",colbb)  %>% 
  dendextend::set("labels_col",colbb)  %>%
  plot(.,axes=FALSE,horiz=TRUE)
M <- ncol(b)
    ## colored_bars(colors = cols[ord[fixed],M:1],dend=dd,sort_by_labels_order=FALSE)
    ## colored_bars(colors = cols[d$order,M:1],dend=dd,sort_by_labels_order=FALSE)
    ## colored_boxes(colors = cols[labels(dd),M:1],dend=dd,sort_by_labels_order=FALSE)
colored_labelled_boxes(colors = cols[labels(dd),M:1],
                       labels=stars[labels(dd),M:1],
                       dend=dd,sort_by_labels_order=FALSE,horiz=TRUE)
rectcol <- "grey30" #tol5qualitative[1]
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

## sub plot for slides

pdf("~/share/as_basis/figures/slides-big-cluster.pdf",height=8,width=16,pointsize=12)
par( mar = c(25,1,0,0))
## par(mar = c(28,2,1,2))
dd2 <- dd  %>% prune(., c(" birdshot chorioretinopathy",
                          ##"UKBB SLE", "NMO  IgGPos",
                          "NMO  combined", 
## "UKBB vitiligo", "JIA RF+", "Vasculitis MPO+", "Myasthenia gravis early onset", 
                          ## "Myasthenia gravis late onset",
                          "Myasthenia gravis  combined", 
## "JIA undiff ", "JIA EO", "JIA RF-", "JIA PO",
"JIA  combined", 
## "Myositis PM", "Vasculitis PR3+", "Myositis JDM", "Myositis DM", 
                                        "Myositis  combined",
## "JIA PsA", "UKBB rheumatoid arthritis", 
"PsA UK", "PsA Spanish", "UKBB diabetic eye disease",
## "UKBB hyperthyroidism thyrotoxicosis", 
## "UKBB hypothyroidism myxoedema", "UKBB pernicious anaemia", "UKBB type 1 diabetes", 
## "UKBB addisons disease", "EGPA MPO+ ",
"EGPA  combined",
##  "EGPA ANCA- ", 
## "UKBB crohns disease", "UKBB colitis not Crohns or UC", "UKBB ulcerative colitis", 
## "NMO  IgGNeg", "UKBB multiple sclerosis", "JIA systemic", "JIA ERA", 
"PsA N American",
"UKBB malabsorption coeliac disease", "Ank.Spond International", 
## "UKBB sarcoidosis", "UKBB ankylosing spondylitis", "Ank.Spond Turkish/Iranian", 
## " Uveitis", "UKBB nasal polyps", "UKBB asthma", "UKBB hayfever allergic rhinitis", 
 "UKBB basal cell carcinoma", "UKBB malignant melanoma", "UKBB allergy hypersensitivity anaphylaxis", 
"UKBB eczema dermatitis", "UKBB psoriasis", " methotrexate", 
"UKBB high cholesterol", "UKBB hypertension", "UKBB diabetes", 
"UKBB emphysema chronic bronchitis"))
dd2 %>% 
    dendextend::set("leaves_pch",19)  %>%
    dendextend::set("labels_cex",1)  %>%
    ## dendextend::set("leaves_col",colbb)  %>% 
  ## dendextend::set("labels_col",colbb)  %>%
  plot(.,axes=FALSE,horiz=FALSE)
M <- ncol(b)
    ## colored_bars(colors = cols[ord[fixed],M:1],dend=dd,sort_by_labels_order=FALSE)
    ## colored_bars(colors = cols[d$order,M:1],dend=dd,sort_by_labels_order=FALSE)
    ## colored_boxes(colors = cols[labels(dd),M:1],dend=dd,sort_by_labels_order=FALSE)
colored_labelled_boxes(colors = cols[labels(dd2),M:1],
                       labels=stars[labels(dd2),M:1],
                       y_shift=-0.37,
                       dend=dd2,sort_by_labels_order=FALSE,horiz=FALSE)
rectcol <- "grey30" #tol5qualitative[1]
rectx <- 0.15 #0.45 #0.25
recth <- 1 #0.85
recthdiff=-0.3
r <- rect.dendrogram(dd2, k=4,border=rectcol,horiz=FALSE,prop_k_height=0.8,lty=2,
                     lower_rect=-0.61) # -1.46 for 6in wide, -0.5 for 8in wide
text(length(r[[1]])*recth + recthdiff,rectx,labels="I",col=rectcol)
text(length(r[[1]]) + length(r[[2]])*recth + recthdiff,rectx,labels="II",col=rectcol)
text(length(unlist(r[1:2])) + length(r[[3]])*recth + recthdiff,rectx,labels="III",col=rectcol)
text(length(unlist(r[1:3])) + length(r[[4]])*recth + recthdiff,rectx,labels="IV",col=rectcol)
dev.off()    

if(interactive())
    system("evince ~/share/as_basis/figures/slides-big-cluster.pdf &")


