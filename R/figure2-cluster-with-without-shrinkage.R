##' MAKE FIGURE 2
##' clusters of significant projected traits

## basis p values
## install_github('ollyburren/cupcake')
library(data.table)
library(magrittr)
library(ggplot2)
library(cupcake)
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

################################################################################

source("R/cw-files.R",echo=TRUE)
source("R/cw-reader.R")
proj <- reader()


NW <- readRDS(NOWEIGHT_BASIS_FILE)
W <- readRDS(BASIS_FILE)

## check pid order
identical(rownames(W$rotation),rownames(NW$rotation))

read.basis(basis.noweight)
NW <- NW[setdiff(rownames(NW),"control"),setdiff(colnames(NW),"PC14")]
rownames(NW)  %<>% make.names()
data <- reader("weight")
W <- dt2mat(proj[trait %in% rownames(NW) & PC %in% colnames(NW)],
            trait ~ PC,value.var="delta")[,colnames(NW)]
Wbasis <- data$basis[,colnames(NW)]
W <- rbind(W,Wbasis)[rownames(NW),]

## BB paired traits
BB_LU13 <- c(
  CD = 'crohns.disease',
  CEL = 'malabsorption.coeliac.disease',
  MS = 'multiple.sclerosis',
  RA = 'rheumatoid.arthritis',
  SLE = 'systemic.lupus.erythematosis.sle',
  T1D = 'type.1.diabetes',
  UC = 'ulcerative.colitis',
VIT = 'vitiligo',
  asthma = 'asthma'
)
traits <- BB_LU13[intersect(names(BB_LU13),rownames(W$x))]  %>% paste0("UKBB_NEALE:SRD:",.)
data <- readraw(traits,rownames(W$rotation))
head(data)
data[,trait:=sub("UKBB_NEALE:SRD:","",trait)]
data <- data[!is.na(trait)]
data[,center:=W$center[pid]][,sbeta:=shrinkage * beta - center]

## name nicely for plot
bb.renames <- c("asthma"="UKBB_asthma",
                "crohns.disease"="UKBB_CD",
                "malabsorption.coeliac.disease"="UKBB_CEL",
                "multiple.sclerosis"="UKBB_MS", 
                "rheumatoid.arthritis"="UKBB_RA",
                "systemic.lupus.erythematosis.sle"="UKBB_SLE", 
                "type.1.diabetes"="UKBB_T1D",
                "ulcerative.colitis"="UKBB_UC",
                "vitiligo"="UKBB_VIT")
data[,trait:=bb.renames[trait]]

## cast
D <- dt2mat(data, trait ~ pid, value.var="sbeta")

## project
P <- rbind(D %*% W$rotation[colnames(D),1:13], W$x[1:13,1:13])
NP <- rbind(D %*% NW$rotation[colnames(D),1:13], W$x[1:13,1:13])

paired <- BB_LU13

BB_LU13 <- c(paired,
             structure(setdiff(rownames(W$x), names(paired)),
                       names=setdiff(rownames(W$x), names(paired))))
BB_LU13 <- BB_LU13[ -which(BB_LU13=="control") ]

source("R/cw-colours.R")
patt=NULL      
pcols  <- structure(c(tol9qualitative, tol9qualitative,
                      rep("black",2*length(BB_LU13) - 2*length(paired))),
                    names=c(names(paired),paste0("UKBB_",names(paired)),
                            setdiff(BB_LU13,paired),
                            setdiff(names(BB_LU13),names(paired))))

## BB vs basis on weight vs noweight
L <- list(weight=P,noweight=NP)

FPATH="figures" # path to figures dir

for(nm in names(L)) {
    b <- L[[nm]]
    png(paste0(nm,".png"),height=6,width=5,units="in",res=300,pointsize=10)
    ## if(horizontal) {
        par(mar=c(2,0,0,18))
    ## } else {
    ##     par(mar = c(15,2,1,1))
    ## }
    dplotter(b,NULL,NULL,what="b",hclust.method="ward.D2",
             pal=grnvi,
             k=1,show.legend=FALSE)#nm=="weight")
    if(nm=="weight")
        legend("topleft",legend="b",bty="n")
    if(nm=="noweight")
        legend("topleft",legend="a",bty="n")
    dev.off()
}

system(paste("montage -mode concatenate noweight.png weight.png -tile 2x",
             file.path(FPATH,"figure2-hclust-shrinkage.png")))
unlink(paste0(names(L),".png"))

system(paste("display",
             file.path(FPATH,"figure2-hclust-shrinkage.png")))
