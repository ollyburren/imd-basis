##' MAKE FIGURE 2
##' clusters of significant projected traits
library(data.table)
library(magrittr)
library(cupcake)

source("R/cw-files.R",echo=TRUE)
source("R/cw-reader.R")
source("R/cw-utils.R")
source("R/cw-dendplot.R")
proj <- reader()


W <- readRDS(BASIS_FILE)
NW <- readRDS(NOWEIGHT_BASIS_FILE)
(load("~/share/as_basis/basis-creation/support/13_trait_basis-degas-comparison.RData"))

BW <- degas.beta
ZW <- degas.z
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
data <- read_raw(traits,rownames(W$rotation))
head(data)
data[,trait:=sub("UKBB_NEALE:SRD:","",trait)]
data <- data[!is.na(trait)]
## data[,center:=W$center[pid]][,sbeta:=shrinkage * beta - center] # center
data[,sbeta:=shrinkage * beta][,z:=beta/seb]
data[,b0:=ifelse(p.value<0.001,beta,0)] #ifelse(p.value>0.001,0,beta)]
data[,z0:=ifelse(p.value<0.001,z,0)] #ifelse(p.value>0.001,0,beta)]
## data[p.value>0.001 | se>0.2, c("lor","z"):=list(0,0)]
## data[indep==TRUE,c("mn","sd"):=list(mean(b0),sd(b0)),by="trait"]
## data[indep==TRUE,zodd:=(beta - mn)/sd,by="trait"]

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

## project
f.project <- function(pca.object, vv) {
    pids <- intersect(names(pca.object$center),data$pid)
    D <- dt2mat(data[pid %in% pids], trait ~ pid, value.var=vv)[,pids]
    D <- rbind(D, ctl=rep(0,ncol(D))) # add control
    D <- D - matrix(pca.object$center[pids],nrow=nrow(D),ncol=ncol(D),byrow=TRUE) # center
    P <- D %*% pca.object$rotation[pids,1:13] # project
    nc <- nrow(P)
    R <- rbind(P[-nc,1:13], pca.object$x[1:13,1:13]) # bind without control
    R - matrix(P[nc,],nrow=nrow(R),ncol=ncol(R),byrow=TRUE) # remove control
}
P <- f.project(W, "sbeta")
NP <- f.project(NW, "beta")
ZP <- f.project(ZW,"z")
BP <- f.project(BW,"beta")

L <- list(noweight=NP,rivas.z=ZP,rivas.beta=BP,weight=P)

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

FPATH="figures" # path to figures dir

## check labels - have we UKBB and not the right way round?

for(nm in names(L)) {
    png(paste0(nm,".png"),height=3.5,width=5,units="in",res=300,pointsize=10)
    ## if(horizontal) {
    par(mar=c(2,0,0,18))
    ## } else {
    ##     par(mar = c(15,2,1,1))
    ## }
    b <- L[[nm]]
    dplotter(b,NULL,NULL,what="b",hclust.method="ward.D2",
             pal=grnvi,
             k=1,show.legend=FALSE)#nm=="weight")
    lab <- letters[which(names(L)==nm)]
    legend("topleft",legend=lab,bty="n")
    dev.off()
}

system(paste("montage -mode concatenate noweight.png rivas.z.png rivas.beta.png weight.png -tile 2x2",
             file.path(FPATH,"figure2-hclust-rivas.png")))

if(interactive())
    system(paste("display", file.path(FPATH,"figure2-hclust-rivas.png")))


q("no")


## sample numbers
library(gsheet)
basis <- gsheet2tbl("https://docs.google.com/spreadsheets/d/1bej1XFq6WwvTlURxhRFmfad2X0BGWZXf7D6IlK6WN0Q/edit#gid=849585942")  %>% as.data.frame()  %>% as.data.table()
neale <- gsheet2tbl("https://docs.google.com/spreadsheets/d/1bej1XFq6WwvTlURxhRFmfad2X0BGWZXf7D6IlK6WN0Q/edit#gid=1742820075")  %>%  as.data.frame()  %>% as.data.table()
neale <- neale[Trait %in% BB_LU13]
head(basis)
head(neale)
basis[,label:=toupper(Trait)  %>%
         sub("LATENT_AUTOIMMUNE_DIABETES_IN_ADULTS","LADA",.)  %>%
         sub("VITILIGO","VIT",.)  %>%
         sub("ASTHMA","asthma",.) %>% 
         sub("IGA_NEPHROPATHY","IgA_NEPH",.)]
neale[,label:=bb.renames[Trait]]
## N <- melt(N,"label")

## take order from L[["weight"]] where everything is nicely paired
## hc <- dist(L[["weight"]])  %>% hclust(.,method="ward.D2")
## N[,label:=factor(label, levels=hc$labels[ hc$order ])]

## order by GWAS sample size
N <- rbind(neale[,.(label,N0,N1,source="UKBB")],basis[,.(label,N0,N1,source="GWAS")])
l <- basis[order(N1)]$label  %>% c(paste0("UKBB_",.),.)  %>% matrix(.,nrow=2,byrow=TRUE)  %>% as.vector()  %>% intersect(., N$label)
N[,label:=factor(label,levels=l)]

theme_set(theme_cowplot())
p <- ggplot(N, aes(x = label, y = N1, fill = source, label = N1)) +
  geom_bar(stat = "identity",position="dodge",colour="grey") +
  geom_text(size = 3,hjust=-0.1) +
  scale_y_continuous("Number of cases",limits=c(0,45000)) +
  scale_x_discrete("") +
  ## scale_y_log10() +
  scale_fill_brewer("Data source",palette="Pastel2") +
  coord_flip() +
  theme(legend.position=c(0.6,0.1),
        axis.text.y=element_text(colour=pcols[levels(N$label)]))
plot_grid(p,labels = "f")
ggsave("samples.png",height=4,width=5,units="in",dpi=300,pointsize=10)


