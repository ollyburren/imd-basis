##' look at eqtl/pqtls and pathways
##' 2019-07-23 projections using only FDR<5% variants result in a load of genes at 0. means we can't use empirical significance.
##' 2019-08-30 add Roederer
##' rerun on 13-basis
##' 2019-10-24 centre betas

## basis p values
## install_github('ollyburren/cupcake')
library(data.table)
library(magrittr)
library(ggplot2)
library(cupcake)
library(parallel)
source("R/cw-reader.R")
source("R/cw-files.R")

proj <- reader()
table(proj$category)

shrink.DT<-readRDS(SHRINKAGE_FILE)

(load(file="~/share/as_basis/sparse-basis/basis-sparse-0.999.RData"))

## basis <- x$basis  %>% reshape2::melt()  %>% as.data.table()
## setnames(basis,c("trait","PC","value"))

## functions to calc consistency and plot
library(wCorr)
wcor <- function(x,y,w,nrep=100) {
    obs=weightedCorr(x,y,method="Spearman",weights=w)
    perm <- replicate(nrep, weightedCorr(x,sample(y),method="Spearman",weights=w))
    z <- (obs-mean(perm))/sd(perm)
    pnorm(abs(z),lower.tail=FALSE) * 2
}
    
## consistency
fcons <- function(dt,res,label) {
    dt2 <- merge(dt[!is.na(beta)],
                 rotuse[,.(pid,centers,rot,PC)],
                 by="pid",allow.cartesian=TRUE)
    dt2[,betas:=beta*shrinkage]
    ## summarise
    SUMM <- dt2[,p.wcors:=wcor(rot,betas,shrinkage), by=c("PC","trait")]
    SUMM <- unique(SUMM[,.(PC,trait,p.wcors)])
    ## add in pc sig
    SUMM <- merge(SUMM,unique(res[,.(PC,trait,p.value,newfdr)]),by=c("PC","trait"))
    SUMM[,cutfdr:=cut(newfdr,c(0,0.001,0.01,0.1,1))]
    SUMM[,cat:=label]
    copy(SUMM)
}

## plotting
library(cowplot); theme_set(theme_cowplot())
plotsumm <- function(...,var="p.wcors") {
    tmp <- rbind(...,fill=TRUE)
    tmp[,cutfdr:=cut(newfdr,c(0,0.01,0.1,0.5,1))]
    levels(tmp$cutfdr)
    levels(tmp$cutfdr) <- c("<0.01","<0.1","<0.5","<1")
    tmp[,logp:=-log10(tmp[[var]])]
    ## tmp[,logp:=(tmp[[var]])]
    ggplot(tmp[!is.na(newfdr)], aes(x=cutfdr,
                    y=logp,
                    fill=cat)) +
        geom_boxplot(notch=TRUE,coef=1, pch=15) +
        facet_grid(.~cat) +
        background_grid("y","none") +
        xlab("FDR category") +
        ylab("-log10 p") +
        scale_fill_discrete("Trait group") +
        theme(legend.position="none")
}

################################################################################
library(Matrix)
mymelt <- function(m,row.name="rn",...) {
    if(!is.matrix(m))
        return(melt(m,...))
    dt <- data.table(rn=rownames(m),m)  %>% melt(.,id="rn",...)
    if(row.name!="rn")
        setnames(dt,"rn",row.name)
    dt
}

                        ## shrinkage and rotation info
use <- mymelt(use.pca,row.name="pid",variable.name="PC",value.name="use")
rot <- mymelt(rot.pca,row.name="pid",variable.name="PC",value.name="rot")
rotuse <- merge(rot,use,by=c("pid","PC"))
rotuse <- merge(rotuse,shrink.DT[,c("pid",'shrinkage'),with=FALSE],by="pid")
## thin - max 1 snp per MB
rotuse <- rotuse[use==TRUE,]
rotuse[,c("chr","pos"):=tstrsplit(pid,":")  %>% lapply(.,as.numeric)]
rotuse <- rotuse[order(PC,chr,pos),]
## find LD-dependent SNPs
f <- function(dt) {
    dt <- dt[order(abs(rot),decreasing=TRUE)]    
    ldt <- LD[dt$pid,dt$pid]
    diag(ldt) <- 0
    ldt[lower.tri(ldt)] <- 0
    dt$maxld <- apply(abs(ldt),2,max)
    dt[maxld<0.1] # ie rsq < 0.01
}
cat("total sparse driver snps")
table(rotuse$PC)
rotuse %<>% split(.,.$PC) %>% lapply(.,f)  %>% rbindlist()
cat("nearly-independent sparse driver snps")
table(rotuse$PC)
rotuse[,centers:=beta.centers[pid]]

################################################################################
## load UKBB data
res1 <- proj[category=="UKBB"]

## significant traits
wanted1 <- unique(res1[newfdr<0.01,]$trait) #   %>% sub("SRD.","",.)
## sample same number of non-significant traits
set.seed(42) # for reproducible plots
wanted0 <- sample(unique(res1[newfdr>0.01,]$trait),length(wanted1)) #%>% sub("SRD.","",.) 
wanted <- unique(c(wanted1,wanted0))
length(wanted)

DT1=readraw(wanted,pids=rownames(use.pca))
SUMM1 <- fcons(DT1,res1,"UKBB")

################################################################################

## all the other projected GWAS
res2 <- proj[fdrcat=="general" & category!="UKBB"]
wanted <- res2$trait  %>% unique()
wanted
DT2=readraw(wanted,pids=rownames(use.pca))
SUMM2 <- fcons(DT2,res2,"GWAS")

################################################################################

## bloods
## wanted <- res3$trait  %>% unique()
## wanted
wanted <- c('pdw','mpv','plt','irf','ret','rdw','hct','mch','mono','baso','eo','neut','lymph')  %>% paste0("ASTLE:",.)
DT3=readraw(wanted,pids=rownames(use.pca))
sdt <- split(DT3,DT3$trait)
res3 <- lapply(sdt, function(x) project.sparse(x$beta, x$seb, x$pid))
for(nm in names(sdt))
    res3[[nm]][,trait:=nm]
res3 %<>% rbindlist()
setnames(res3,"p","p.value")
res3[,newfdr:=p.adjust(p.value,method="BH"),by="PC"]

SUMM3 <- fcons(DT3,res3,"Blood counts")

plotsumm(SUMM1,SUMM2,SUMM3)
SUMM3[p.wcors<1e-3]

################################################################################

## cytokines
res4 <- proj[grepl("Cytokine",category.label)]
wanted <- res4$trait  %>% unique()
wanted
DT4=readraw(wanted, pids=rownames(use.pca))
SUMM4 <- fcons(DT4,res4,"Cytokines")

################################################################################

## Roederer
traits <- list.files(DATA_DIR,pattern="roederer.*_source.RDS")  %>%
  sub("_source.RDS","",.)
DT5=readraw(traits,pids=rownames(use.pca))
sdt <- split(DT5,DT5$trait)
res5 <- lapply(sdt, function(x) project.sparse(x$beta, x$seb, x$pid))
for(nm in names(sdt))
    res5[[nm]][,trait:=nm]
res5 %<>% rbindlist()
setnames(res5,"p","p.value")
res5[,newfdr:=p.adjust(p.value,method="BH"),by="PC"]
wanted <- res5$trait  %>% unique()
wanted
DT5=readraw(wanted, pids=rownames(use.pca))
SUMM5 <- fcons(DT5,res5,"Flow")

################################################################################
## pqtl?

################################################################################
## plot
plotsumm(SUMM1,SUMM2,SUMM3,SUMM4,SUMM5,var="p.wcors") 
ggsave("~/share/as_basis/figures/suppfig-consistency.pdf",height=6,width=8)

