## install.packages("gtx")  ## summary statistic GRS from Toby Johnson
library(data.table)
library(magrittr)
library(gtx)

## summary basis
source("R/cw-reader.R") # also loads utils and files
proj=reader()

## load source data
todo <- proj[PC %in% c("PC3","PC13") &
             ((fdr.overall<0.01 & newfdr<0.01) # significant
                 | trait %in% c("SRD.osteoarthritis","SRD.rheumatoid.arthritis","SRD.type.1.diabetes") # comp to external data
                 | grepl("jia",trait) # comp to external data
             ) &
             (fdrcat %in% c("general","Cytokines") | trait=="ASTLE:eo")]
data <- read_raw(unique(todo$trait))

## sparse drivers
(load(SPARSE_BASIS_FILE)) # use.pca, stats
## source("~/Projects/auto-basis/scripts/annot-snps.R")
## source("~/Projects/auto-basis/scripts/plots.R")
pids.pc3 <- rownames(use.pca)[use.pca[,"PC3"]==TRUE]
pids.pc13 <- rownames(use.pca)[use.pca[,"PC13"]==TRUE]
pids.basis <- rownames(use.pca)

## summary data for IP10/MIG/eo
tmp <- data[!is.na(p.value) & seb!=0 & !is.na(beta) & !is.na(seb)]
S <- split(tmp, tmp$trait)
P <- S[c("CK:IP10","CK:MIG","ASTLE:eo")]
T <- S[setdiff(names(S),c("CK:IP10","CK:MIG","ASTLE:eo"))]

## build 5 inputs for a P-T pair:
## PC snps
## p < 1e-3 for P
## p < 1e-4 for P
## p < 1e-3 for T
## p < 1e-4 for T

man.DT <- readRDS(SNP_MANIFEST_FILE)
man.DT[,chr:=sub(":.*","",pid)]

pnm <- "ASTLE:eo"; tnm <- "mpo.anca.positive.egpa_lmm"
library(MendelianRandomization)
library(Matrix)
make_inp <- function(pids.use, LD, tnm, pnm, RTHR=sqrt(0.8)) {
    xdata <- S[[tnm]][order(p.value),]
    ydata <- S[[pnm]]
    myld <- LD[pids.use,pids.use]
    ## print(dim(myld))
    pids.use <- Reduce(intersect,list(pids.use,ydata$pid,xdata$pid))
    ivld <- indep_snps(xdata[pid %in% pids.use],
                       myld,
                       rthr=RTHR)
    ## setkey(ivld,pid)
    if(length(pids.use)==1)
        return(NULL)
    inp=mr_input(bx=xdata[match(ivld$pid,pid)]$beta,
                 bxse=xdata[match(ivld$pid,pid)]$seb,
                 by=ydata[match(ivld$pid,pid)]$beta,
                 byse=ydata[match(ivld$pid,pid)]$seb,
                 correlation=as.matrix(myld[ivld$pid,ivld$pid]),
                 exposure=pnm,
                 outcome=tnm,
                 snps=ivld$pid)
}

library(parallel)
builder <- function(tnm, pnm, RTHR=sqrt(0.8)) {
    message(tnm, " - vs - ", pnm)
    ## which pc
    pc <- if(pnm=="ASTLE:eo") { "pc13" } else { "pc3" }
    pids <- list(pc=if(pnm=="ASTLE:eo") { pids.pc13 } else { pids.pc3 },
                 p3 = S[[pnm]][p.value<1e-3]$pid,
                 p4 = S[[pnm]][p.value<1e-4]$pid,
                 t3 = S[[tnm]][p.value<1e-3]$pid,
                 t4 = S[[tnm]][p.value<1e-4]$pid)
    pids <- pids[ sapply(pids,length) >= 2 ]
    allpids <- unlist(pids) %>% unique()
    
    ## get LD to allow thinning to make GRS
    tmp <- man.DT[pid %in% allpids]
    s.DT <- split(tmp, tmp$chr)
    BIGLD <- lapply(names(s.DT), function(chr) {
        ss.file <- file.path(REF_GT_DIR, sprintf("%s.RDS", chr))
        sm <- readRDS(ss.file)
        pids <- colnames(sm)
        dup.idx <- which(duplicated(pids))
        if (length(dup.idx) > 0) {
            sm <- sm[, -dup.idx]
            pids <- pids[-dup.idx]
        }
        sm.map <- match(s.DT[[chr]]$pid, pids)
        r <- ld(sm[, sm.map], sm[, sm.map], stats = "R")
        r[is.na(r)] <- 0
        r
    })  %>% bdiag_with_dimnames(.)
    ## print(dim(BIGLD))

    ## print(table(allpids %in% rownames(BIGLD)))

    ## print(str(pids))
    ## inp <- vector("list",length(pids))
    ## names(inp) <- names(pids)
    ## for(nm in names(pids)) {
    ##     message(nm)
    ##     inp[[nm]] <- make_inp( pids[[nm]], BIGLD, tnm, pnm)
    ##     }

    ## return(inp)

    inp <- mclapply(pids, make_inp, BIGLD, tnm, pnm, RTHR=RTHR, mc.cores=length(pids))
}

traits.pc3 <- todo[PC=="PC3" & fdr.overall<0.01 & newfdr< 0.01 &
                   !(trait %in% c("ASTLE:eo","CK:IP10","CK:MIG"))]$trait
traits.pc13 <- todo[PC=="PC13" & fdr.overall<0.01 & newfdr< 0.01 &
                    !(trait %in% c("ASTLE:eo","CK:IP10","CK:MIG"))]$trait

INP.eo <- lapply(traits.pc13, builder, pnm="ASTLE:eo")
names(INP.eo) <- paste("ASTLE:eo",traits.pc13,sep="_")
INP.ip10 <- lapply(traits.pc3, builder, pnm="CK:IP10")
names(INP.ip10) <-  paste("CK:IP10",traits.pc3, sep="_")
INP.mig <- lapply(traits.pc3, builder, pnm="CK:MIG")
names(INP.mig) <- paste("CK:MIG",traits.pc3,sep="_")

INP <- c(INP.eo, INP.ip10, INP.mig)

runmr <- function(x) {
    if(is.list(x)) {
        ret <- lapply(x, runmr)  %>% do.call("rbind",.)
        ret$limit <- names(x)
        return(ret)
    }
    ivw <- mr_ivw(x)
    egg <- mr_egger(x)
    data.frame(exposure=ivw@Exposure,outcome=ivw@Outcome,
               ivw.est=ivw@Estimate,
               ivw.se=ivw@StdError,
               ivw.p=ivw@Pvalue,
               nsnps=ivw@SNPs,
               egger.est=egg@Estimate,
               egger.se=egg@StdError.Est,
               egger.p=egg@Pvalue.Est,
               int.est=egg@Intercept,
               int.se=egg@StdError.Int,
               int.p=egg@Pvalue.Int)
}

MR <- lapply(INP, runmr)  #%>% do.call("rbind",.)
MR[[1]]



plotraw_df=function(object) {
    data.table(snp=object@snps,
               Bx = object@betaX,
               By = object@betaY,
               Bxse = object@betaXse,
               Byse = object@betaYse)
}

geom_cross <- function(...) {
    varnames=lapply(substitute(list(...))[-1], deparse)
    L <- list(...)
    ## L[["alpha"]] <- 0.5
    Lp <- Ll <- L
    Ll[["size"]] <- NULL
    Ll[["shape"]] <- NULL
    list(do.call(geom_point,Lp),
         do.call(geom_errorbar,Ll),
         do.call(geom_errorbarh,Ll))
}

c4="#ff553f" # exposure
c2="#EE442f" # exposure
c1 <- "#63acbe" # outcome
c3 <- "grey10"

plotraw <- function(mres) {
    ## datasets
    nm <- paste(mres$exposure,mres$outcome,sep="_")[1]
    
    df <- lapply(INP[[nm]], plotraw_df)
    df$t3 <- df$t3[ !(snp %in% df$t4$snp) ]
    df$p3 <- df$p3[ !(snp %in% df$p4$snp) ]
    ## df <- do.call("rbind",df2)
    ## df[,basis:=df$snp %in% INP[[nm]]$pc@snps]
    ## df[,t4:=df$snp %in% INP[[nm]]$pc@snps]

    ## df1 <- plotraw_df(o1)
    ## df2 <- plotraw_df(o2)
    ## df <- rbind(df1[snp %in% intersect(df1$snp,df2$snp),][,source:="both"],
    ##             df1[snp %in% setdiff(df1$snp,df2$snp),][,source:="exposure"],
    ##             df2[snp %in% setdiff(df2$snp,df1$snp),][,source:="PC"])
    
    ## Create the initial plot
    p <- ggplot(data=NULL, aes(x = Bx, y = By,
                               ymin = By - qnorm(0.975)*Byse,
                               ymax = By + qnorm(0.975)*Byse,
                               xmin = Bx - qnorm(0.975)*Bxse,
                               xmax = Bx + qnorm(0.975)*Bxse)) +
      geom_cross(size = 1,alpha=0.5,data=df$p3,colour=c1,shape=5) +
      geom_cross(size = 3,alpha=0.5,data=df$p4,colour=c1,shape=5) +
      geom_cross(size = 1,alpha=0.5,data=df$t3,colour=c2,shape=4) +
      geom_cross(size = 3,alpha=0.5,data=df$t4,colour=c2,shape=4) +
      geom_cross(size = 3,alpha=0.5,data=df$pc,colour=c3,shape=15) +
      geom_abline(intercept = 0, slope = mr_ivw(INP[[nm]]$p3)$Estimate, color = c1,linetype=2) +
      geom_abline(intercept = 0, slope = mr_ivw(INP[[nm]]$p4)$Estimate, color = c1,linetype=1) +
      geom_abline(intercept = 0, slope = mr_ivw(INP[[nm]]$t3)$Estimate, color = c2,linetype=2) +
      geom_abline(intercept = 0, slope = mr_ivw(INP[[nm]]$t4)$Estimate, color = c2,linetype=1) +
      geom_abline(intercept = 0, slope = mr_ivw(INP[[nm]]$pc)$Estimate, color = c3)  +
      labs(x=mres$exposure[1], y=mres$outcome[1]) 
    p
}


plotraw(MR[[4]])

################################################################################

## overall results plot

MRdt <- rbindlist(MR)
MRdt  %<>%  merge(., unique(todo[,.(trait, trait.label, category.label)]),
                  by.x="outcome",by.y="trait")
MRdt[,limit:=factor(limit, levels=c("t3","t4","p3","p4","pc"))]


c4="#ff775f" # exposure
c2="#cc220f" # exposure
c1 <- "#63acbe" # outcome
c3 <- "grey10"
what="est"; x=MRdt[exposure=="ASTLE:eo"];limits.do=c("p3","p4","pc")

library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
plotter <- function(x,what=c("est","int","grs"),
                    limits.do=c("pc","p3","t3")) {
    what <- match.arg(what)

    x[,x:=paste(trait.label,limit,category.label)]
    ox <- with(x[limit=="pc"], x[order(ivw.est)])  %>%
      lapply(., function(s) sapply(sort(limits.do), function(l) sub("pc",l,s)))  %>%
      unlist()  %>%
      unique()
    x[,x:=factor(x,levels=ox)]
    p <- switch(what,
        est=ggplot(x[limit %in% limits.do],
                    aes(x=x,y=ivw.est,ymin=ivw.est-1.96*ivw.se,ymax=ivw.est+1.96*ivw.se,
                        shape=limit,colour=limit)),
        int=ggplot(x[limit %in% limits.do],
                    aes(x=limit,y=int.est,ymin=int.est-1.96*int.se,ymax=int.est+1.96*int.se,
                        shape=limit,colour=limit)),
        grs=ggplot(x[limit %in% limits.do],
                    aes(x=limit,y=ahat,ymin=ahat-1.96*aSE,ymax=ahat+1.96*aSE,
                        shape=limit,colour=limit))
        )
    p <- p + geom_vline(xintercept=seq(0.5,length(ox)+0.5,by=3),
                        size=0.5,colour="grey85") +
      geom_point() +
      geom_linerange() + 
      geom_hline(yintercept=0,col="grey",linetype=2)
    if(length(unique(x$exposure))>1) {
        p <- p + facet_grid(category.label ~ exposure, scales="free",switch="y",space="free") 
          } else {
              p <- p + facet_grid(category.label ~ ., scales="free",switch="y",space="free")
          }
  p <- p + coord_flip() +
  background_grid(major="none") +
    scale_shape_manual("SNP selection",
                       values=c(p3=5,p4=4,t3=4,t5=4,pc=15),
                       labels=function(s) c(p3="p<10-3",p4="p<10-4",pc="PC driver SNPs")[s]) +
    scale_colour_manual("SNP selection",
                        values=c(p3=c4,p4=c2,t3=c1,t5=c1,pc=c3),
                        labels=function(s) c(p3="p<10-3",p4="p<10-4",pc="PC driver SNPs")[s]) +
    scale_x_discrete(labels=function(a) ifelse(grepl(median(limits.do),a),
                                               sub(paste0(" ",median(limits.do),".*"),"",a),
                                               "")) +
    scale_y_continuous("Estimated causal effect") +
  theme(plot.title=element_text(hjust=0),
        strip.placement = "outside",
        strip.background=element_rect(fill="grey90"),
        ## axis.ticks.x=element_blank(),
        ## axis.text.x=element_blank(),
        legend.position="top",
        axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        strip.text.y = element_text(angle=180))
    p


}

## plotter(MRdt[exposure=="ASTLE:eo", ])
## plotter(MRdt[exposure=="CK:IP10", ])
## plotter(MRdt[exposure=="CK:MIG", ])

plotter(x=MRdt[exposure=="ASTLE:eo"],limits.do=c("p3","p4","pc"))
ggsave("figures/suppfig-mr-astle.pdf",height=8,width=6,scale=1.5)
## plotter(x=MRdt[exposure=="CK:IP10"],limits.do=c("p3","p4","pc"))
## plotter(x=MRdt[exposure=="CK:MIG"],limits.do=c("p3","p4","pc"))
plotter(x=MRdt[exposure %in% c("CK:IP10","CK:MIG")],limits.do=c("p3","p4","pc"))
ggsave("figures/suppfig-mr-cytokines.pdf",height=8,width=8,scale=1.5)




################################################################################

## GRS

make_score <- function(ivdata) {
    if(is.list(ivdata)) {
        ret <- lapply(ivdata, make_score)  %>% do.call("rbind",.)
        ret$limit <- names(ivdata)
        return(ret)
    }
    res <- grs.summary(w=ivdata@betaX,
                b=ivdata@betaY,
                s=ivdata@betaYse,
                n=5000)   %>% do.call(data.frame,.)
    res$exposure=ivdata@exposure
    res$outcome=ivdata@outcome
    res
}



THINP.eo <- lapply(traits.pc13, builder, pnm="ASTLE:eo", RTHR=sqrt(0.2))
names(THINP.eo) <- paste("ASTLE:eo",traits.pc13,sep="_")
THINP.ip10 <- lapply(traits.pc3, builder, pnm="CK:IP10")
names(THINP.ip10) <-  paste("CK:IP10",traits.pc3, sep="_")
THINP.mig <- lapply(traits.pc3, builder, pnm="CK:MIG")
names(THINP.mig) <- paste("CK:MIG",traits.pc3,sep="_")

THINP <- c(THINP.eo, THINP.ip10, THINP.mig)



scores <- lapply(THINP, make_score)  %>% rbindlist()
scores  %<>%  merge(., unique(todo[,.(trait, trait.label, category.label)]),
                  by.x="outcome",by.y="trait")
scores[,limit:=factor(limit, levels=c("t3","t4","p3","p4","pc"))]

plotter(scores,what="grs",limits.do=c("p3","p4","pc"))


## experimental data on IP10

plasma.ip10 <- data.table(trait=c("jia_PO_19","jia_RFneg_19","jia_sys_19","SRD.rheumatoid.arthritis","SRD.type.1.diabetes","control"),
                     ahat=c(719,	388,	106,	383,		242,	98),
                     sd=c(1016,421, 108, 279, 105, 36),
                     n=c(30,20,15,9,9,20))
plasma.mig <- data.table(trait=c("jia_PO_19","jia_RFneg_19","jia_sys_19","SRD.rheumatoid.arthritis","SRD.type.1.diabetes","control"),
                     ahat=c(111,	73,	23,	26,		6.3,8.4),
                     sd=c(250,150,53,26,10,3),
                     n=c(30,20,15,9,9,20))
sf.ip10 <- data.table(trait=c("jia_PO_19","jia_RFneg_19","jia_sys_19","SRD.rheumatoid.arthritis","SRD.osteoarthritis"),
                 ahat=c(4001,	5128,	2772,	4289,	197),
                 sd=c(959,	  1759, 1443,	 3149,231),
                 n=c(19,9,5,10,5))
sf.mig <- data.table(trait=c("jia_PO_19","jia_RFneg_19","jia_sys_19","SRD.rheumatoid.arthritis","SRD.osteoarthritis"),
                 ahat=c(2817,2266,85,673,13),
                 sd=c(844,123,19,973,7.8),
                 n=c(19,9,5,10,5))
edata <- rbind(plasma.ip10[,c("cytokine","fluid"):=list("IP10","plasma")],
               sf.ip10[,c("cytokine","fluid"):=list("IP10","sf")],
               plasma.mig[,c("cytokine","fluid"):=list("MIG","plasma")],
               sf.mig[,c("cytokine","fluid"):=list("MIG","sf")])
edata[,aSE:=sd/sqrt(n)]

use.traits <- unique(edata$trait)
table(use.traits %in% MRdt$trait)



### HERE


pc3 <- x$proj[trait %in% grs.pred$trait & variable=="PC3",.(trait,ahat=-delta,aSE=sqrt(variance))]

pc3[,what:="PC3"]
plasma[,what:="plasma.pg.ml"]
sf[,what:="synovial_fluid.pg.ml"]






df=rbindlist(list(pc3[trait!="MIG",],plasma[trait!="control"],sf,grs.pred),fill=TRUE)[trait %in% use.traits]
o <- order(plasma$ahat)
df[,trait:=factor(trait,levels=c("SRD.osteoarthritis",plasma$trait[o]))]
levels(df$trait)  %<>% sub("jia_sys_19","jia_systemic",.)  %>%
  sub("SRD.osteoarthritis","osteoarthritis",.) %>%
  sub("SRD.type.1.diabetes","type_1_diabetes",.)  %>%
  sub("SRD.rheumatoid.arthritis","rheumatoid_arthritis",.)  %>%
  sub("jia_PO_19","jia_oligoarticular",.)  %>%
  sub("jia_RFneg_19","jia_polyarticular",.)  %>% 
  sub("jia_RFpos_19","jia_polyarticular",.)


p.grs <- ggplot(df,aes(x=trait,y=ahat,ymin=ahat-1.96*aSE,ymax=ahat+1.96*aSE)) +
  geom_pointrange() +
  facet_wrap(~what,scales="free_y",ncol=1) +
  background_grid() +
  geom_hline(aes(yintercept=0)) +
  theme(axis.text.x=element_text(angle=45,hjust=1))

ggsave("~/IP10.pdf",plot=p.grs,height=6,width=8)

## compare weights and IP10 gwas coefficients
SHRINKAGE_FILE <- '/home/ob219/share/as_basis/GWAS/support/ss_shrinkage_gwas_13_traits_0919.RDS'
  shrink.DT <- readRDS(SHRINKAGE_FILE)
setnames(shrink.DT,'ws_emp_shrinkage','shrinkage')
tmp <- ip10[pid %in% rownames(rot.pca)[use.pca[,"PC3"]==TRUE]]
tmp[,pc3.coef:=rot.pca[pid,"PC3"]]
tmp <- merge(tmp,shrink.DT[,.(pid,shrink=shrinkage)],by="pid")
## tmp[,shrink:=as.numeric(shrink)]
p.coef <- ggplot(tmp,aes(x=-log(or),y=shrink*pc3.coef)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(x="IP10 GWAS coefficient", y="PC3 coefficient") +
  background_grid()

plot_grid(p.pc3 + theme(axis.text.x=element_blank()) + ylab("PC3") + scale_y_reverse(),
          plot_grid(p.coef,p.grs,nrow=2,rel_heights=c(.3,.7)),
          rel_widths=c(.6,.4))
ggsave("~/IP10.pdf",height=10,width=12)


with(tmp,cor(-log(or),pc3.coef))

0

    ## scale_y_log10()

## A named list with the following elements:
## m is the number of SNPs used in the risk score.
## n is the input sample size.
## X2m is the chi squared test statistic for an m d.f. test in the testing dataset (all SNPs have independent effects).
## R2m is the (pseudo) variance explained by the m d.f. model in the testing dataset.
## ahat is the estimated coefficent for regressing the response onto the m SNP risk score.
## aSE is the standard error.
## X2rs is the chi squared test statistic for a 1 d.f. test for the risk score in the testing dataset.
## R2rs is the (pseudo) variance explained by the risk score model in the testing dataset.
## pval is the P-value for the 1 d.f. test.
## Qrs is the (m-1) d.f. heterogeneity test statistic.
## phet is the heterogeneity test P-value.

## all predictive, except hypothyroidism and crohn's, all heterogenous

## do same with PC4 driver snps only

library(ggplot2)
m <- merge(scores.basis,scores.pc4,by="nm",suffixes=c(".basis",".pc4"))
m <- merge(m,scores.any,by="nm")
p1=ggplot(m,aes(x=-log10(pval),y=-log10(pval.pc4))) + geom_point() + geom_abline()
p2=ggplot(m,aes(x=-log10(pval),y=-log10(pval.basis))) + geom_point() + geom_abline()
p3=ggplot(m,aes(x=-log10(pval.basis),y=-log10(pval.pc4))) + geom_point() + geom_abline()

library(cowplot)
plot_grid(p1,p2,p3)

m <- merge(scores.any,scores.pc4,by="nm",suffixes=c(".any",".pc4"))
ggplot(m,aes(x=-log10(pval.any),y=-log10(pval.pc4))) + geom_point() + geom_abline()

## check input

m=merge(data[trait=="CK.IP10",.(pid,beta,p.value,trait)],
        data[trait!="CK.IP10",.(pid,beta,p.value,trait)],
        by=c("pid"),suffixes=c(".ip10",""))

ggplot(m[p.value<0.01 | p.value.ip10<0.01],aes(x=beta.ip10,y=beta,col=pid %in% pids.pc4)) + geom_point() + geom_abline() + geom_smooth(se=FALSE,method="lm") + facet_wrap(~trait,scales="free")

ggplot(m[p.value<0.01 | p.value.ip10<0.01],aes(x=-log10(p.value.ip10),y=-log10(p.value),col=pid %in% pids.pc4)) + geom_point() + geom_abline() + geom_smooth(se=FALSE,method="lm") + facet_wrap(~trait,scales="free")

ggplot(m[pid %in% pids.pc4 & (p.value<0.001 | p.value.ip10<0.001)],aes(x=beta.ip10,y=beta,col=trait)) + geom_point() + geom_abline() + geom_smooth(se=FALSE,method="lm") + facet_wrap(~trait)
