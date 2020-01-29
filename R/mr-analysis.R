## install.packages("gtx")  ## summary statistic GRS from Toby Johnson
library(data.table)
library(magrittr)
library(gtx)

## summary basis
source("R/cw-reader.R") # also loads utils and files
source("R/cw-palette.R")
proj=reader()

## load source data
todo <- proj[PC %in% c("PC3","PC13") &
             ((fdr.overall<0.01 & newfdr<0.01) # significant
                 | trait %in% c("SRD.osteoarthritis","SRD.rheumatoid.arthritis","SRD.type.1.diabetes") # comp to external data
                 | grepl("jia",trait) # comp to external data
             ) &
             (fdrcat %in% c("general","Cytokines") | trait=="ASTLE:eo")]
data <- read_raw(c(unique(todo$trait),
                   c("UKBB_NEALE:SRD:type.1.diabetes",    
                     "UKBB_NEALE:SRD:osteoarthritis")))

## sparse drivers
(load(SPARSE_BASIS_FILE)) # use.pca, stats
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
    xdata <- S[[pnm]][order(p.value),] # exposure / predictor
    ydata <- S[[tnm]] # outcome / trait
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
builder <- function(tnm, pnm, RTHR=sqrt(0.01)) {
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

    inp <- mclapply(pids, make_inp, BIGLD, tnm=tnm, pnm=pnm,
                    RTHR=RTHR, mc.cores=length(pids))
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

runmr <- function(x,reverse=FALSE) {
    if(is.list(x)) {
        ret <- lapply(x, runmr, reverse=reverse)  %>% do.call("rbind",.)
        ret$limit <- names(x)
        return(ret)
    }
    message(x@exposure,"\t",x@outcome)
    if(reverse) {
        x <- make_inp(pids.use=x@snps, LD=x@correlation, tnm=x@exposure, pnm=x@outcome)
    }
    ivw <- mr_ivw(x)
    ## egg <- mr_egger(x2)
    data.frame(exposure=ivw@Exposure,outcome=ivw@Outcome,
               ivw.est=ivw@Estimate,
               ivw.se=ivw@StdError,
               ivw.p=ivw@Pvalue,
               nsnps=ivw@SNPs)
               ## egger.est=egg@Estimate,
               ## egger.se=egg@StdError.Est,
               ## egger.p=egg@Pvalue.Est,
               ## int.est=egg@Intercept,
               ## int.se=egg@StdError.Int,
               ## int.p=egg@Pvalue.Int)
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
      ## geom_cross(size = 1,alpha=0.5,data=df$p3,colour=c1,shape=5) +
      geom_cross(size = 3,alpha=0.3,data=df$pc,colour=c3,shape=15) +
      geom_cross(size = 3,alpha=0.5,data=df$p4,colour=c1,shape=5) +
      ## geom_cross(size = 1,alpha=0.5,data=df$t3,colour=c2,shape=4) +
      geom_cross(size = 3,alpha=0.5,data=df$t4,colour=c2,shape=4) +
      ## geom_abline(intercept = 0, slope = mr_ivw(INP[[nm]]$p3)$Estimate, color = c1,linetype=2) +
      geom_abline(intercept = 0, slope = mr_ivw(INP[[nm]]$p4)$Estimate, color = c1,linetype=1) +
      ## geom_abline(intercept = 0, slope = mr_ivw(INP[[nm]]$t3)$Estimate, color = c2,linetype=2) +
      geom_abline(intercept = 0, slope = mr_ivw(INP[[nm]]$t4)$Estimate, color = c2,linetype=1) +
      geom_abline(intercept = 0, slope = mr_ivw(INP[[nm]]$pc)$Estimate, color = c3)  +
      labs(x=mres$exposure[1], y=mres$outcome[1]) 
    p
}

library(ggplot2)
plotraw(MR[[4]])
plotraw(MR[["CK:IP10_UKBB_NEALE:SRD:hyperthyroidism.thyrotoxicosis" ]])
plotraw(MR[[ "ASTLE:eo_UKBB_NEALE:SRD:asthma"]])
plotraw(MR[[ "ASTLE:eo_UKBB_NEALE:SRD:crohns.disease"]])

MREV <- lapply(INP[c("CK:IP10_UKBB_NEALE:SRD:hyperthyroidism.thyrotoxicosis",
                     "ASTLE:eo_UKBB_NEALE:SRD:asthma")], runmr, reverse=TRUE)
MREV
MR[names(MREV)]
################################################################################

## overall results plot

MRdt <- rbindlist(MR)
MRdt  %<>%  merge(., unique(todo[,.(trait, trait.label, category.label)]),
                  by.x="outcome",by.y="trait")
MRdt[,limit:=factor(limit, levels=c("t3","t4","p3","p4","pc"))]
MRdt[category.label=="Myasthenia gravis", category.label:="Myasthenia\ngravis"]


c4=mygreen #"#ff775f" # exposure
c2=mygreen #"#cc220f" # exposure
c1 <- "#63acbe" # outcome
c3 <- myblue #"OliveDrab" #"grey10"
what="est"; x=MRdt[exposure=="ASTLE:eo"];limits.do=c("p3","p4","pc")

library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
plotter <- function(x,what=c("est","int","grs"),
                    limits.do=c("pc","p3","t3")) {
    what <- match.arg(what)

    ynm <- switch(what, est="ivw.est",int="int.est",grs="ahat")
    x[,x:=paste(trait.label,limit,category.label)]
    ## ox <- with(x[limit=="pc"], x[order(x[[ynm]])])  %>%
    ox <- with(x[limit=="pc"], x[order(category.label,ivw.est)])  %>%
      lapply(., function(s) sapply(sort(limits.do), function(l) sub("pc",l,s)))  %>% ## add other limits
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
    p <- p + geom_vline(xintercept=seq(0.5,length(ox)+0.5,by=length(limits.do)),
                        size=0.5,colour="grey85") +
      geom_point(size=3) +
      geom_linerange() + 
      geom_hline(yintercept=0,col="grey",linetype=2)
    if(length(unique(x$exposure))>1) {
        p <- p + facet_grid(category.label ~ exposure, scales="free",switch="y",space="free_y") 
          } else {
              p <- p + facet_grid(category.label ~ ., scales="free",switch="y",space="free_y")
          }
  p <- p + coord_flip() +
  background_grid(major="none") +
    scale_shape_manual("SNP selection",
                       values=c(p3=5,p4=4,t3=4,t5=4,pc=15),
                       labels=function(s) c(p3="Exposure p<10-3",p4="Exposure p<10-4",pc="PC driver SNPs")[s]) +
    scale_colour_manual("SNP selection",
                        values=c(p3=mygreen,p4=c2,t3=c2,t5=c2,pc=c3),
                        labels=function(s) c(p3="Exposure p<10-3",p4="Exposure p<10-4",pc="PC driver SNPs")[s]) +
    scale_x_discrete(labels=function(a) ifelse(grepl(limits.do[2],a),
                                               sub(paste0(" ",paste(limits.do,collapse="|"),".*"),"",a) ,
                                               "")) +
    scale_y_continuous("Estimated causal effect") +
  theme(plot.title=element_text(hjust=0),
        strip.placement = "outside",
        strip.background=element_rect(fill="grey90"),
        ## axis.ticks.x=element_blank(),
        ## axis.text.x=element_blank(),
        legend.position="top",
        axis.line=element_blank(),
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
ggsave("figures/slides-mr-astle.pdf",height=6,width=8,scale=1.5)
## plotter(x=MRdt[exposure=="CK:IP10"],limits.do=c("p3","p4","pc"))
## plotter(x=MRdt[exposure=="CK:MIG"],limits.do=c("p3","p4","pc"))
plotter(x=MRdt[exposure %in% c("CK:IP10","CK:MIG")],limits.do=c("p3","p4","pc"))
ggsave("figures/suppfig-mr-cytokines.pdf",height=8,width=8,scale=1.5)


## add forest
theme_set(theme_cowplot(10))
(load("~/basis-pc3-forest.RData"))

fdata <- as.data.table(p.pc3$data)[newfdr<0.01][,category.label:=as.character(category.label)]
tmp <- rbind(MRdt[exposure %in% c("CK:IP10","CK:MIG"),
                  .(exposure=sub("CK:","MR exposure=",exposure),
                    category.label,trait.label,outcome, ivw.est,ivw.se,limit)],
             fdata[,.(category.label,trait.label,outcome=trait,exposure=PC,
                      ivw.est=-delta,ivw.se=sqrt(var.proj),limit="pc")],
             data.table(category.label="Cytokines",trait.label=c("IP10","MIG"),
                        outcome=c("CK:IP10","CK:MIG"),limit="p3",exposure="PC3"),
             fill=TRUE)
tmp$exposure %<>% as.factor()  %>% relevel(.,"PC3")
tmp[,category.label:=factor(category.label,
                            levels=c("Cytokines","JIA","Myasthenia\ngravis",
                                     "Ank.Spond","PsA", "UKBB"))]
plotter(x=tmp, limits.do=c("p3","pc")) +
  scale_y_continuous(breaks=0) +
  theme(legend.position="bottom",legend.justification="right",
        axis.title.x=element_blank())
ggsave("figures/fig6-pc3-mr.pdf",height=6,width=8)
ggsave("figures/slides-pc3-mr.pdf",height=6,width=12,scale=1,pointsize=16)


(load("~/basis-pc13-forest.RData"))

fdata <- as.data.table(p.pc13$data)[newfdr<0.01][,category.label:=as.character(category.label)]
tmp <- rbind(MRdt[exposure %in% c("ASTLE:eo"),
                  .(exposure=sub("ASTLE:eo","MR exposure=eos count",exposure),
                    category.label,trait.label,outcome, ivw.est,ivw.se,limit)],
             fdata[,.(category.label,trait.label,outcome=trait,exposure=PC,
                      ivw.est=delta,ivw.se=sqrt(var.proj),limit="pc")],
             data.table(category.label="Blood counts",trait.label=c("eosinophils"),
                        outcome=c("ASTLE:eo"),limit="p3",exposure="PC13"),
             fill=TRUE)
tmp$exposure %<>% as.factor()  %>% relevel(.,"PC13")
tmp[,category.label:=factor(category.label,
                            levels=c("Blood counts", "EGPA","Myositis",
                                     "JIA",  "Ank.Spond",
                                     "UKBB"))]
plotter(x=tmp, limits.do=c("p3","pc")) +
  scale_y_continuous(breaks=0) +
  theme(legend.position="bottom",legend.justification="right",
        axis.title.x=element_blank())
ggsave("figures/fig5-pc13-mr.pdf",height=6,width=8,scale=1)
ggsave("figures/slides-pc13-mr.pdf",height=6,width=12,scale=1,pointsize=16)


################################################################################

