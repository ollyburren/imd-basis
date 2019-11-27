library(data.table)
library(magrittr)
library(ggplot2)

## I think that you computed latest summary stats for JIA (and only included basis SNPs). I use this file ~/share/Data/GWAS/jia-mar-2019/summary-stats-mar2019.csv

source("~/A/R/table1-functions.R")
################################################################################

## sample sizes
N <- list(egpa=list(N1=534,#c(97,147,5,48,68+119,29,21),
                    N0=6688),#c(5465,273,130,118,0+226,93,343)),
          egpa.neg=list(N1=352,N0=6688),
          ra=list(N1=14361,N0=42923),
          jia.overall = list(N0 = 9196, N1 = 3490L), 
          jia.rfneg=list(N0=9196,N1=204),
          jia.systemic = list(N0 = 9196, N1 = 219L), 
          jia.persistent_oligoarthritis = list(N0 = 9196, N1 = 907L), 
          jia.extended_oligoarthritis = list(N0 = 9196, N1 = 537L), 
          "jia.RF_neg polyarthritis" = list(N0 = 9196, N1 = 887L), 
          "jia.RF_pos oligoarthritis" = list(N0 = 9196, N1 = 204L), 
          jia.ERA = list(N0 = 9196, N1 = 267L),
          jia.PSA = list(N0 = 9196, N1 = 233L),
          mg.late=list(N0=1977,N1=737),
          mg=list(N0=1977,N1=972),
          psc=list(N0=12019,N1=2871),
          pbc=list(N0=10475,N1=2764),
          vit=list(N0=37405,N1=2853),
          t1d=list(N0=8829,N1=5913),
          cro=list(N0=28072,N1=12194),
          uc=list(N0=33609,N1=12366),
          ## asthma=list(N0=107715,N1=19954))
asthma=list(N0=c(2110, 266, 94,257,96,373,213,797,364,123,216,555,267,453,179,1675,70,113,201,49,317,328,986,843,643,482,194,607,279,385,239,188,210,841,112,179,177,630,172,237,371,197,62,269,61,256,112,136,25,290,451,143,67,79,463,452),
            N1=c(3857,156,250,267,402,390,223,6463,3502,910,2005,6923,2381,9713,3058,33408,240,850,1749,999,545,870,5505,580,959,598,654,609,620,926,246,390,419,851,116,161,399,572,187,356,557,91,331,209,139,427,107,127,37,974,2416,4662,1616,1727,469,517)))







targets <- list(
## EGPA combined rs12405671 1 117263868 2.99E-06 3.70E-04 CD2, CD28 RA 1e-07
## EGPA ANCA- rs12405671 1 117263868 4.06E-05 3.04E-03 CD2, CD28 RA 1e-07
rs12405671=list(pos=117263868,chr=1,test=c("egpa","egpa.neg"),cond="ra"),

## EGPA combined rs1457115 5 110567598 3.21E-05 1.98E-03 TSLP, WDR36, CAMK4 asthma Unlinked to rs1837253 REF (r2=0.01)
## EGPA ANCA- rs1457115 5 110567598 2.16E-04 8.01E-03 TSLP, WDR36, CAMK4 asthma Unlinked to rs1837253 REF (r2=0.01)
rs1457115=list(pos=110567598,chr=5,test=c("egpa","egpa.neg"),cond="asthma"),

## EGPA combined rs10876864 12 56401085 1.19E-04 4.42E-03 SUOX, IKZF4 T1D VIT
rs10876864=list(chr=12,pos=56401085,test="egpa",cond=c("t1d","vit")),

## JIA rf-
rs9594746=list(chr= 13,pos=42989660, test=c("jia.RF_neg polyarthritis","jia.overall"),
               cond="pbc"),

## Myasthenia gravis combined
rs2188962=list(chr= 5,pos=131770805,test=c("mg","mg.late"),cond=c("asthma","cro"))

                )

## makes datasets smaller
## getchr <- sapply(targets,"[[","chr")  %>% unlist()  %>% unique()
## minpos <- sapply(targets,"[[","pos")  %>% unlist()  %>% as.numeric() %>% min()
## maxpos <- sapply(targets,"[[","pos")  %>% unlist()  %>% as.numeric() %>% max()
################################################################################

## 1kg alleles
if(!file.exists("~/alleles.RData")) {
    alleles <- lapply(targets,geta)
    alleles %<>% rbindlist()
    save(alleles, file="~/alleles.RData")
} else {
    load("~/alleles.RData")
}


################################################################################

## load disease data
if(FALSE) {

    data <- list()
fsave <- function() {
    sapply(data, nrow)  %>% print()
    ## save(data,  file="~/scratch/data2.RData")
}

edir <- "~/share/as_basis/egpa_lmm_with_pids"
files <- list.files(edir,full=TRUE)
files

    
    
message("--- egpa.neg ---")
neg <- readRDS(files[2])
neg[,c("chr","position"):=tstrsplit(pid,":")]
neg <- neg[ limit_chr_pos(chr, position) ]
neg[,sw:=switch_indicator(pid,a1,a2)]
neg[,position:=as.numeric(position)]
neg[,c("beta","vbeta"):=orp2bv(or,p.value)]
head(neg)
data$egpa.neg <- neg[,.(pid,chr,position,sw,beta,vbeta,p.value,trait="egpa.neg")]
fsave()

message("--- egpa ---")
comb <- readRDS(files[1])
comb[,c("chr","position"):=tstrsplit(pid,":")]
comb <- comb[ limit_chr_pos(chr, position) ]
comb[,sw:=switch_indicator(pid,a1,a2)]
comb[,position:=as.numeric(position)]
comb[,c("beta","vbeta"):=orp2bv(or,p.value)]
data$egpa <- comb[,.(pid,chr,position,sw,beta,vbeta,p.value,trait="egpa")]
fsave()

## MG imputed summary stats ~/share/Data/GWAS-summary/MyastheniaGravis_Renton_JAMA_Neurol_2015/ssimp_imputed
f.z2b <- function(z, af, n) {
    # z = imputed z statistics
    # af = allele frequency
    # n = sample size (effective)
    se.b <- 1/sqrt(2* af * (1-af) * n)
    b <- z * se.b
    return(list(b, se.b^2))
}

message("--- mg ---")
mg <- fread("~/share/Data/GWAS-summary/MyastheniaGravis_Renton_JAMA_Neurol_2015/ssimp_imputed/MyastheniaGravis_Overall_Renton_JAMA_Neurol_2015_ssimp_imputed.tab")[r2.pred>0.3]
mg <- mg[limit_chr_pos(chr,pos)]
setnames(mg,c("pos","P.imp"),c("position","p.value"))
mg[,pid:=paste(chr,position,sep=":")]
head(mg)
mg[,sw:=switch_indicator(pid,Allele1,Allele2)]
mg[,c("beta","vbeta"):=f.z2b(z_imp, maf, N.imp)]
data$mg <- mg[,.(pid,chr,position,sw,
                 beta,vbeta,p.value,N.imp,trait="mg")]
fsave()

message("--- mg.late ---")
mg.late <- fread("~/share/Data/GWAS-summary/MyastheniaGravis_Renton_JAMA_Neurol_2015/ssimp_imputed/MyastheniaGravis_LateOnset_Renton_JAMA_Neurol_2015_ssimp_imputed.tab")[r2.pred>0.3]
mg.late <- mg.late[limit_chr_pos(chr,pos)]
setnames(mg.late,c("pos","P.imp"),c("position","p.value"))
mg.late[,pid:=paste(chr,position,sep=":")]
head(mg.late)
mg.late[,sw:=switch_indicator(pid,Allele1,Allele2)]
mg.late[,c("beta","vbeta"):=f.z2b(z_imp, maf, N.imp)]
data$mg.late <- mg.late[,.(pid,chr,position,sw,
                           beta,vbeta,p.value,N.imp,trait="mg.late")]
fsave()

## JIA - TODO - all SNPs
## modify ~cew54/Projects/auto-basis/gwas-jia.R
message("--- jia.rfneg ---")
jia <- fread("/home/cew54/share/Data/GWAS/jia-mar-2019/summary-stats-chr13-rfneg-nov19.csv")
head(jia)
setnames(jia,c("b","v"),c("beta","vbeta"))
jia[,pid:=paste(chromosome,position,sep=":")]
jia[,sw:=switch_indicator(pid,allele.2,allele.1)]
data$jia.rfneg <- jia[,.(pid=paste(chromosome,position,sep=":"),
                         chr=chromosome,position,beta,vbeta,
                         sw,
                         p.value=pchisq(beta^2/vbeta,df=1,lower=FALSE),
                         trait="jia.rfneg")]

files <- list.files("/home/cew54/share/Data/GWAS/jia-mar-2019",
                    pattern="summary-stats-chr13-[0-7]-nov19.csv",
                    full=TRUE)
files

ilar <- read.table(text="0=overall=3490
1=systemic=219
2=persistent oligoarthritis=907
3=extended oligoarthritis=537
4=RF neg polyarthritis=887
5=RF pos oligoarthritis=204
6=ERA=267
7=PSA=233",sep="=")
colnames(ilar) <- c("n","w","N1")
ilar$w  %<>% sub(" ","_",.)

for(f in files) {
    jia <- fread(f)
    head(jia)
    setnames(jia,c("b","v"),c("beta","vbeta"))
    jia[,pid:=paste(chromosome,position,sep=":")]
    jia[,sw:=switch_indicator(pid,allele.2,allele.1)]
    i <- strsplit(f,"-")[[1]][6]
    m <-  match(i,ilar$n)
    nm <- paste0("jia.",ilar$w[m])
    N[[nm]] <- list(N0=9196,N1=ilar$N[m])
    data[[nm]] <- jia[,.(pid=paste(chromosome,position,sep=":"),
                         chr=chromosome,position,beta,vbeta,
                         sw,
                         p.value=pchisq(beta^2/vbeta,df=1,lower=FALSE),
                         trait=nm)]
}
## targets[[4]]$test <- paste0("jia.",ilar$w)

fsave()

## Basis traits
## CD GWAS-summary/asthma/cd_build37_40266_20161107.txt.gz
message("--- cro ---")
cro <- fread("zcat ~/share/Data/GWAS-summary/cd_build37_40266_20161107.txt.gz")
cro[,chr:=sub(":.*","",MarkerName)]
cro[,pos:=as.numeric(gsub(".*:|_.*","",MarkerName))]
cro <- cro[limit_chr_pos(chr,pos)]
cro[,pid:=paste(chr,pos,sep=":")]
head(cro)
cro[,sw:=switch_indicator(pid,toupper(Allele1),toupper(Allele2))]
data$cro <- cro[,.(pid=paste(chr,pos,sep=":"),chr,position=pos,sw,
                   beta=Effect,vbeta=StdErr^2,p.value=P.value,trait="cro")]
fsave()

## RA GWAS-summary/RA_GWASmeta_European_v2.txt.gz
message("--- ra ---")
ra <- fread("cat ~/share/Data/GWAS-summary/RA_GWASmeta_European_v2.txt.gz")
setnames(ra,make.names(names(ra)))
ra <- ra[limit_chr_pos(Chr,Position.hg19.)]
setnames(ra,c("Chr","Position.hg19.","OR.A1.","P.val"),
         c("chr","position","or","p"),skip_absent=TRUE)
ra[,pid:=paste(chr,position,sep=":")]
head(ra)
ra[,sw:=switch_indicator(pid,A1,A2)]
ra[,beta:=log(or)]
ra[,se:=(log(OR_95.CIup) - beta)/1.96]
data$ra <- ra[,.(pid=paste(chr,position,sep=":"),chr,position,sw,
                 beta,vbeta=se^2,p.value=p,trait="ra")]
fsave()

## Myogen imputed summary stats ~/share/Data/GWAS-summary/MYOGEN/imputed_summary_stats
## myogen <-
  
## Asthma GWAS-summary/asthma/TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv.gz
message("--- asthma ---")
asthma <- fread("zcat  ~/share/Data/GWAS-summary/asthma/TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv.gz")
setnames(asthma,make.names(names(asthma)))
asthma <- asthma[limit_chr_pos(chr,position)]
setnames(asthma,c("European_ancestry_beta_fix","European_ancestry_se_fix"),
         c("beta","se"))
asthma[,pid:=paste(chr,position,sep=":")]
head(asthma)
asthma[,sw:=switch_indicator(pid,reference_allele,alternate_allele)]
data$asthma <- asthma[,.(pid=paste(chr,position,sep=":"),chr,position,sw,
                         beta,vbeta=se^2,p.value=European_ancestry_pval_fix,
                         trait="asthma")]
fsave()

## VIT GWAS-summary/vitiligo-jin/raw/
message("--- vit ---")
vit <- fread("zcat ~/share/Data/GWAS-summary/vitiligo-jin/raw/GWAS123chr12cmh.txt.gz")
head(vit)
vit <- vit[limit_chr_pos(CHR,BP)]
vit[,pid:=paste(CHR,BP,sep=":")]
vit[,sw:=switch_indicator(pid,A1,A2)]
data$vit <- vit[,.(pid=paste(CHR,BP,sep=":"), chr=CHR,position=BP,sw,
                   beta=log(ORX),vbeta=SE^2, p.value=P,trait="vit")]
fsave()

## T1D GWAS-summary/t1d_cooper_2017.txt
t1d <- fread("zcat ~/share/Data/GWAS-summary/t1d_cooper_2017.txt.gz")
t1d <- t1d[limit_chr_pos(chromosome,position)]
t1d[,pid:=paste(chromosome,position,sep=":")]
head(t1d)
t1d[,sw:=switch_indicator(pid,a0,a1,TRUE)]
data$t1d <- t1d[!is.na(beta.meta),.(pid,chr=chromosome,position,sw,
                   beta=beta.meta,
                   vbeta=se.meta^2, p.value=p.meta^2,
                   trait="t1d")]
head(t1d)
fsave()


## PBC GWAS-summary/PBC-gwas-cordell-2015/GCMETA_fixedeffects.gz
pbc <- fread("zcat ~/share/Data/GWAS-summary/hg19_gwas_pbc_cordell_4_20_0.tab.gz")
pbc <- pbc[limit_chr_pos(Chr,Position)]
setnames(pbc,make.names(names(pbc)))
pbc[,beta:=log(OR.MinAllele.)]
pbc[,se:=(log(UpperOR) - beta)/1.96]
pbc[,pid:=paste(Chr,Position,sep=":")]
head(pbc)
pbc[,c("a1","a2"):=tstrsplit(Alleles.Maj.Min.,">")]
pbc[,sw:=switch_indicator(pid,a1,a2)]
data$pbc <- pbc[,.(pid=paste(Chr,Position,sep=":"),
                   chr=Chr,position=Position,sw,
                   beta=beta,vbeta=se^2,
                   p.value=PValue,trait="pbc")]
fsave()


## check
for(i in 1:(length(data)-1)) {
    for(j in (i+1):length(data)) {
        nint <- length(intersect(data[[i]]$pid, data[[j]]$pid))
        if(nint < 100) 
        cat(names(data)[i], names(data)[j], nrow(data[[i]]), nrow(data[[j]]),
            nint,"\n",sep="\t")
    }
}


################################################################################

## add data to targets
for(nm in names(targets)) {
    target <- targets[[nm]]
    diseases <- unlist(target[c("test","cond")])
    x=lapply(data[diseases], function(x) x[chr==target$chr &
                                           position > target$pos - 1e+6 &
                                           position < target$pos + 1e+6])
    names(x) <- diseases
    targets[[nm]]$data <- x
}

targets <- lapply(targets, fcoloc)
save(targets,file="~/scratch/targets.RData")

}
################################################################################

################################################################################

## analysis

(load("~/scratch/targets.RData"))
    
(load("~/alleles.RData"))

library(data.table)
library(magrittr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(egg)
library( 'gridExtra' )
library(grid)
library(gtable)

col1 <- "darkslateblue"

## save(targets,file="~/targets.RData")
## make manhattans
tt2 <- ttheme_minimal(core=list(fg_params=list(hjust=1, x=0.9),
                                bg_params = list(fill = blues9[2], col=NA)),
                      rowhead=list(fg_params=list(fontface=1,hjust=1, x=0.95)),
                      colhead=list(fg_params=list(fontface=1,hjust=1, x=0.95)))


mann <- lapply_with_names(names(targets),plotmann,pp4=TRUE)

## targeta <- targets[[4]]
## targeta$test <- "jia.RF_neg polyarthritis"
## targeta$data <- targeta$data[c("jia.RF_neg polyarthritis","pbc")]
## targeta$coloc <- targeta$coloc["jia.RF_neg polyarthritis"]
## targets$sub <- targeta
## plotmann("sub",ti=names(targets)[4])
## plotmann(names(targets)[4],addsnps=FALSE,pp4=FALSE) 

mann[[1]]
mann[[2]] # 2 signals?
mann[[3]]
mann[[4]] # JIA
mann[[5]] # 2 signals

if(!file.exists("~/LD.RData")) {
    LD <- vector("list",length(targets))
    names(LD) <- names(targets)
    LD[[1]]  <- getld(targets[[1]])
    LD[[2]]  <- getld(targets[[2]])
    LD[[3]]  <- getld(targets[[3]])
    LD[[4]]  <- getld(targets[[4]])
    LD[[5]]  <- getld(targets[[5]])
    save(LD, file="~/LD.RData")
} else {
    (load("~/LD.RData"))
    LD <- LD[names(targets)]
}

splits <- lapply_copy_names(LD, function(x) splitter(x$ld^2,thr=0.2))
sapply(splits,length)
save(splits, file="~/splits.RData")

(load("~/splits.RData"))

for(nm in names(targets)[1:5]) {
    targets[[nm]]$scoloc <- fcoloc_ld(nm)
}

getc <- function(what) {
    tmp <- lapply(targets,"[[",what)
    lapply(tmp,function(x) lapply(x, as.data.frame)  %>% do.call("rbind",.))
}

getc("coloc")

## $rs12405671
##             ra.PP4 ra.best.snp ra.best.pp
## egpa     0.7791986 1:117260856  0.6431025
## egpa.neg 0.6132652 1:117263790  0.7222521

## $rs1457115
##          asthma.PP4 asthma.best.snp asthma.best.pp
## egpa      0.9666964     5:110558064      0.9963530
## egpa.neg  0.2721288     5:110558064      0.9603106

## $rs10876864
##        t1d.PP4 t1d.best.snp t1d.best.pp   vit.PP4 vit.best.snp vit.best.pp
## egpa 0.4290869  12:56473808    0.287724 0.5648692  12:56396768   0.2122132

## $rs9594746
##                            pbc.PP4 pbc.best.snp pbc.best.pp
## jia.RF_neg polyarthritis 0.9682804  13:42989660   0.9096884
## jia.overall              0.9367632  13:42989660   0.9437064

## $rs2188962
##         asthma.PP4 asthma.best.snp asthma.best.pp   cro.PP4 cro.best.snp
## mg      0.06301913     5:131995964      0.6528129 0.7529088  5:131770805
## mg.late 0.03341388     5:131995964      0.5502870 0.7963004  5:131770805
##         cro.best.pp
## mg        0.6378797
## mg.late   0.7008283

getc("scoloc")
## $rs12405671
##             ra.PP4 ra.best.snp ra.best.pp
## egpa     0.6345481 1:117263790  0.4517748
## egpa.neg 0.6594816 1:117263790  0.8996040

## $rs1457115
##          asthma.PP4 asthma.best.snp asthma.best.pp
## egpa     0.08808766     5:110401872      0.9232543
## egpa.neg 0.09247388     5:110401872      0.8558965

## $rs10876864
##        t1d.PP4 t1d.best.snp t1d.best.pp   vit.PP4 vit.best.snp vit.best.pp
## egpa 0.4987594  12:56473808    0.303807 0.6504957  12:56396768    0.239566

## $rs9594746
##                            pbc.PP4 pbc.best.snp pbc.best.pp
## jia.RF_neg polyarthritis 0.9701729  13:42989660   0.9097862
## jia.overall              0.9529662  13:42989660   0.9487063

## $rs2188962
##         asthma.PP4 asthma.best.snp asthma.best.pp   cro.PP4 cro.best.snp
## mg       0.7863205     5:131797547      0.7116299 0.7580343  5:131770805
## mg.late  0.7756247     5:131797547      0.6138965 0.8043486  5:131770805
##         cro.best.pp
## mg        0.6462316
## mg.late   0.7089781


mann_ld <- lapply_with_names(names(targets)[1:5],plotmann_ld,pp4=TRUE)

pdf("~/mann1.pdf",height=8,width=8)
plot_grid(mann[[1]],nrow=1)
plot_grid(mann[[2]],nrow=1)
plot_grid(mann[[3]],nrow=1)
## plot_grid(mann[[4]],nrow=1)
plot_grid(mann[[5]],nrow=1)
dev.off()


pdf("~/mann2.pdf",height=8,width=8)
plot_grid(mann_ld[[1]],nrow=1)
plot_grid(mann_ld[[2]],nrow=1)
plot_grid(mann_ld[[3]],nrow=1)
## plot_grid(mann_ld[[4]],nrow=1)
plot_grid(mann_ld[[5]],nrow=1)
dev.off()


pdf("~/mann.pdf",height=8,width=8)
plot_grid(mann[[1]], mann_ld[[1]],nrow=1)
plot_grid(mann[[2]], mann_ld[[2]],nrow=1)
plot_grid(mann[[3]], mann_ld[[3]],nrow=1)
plot_grid(mann[[4]], mann_ld[[4]],nrow=1)
plot_grid(mann[[5]], mann_ld[[5]],nrow=1)
dev.off()

save(targets, file="~/targets+scoloc.RData")

################################################################################

## done

if(!interactive())
    q("no")

################################################################################

## sandbox
lapply(targets,"[[","coloc")

## ra
13:42989660


    if(target.chr==5)
    man <- man + geom_vline(xintercept=110401872,linetype="dashed",col="red")
man




(load("~/genes.RData"))
genes %<>% as.data.table()
gg <- genes[chromosome_name==target.chr &
            ( abs(start_position-target.pos) < 1e+5 | abs(end_position-target.pos) < 1e+5) ]
gg <- gg[order(start_position)]
gg[,gene:=factor(external_gene_name,levels=gg$external_gene_name)]


g <- ggplot(gg, aes(x=start_position,xmin=start_position,xmax=end_position,y=gene)) +
  geom_errorbarh()

library(cowplot)
plot_grid(man,g,nrow=2,rel_heights=c(3,1))

## my.res <- coloc.abf(dataset1=list(beta=b1$beta, varbeta=b1$varbeta, N=nrow(X1),sdY=sd(Y1),type="quant"),
##                     dataset2=list(beta=b2$beta, varbeta=b2$varbeta, N=nrow(X2),sdY=sd(Y2),type="quant"),
##                     MAF=maf)
library(coloc)
egpa.ra <- coloc.abf(list(beta=x2$beta,
                          varbeta=x2$vbeta,
                          snp=x2$pid,
               type="cc",
               N=6688+534,
               s=534/(6688+534)),
          list(beta=y1$beta,
               varbeta=y1$vbeta,
               snp=y1$pid,
               type="cc",
               N=(42923+14361),
               s=14361/(42923+14361)))
          
neg.ra <- coloc.abf(list(beta=x1$beta,
                         varbeta=x1$vbeta,
                         snp=x1$pid,
               type="cc",
               N=6688+352,
               s=352/(6688+352)),
          list(beta=y1$beta,
               varbeta=y1$vbeta,
               snp=y1$pid,
               type="cc",
               N=(42923+14361),
               s=14361/(42923+14361)))
          

               
