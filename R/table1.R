library(data.table)
library(magrittr)

## load results
source("R/cw-reader.R")

## supp tab 10
x <- reader()
fwrite(x[,.(category.label,trait.label,PC,delta,var.delta=var.proj,fdr.delta=newfdr,p.overall,fdr.overall)],
file="figures/supptable-10-projections.csv.gz")

## source("R/cw-renamer.R")
## res <- reader(J=2)
## res$proj[,p.adj:=p.adjust(p.value,method="BH"),by=c("fdrcat","variable")]
## proj11 <- res$proj[(category %in% c("myogen","lyons_egpa",
##                                                "renton_mg","estrada_NMO",
##                                                "bowes_jia_2019")) &
##                  !(trait %in% c("myositis_myogen","jia_case_19"))]
## subset to small datasets (<2000 cases)
proj <- reader()[category %in% c("Ank.Spond", "birdshot_retinopathy", "EGPA", "JIA", "methotrexate", "Myasthenia gravis", "NMO", "PsA", "Uveitis", "Vasculitis") &
                 !(trait %in% c("bowes_psa","ank_spond","jia_case_19")) &
                 fdr.overall < 0.01 &
                 newfdr < 0.01]

with(proj,table(trait.label, category.label))

# list for Methods section
with(proj,cat(unique(paste(category.label, trait.label, sep=" ")),sep="\n"))

## number of tests in parent family
## N.parents=nrow(proj13)
## proj11 <- proj11[p.adj<0.01]
traits.considered=reader()[category %in% c("Ank.Spond", "birdshot_retinopathy", "EGPA", "JIA", "methotrexate", "Myasthenia gravis", "NMO", "PsA", "Uveitis", "Vasculitis") &
                             !(trait %in% c("bowes_psa","ank_spond","jia_case_19"))]
N.traits.considered <- length(unique(traits.considered$trait))
N.traits.signif <- length(unique(proj$trait))

unique(traits.considered[,.(category.label, trait.label, fdr.overall)])

sort(table(proj$trait))

dim(proj) ## 22 trait/pc pairs
length(unique(proj$trait)) # 12 unique traits

## load input data
(load(SPARSE_BASIS_FILE))
read_with_fdr <- function(tr) {
    d <- read_raw(tr)
    d[,fdr.gw:=p.adjust(p.value,method="BH")]
    d[pid %in% rownames(rot.pca)]
}
data <- lapply(unique(proj$trait), read_with_fdr)  %>%  rbindlist(.,use.names=TRUE)

## basis drivers
source("R/cw-annot-snps.R")
## source("~/Projects/auto-basis/scripts/plots.R")
library(ggplot2)
library(cowplot) ; theme_set(theme_cowplot())

pids.use <- rownames(use.pca)
length(pids.use) # 566 snps in sparse basis

m <- match(data$pid,rownames(use.pca))
for(j in colnames(use.pca)) {
    data[[j]] <- ifelse(is.na(m), FALSE, use.pca[m,j])
}

## look at FDR
FDR <- lapply(1:nrow(proj), function(r) {
    pc <- as.character(proj$PC[[r]])
    dsub <- data[ trait==proj$trait[[r]] & data[[pc]]==TRUE & !is.na(p.value), ]
    if(!nrow(dsub))
        return(NULL)
    dsub[,fdr.pc:=p.adjust(p.value,method="BH")] # FDR within driver snps
    dsub[,p.bonf:=p.value*.N,by=pc]
    dsub[,.(row=r,trait,pc=pc,pid,p.value,fdr.gw,fdr.pc,p.bonf)]
})  %>% rbindlist()

## annotate
length(unique(FDR$pid))
fdr.ann <- annot(unique(FDR$pid),snp=FALSE)
fdr.ann <- as.data.table(fdr.ann)

## reduce genes
fdr.red <- fdr.ann[,.(gene=paste(external_gene_name,collapse="/"),
                      egene=paste(ensembl_gene_id,collapse="/")),
                   by=c("pid")]

## add rsid
dbsnp <- fread(cmd="zcat /home/cew54/share/Data/reference/dbsnp-common_all.vcf.gz | grep -v '##' | cut -f1-5")
head(dbsnp)
pids  <- tstrsplit(fdr.red$pid,":")  %>% as.data.table()
setnames(pids,c("chr","pos"))
pids[,pos:=as.integer(pos)]
setnames(dbsnp,c("#CHROM","POS"),c("chr","pos"))
dim(pids)
pids <- merge(pids,dbsnp,all.x=TRUE,by=c("chr","pos"))
which(duplicated(pids$ID))
head(pids)


## library(biomaRt)
## ## variation = useEnsembl(biomart="snp")
## ## listDatasets(variation)
## variation = useEnsembl(biomart="snp", dataset="hsapiens_snp")
## listFilters(variation)
## listAttributes(variation)

## rs1333049 <- getBM(attributes=c('refsnp_id','refsnp_source','chr_name','chrom_start','chrom_end','minor_allele','minor_allele_freq','minor_allele_count','consequence_allele_string','ensembl_gene_stable_id','ensembl_transcript_stable_id'),
## filters = 'snp_filter', values ="rs1333049",
## mart = variation)
## rs1333049

## ss <- tstrsplit(fdr.red$pid[1:3],":")
## rsid <- getBM(attributes=c('refsnp_id','refsnp_source','chr_name','chrom_start','consequence_allele_string'),
## filters = c("chr_name","start"),
## values =ss,
## mart = variation)



FDR <- merge(FDR,fdr.red,by="pid")
FDR <- merge(FDR,
             unique(proj[,.(trait,category.label,trait.label)]),
             by="trait")
uFDR <- unique(FDR[,.(trait,pid,p.value)])
uFDR[,fdr.drivers:=p.adjust(p.value,method="BH"),by="trait"]
FDR <- merge(FDR,uFDR[,.(trait,pid,fdr.drivers)],by=c("trait","pid"))
FDR[,lab:=paste(category.label,trait.label,pc,sep=" / ")]

my_stat_qq = function(data, ...) {
}

FDR[,sigcol:=fdr.drivers<0.01]
my_stat_qq(FDR, "sigcol")

data <- copy(FDR)
data <- data[order(data$p.value,decreasing=TRUE),][,sample:=-log(p.value)]
data[,theoretical:=stats::qexp((1:.N)/(.N+1)),by=lab]

## ggplot(data) + 
##   geom_point(aes_string(x="theoretical", y="sample", colour="sigcol")) +
##   facet_wrap(~lab) +
##   geom_abline()

## pqq <- ggplot(FDR,aes(sample=-log(p.value),group=lab)) +
##   stat_qq(aes(col=fdr.drivers<0.01), distribution = stats::qexp) +
##   stat_qq_line(distribution = stats::qexp) +
##   facet_wrap(.~lab,scales="free_y") +
##   labs(x="Expected",y="Observed")
## pqq

pqq <- ggplot(data,aes(x=theoretical,y=sample,group=lab)) +
  geom_point(aes(col=sigcol)) +
  geom_abline() +
  facet_wrap(.~lab,scales="free_y",nrow=6) +
  scale_colour_manual("FDR<0.01",values=c("FALSE"="grey","TRUE"="black")) +
  labs(x="Expected",y="Observed") +
  theme(legend.position="bottom")
pqq

ggsave("figures/suppfig-sparsesig-qqplots.pdf",height=8,width=8,scale=1.5)
## pdf("figures/suppfig-sparsesig-qqplots.pdf",height=12,width=12)
## print(pqq)
## dev.off()


## phist <- ggplot(FDR,aes(x=p.value,group=lab)) +
##   ## geom_density() +
##   ## facet_wrap(~lab) 
##   geom_histogram(col="grey",fill="grey20",binwidth=0.05) +
##   facet_wrap(.~lab,scales="free_y") 
## phist

pfdr <- ggplot(FDR,aes(x=-log10(p.value),y=-log10(fdr.gw))) + geom_point()  +
  geom_point(aes(y=-log10(fdr.pc)),col="red") +
  geom_point(aes(y=-log10(fdr.drivers)),col="green") +
  geom_vline(xintercept=-log10(5e-8)) +
  geom_hline(yintercept=-log10(0.01)) +
  geom_smooth() +  background_grid() +
  xlim(0,13) + ylim(0,11)
## plot_grid(phist,pfdr,nrow=1,labels=c("a","b"))
pfdr
## FDR 0.01 approx GW sig across all SNPs

unique(FDR[p.value<5e-8,.(category.label,trait.label,pid)]) #,RefSNP_id)])

sig <- FDR[fdr.drivers<0.01,.(p.value=p.value[which.min(fdr.drivers)],
                              fdr=min(fdr.drivers),
                              fdr.gw,
                       components=paste(sort(unique(pc)),collapse=":")),
           by=c("pid","category.label","trait.label","trait"#,"RefSNP_id"
)]  %>%
    unique(.,by=c("pid","trait"))
## compare to calling per-PC
unique(FDR[fdr.drivers>0.01 & fdr.pc<0.01 & p.value > 5e-8,.(category.label,trait.label,pid,#RefSNP_id,
p.value)])
sig[,c("chr","pos"):=lapply(tstrsplit(pid,":"),"as.numeric")]


## when 2+ snps are significant near each other, drop less significant ones
ssig <- split(sig,sig$trait)
x=ssig[[1]]
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param x
##' @return
##' @author Chris Wallace
f <- function(x) {
    x <- x[order(p.value)]
    xld <- as.matrix(LD[x$pid,x$pid])
    diag(xld) <- 0
    xld[lower.tri(xld)] <- 0
    x[,maxld:=apply(abs(xld),2,max)] # this will be > 0 if there is another SNP with smaller fdr in LD
    x[abs(maxld)<0.2]
}
thinsig <- split(sig,sig$trait)  %>% lapply(., f)  %>% rbindlist()

dim(sig)
dim(thinsig)

## add basis p values
if(!exists("basis.DT")) {
    source("~/Projects/auto-basis/cupcake/R/cupfunc.R")
    TRAIT_MANIFEST <- '~/share/as_basis/GWAS/trait_manifest/as_manifest_gwas.tab'
    GWAS_DATA_DIR <- '~/share/as_basis/GWAS/sum_stats'
    SNP_MANIFEST_FILE <-'~/share/as_basis/GWAS/snp_manifest/gwas_june_19_w_vitiligo.tab'
    basis.DT<-get_gwas_data(TRAIT_MANIFEST,SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE)
    }
gwsig <- basis.DT[pid %in% sig$pid & p.value<5e-8,][,.(basis.diseases=paste(sort(trait),collapse=",")),by="pid"]
sig <- merge(thinsig,gwsig,by="pid",all.x=TRUE)
ugwsig <- basis.DT[pid %in% setdiff(sig$pid,gwsig$pid),][,.(basis.diseases2=paste(trait[which.min(p.value)],format.pval(min(p.value),digits=2))),by="pid"]
sig <- merge(sig,ugwsig,by="pid",all.x=TRUE)
sig[is.na(basis.diseases),basis.diseases:=basis.diseases2]
sig[,basis.diseases2:=NULL]

## get snp comments from
comments=fread("~/basis-drivers-comments.csv")
comments  %<>% rename.traits()
## comments[,RefSNP_id:=NULL]
setnames(comments,make.names(names(comments)))
comments[,basis.traits.with.GWsig.assoc:=NULL]
table(comments$trait.label)
table(sig$trait.label)
comments[,trait.label:=sub("EGPA","",trait.label)]

g <- unique(comments[,.(pid,gene)])
nog <- comments[,.(trait.label,pid,RefSNP_id,comments,gwsig,other.evidence,novel)]
sig <- merge(sig,nog,by=c("pid","trait.label"),all.x=TRUE)
sig <- merge(sig,g,by=c("pid"),all.x=TRUE)
## url <- "https://docs.google.com/spreadsheets/d/1W6tgXaYkyyKnjo3vVTGc7Gfex2w373T86h_Wa6qSAc0/edit?usp=sharing"
## library(googledrive)
## drive_download(
##     url,
##     path = "comments.csv",
##     overwrite = TRUE
## )
sig[is.na(RefSNP_id)]

dim(sig)
## update standard blank fields
known <- sig[!is.na(gene)]$pid
sig[pid %in% known, gene:=setdiff(gene,NA),by="pid"]
sig[is.na(comments),comments:=""]
sig[is.na(gene),gene:=""]
setnames(sig,"gene","genelist")

sig[is.na(gwsig)]
sig <- sig[order(gwsig,other.evidence,trait)]

#18:12847136

## final table - anything missing analysis?
sig[is.na(genelist),.(category.label,trait.label,pid,RefSNP_id,p.value,fdr,genelist,basis.diseases,gwsig,other.evidence,novel,comments)]  %>% tail
sig[,.(category.label,trait.label,pid,RefSNP_id,p.value,fdr,genelist,basis.diseases,gwsig,other.evidence,novel,comments)]  

sig[order(p.value),.(category.label,trait.label,pid,RefSNP_id,p.value,fdr,genelist,basis.diseases,gwsig,other.evidence,novel,comments)]  


sig[category.label=="JIA",.(category.label,trait.label,pid,p.value,fdr,basis.diseases)]

fwrite(sig[,.(category.label,trait.label,pid,RefSNP_id,p.value,fdr,genelist,basis.diseases,gwsig,other.evidence,novel,comments)],
       file="figures/table-drivers-significant.csv") 

sig[,chr:=sub(":.*","",pid)]
sig[order(chr,pid),.(pid,category.label,trait.label,RefSNP_id)]
sig[category.label=="JIA"][order(chr,pid),.(pid,category.label,trait.label,RefSNP_id,p.value,fdr)]

## stats
cat(c("number of traits considered: ",N.traits.considered,"\n",
      "number of traits significant: ",N.traits.signif,"\n",
      "number of trait/comp sig pairs: ",nrow(proj),"\n",
      "number of traits with sig comp: ",length(unique(proj$trait)),"\n",
      "number of traits with sig snp: ",length(unique(sig$trait)),"\n",
      "number of sig assocs: ",nrow(sig),"\n",
      "number of GW sig assocs: ",sum(sig$p.value< 5e-8),"\n"),
    sep="",file="stats-sparse-drivers-significance.txt")


if(!interactive())
    q("no")


sig[p.value<5e-8,.(fdr.gw)]

sig[,.(trait,pid,pc,p.value,gene)]

dim(sig) # 65 
with(sig,table(p.value<5e-8)) # 10 gw sig, 55 new
with(unique(sig,by=c("trait","pid")),table(p.value<5e-8)) # 4 gw sig, 44 new
with(unique(sig,by=c("pid")),table(p.value<5e-8)) # 2 gw sig, 26 new
length(unique(proj$trait)) # 12
length(unique(sig$trait)) # 8
summary(sig$p.value)

files=list.files(d,pattern="jia")[-1]
jia.traits <- sub("_19_source.RDS","",files)
jia <- lapply(jia.traits,function(tr) {
    tmp <- readRDS( file.path(d, paste0(tr,"_19_source.RDS")) )
    tmp[,trait:=tr]
})  %>% rbindlist(.,fill=TRUE)

alldata <- lapply(unique(sig$trait), function(tr) {
    message("loading ",tr)
    tr1  <-   sub("renton_mg_combined","renton_mg",tr)
    tmp <- readRDS( file.path(d, paste0(tr1,"_source.RDS")) )
    tmp[,fdr:=p.adjust(p.value,method="BH")]
    tmp[,trait:=tr]
})  %>% rbindlist(.,fill=TRUE)
head(alldata)

alldata[,c("chr","pos"):=lapply(tstrsplit(pid,":"),"as.numeric")]
alldata <- alldata[order(trait,chr,pos),]
alldata[,cpos:=cumsum(c(0,pmax(0,diff(pos)))),by="trait"]
alldata <- merge(alldata,sig[,.(trait,pid,sfdr=TRUE)],by=c("trait","pid"),all.x=TRUE)

ggplot(alldata[p.value<1e-2,],aes(x=cpos,y=-log10(p.value),col=factor(chr %% 2))) +
  geom_point(size=2,alpha=0.5) +
  geom_point(data=alldata[p.value<1e-3 & sfdr==TRUE & p.value>5e-8],size=4,col="black",pch=22) +
  facet_grid(trait~.,scales="free_y")


TRAIT_MANIFEST <- '~/share/as_basis/GWAS/trait_manifest/as_manifest_gwas.tab'
GWAS_DATA_DIR <- '~/share/as_basis/GWAS/sum_stats'
SNP_MANIFEST_FILE <-'~/share/as_basis/GWAS/snp_manifest/gwas_june_19_w_vitiligo.tab'
basis.DT<-get_gwas_data(TRAIT_MANIFEST,SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE) 
basis.DT[pid %in% sig$pid & p.value<5e-8,]

basis.DT[pid=="5:110567598" & p.value<5e-8]
basis.DT[pid=="5:131796922" & p.value<5e-8]
basis.DT[,c("chr","pos"):=lapply(tstrsplit(pid,":"),"as.numeric")]

## split by chr, trait
ss=split(alldata,alldata[,.(chr,trait)])
length(ss)

## function to calc min dist to sfdr=TRUE SNP
f=function(a) {
    w=which(a$sfdr==TRUE)
    if(length(w)) {
        tmp=lapply(w, function(i) { abs(a$pos-a$pos[i]) })  %>% do.call("cbind",.)
        tmpw=apply(tmp,1,which.min)
        tmpd=apply(tmp,1,min)
        a[,nearest:=pid[w[tmpw]]][,nearest.dist:=tmpd]
    } else {
        a[,nearest:=""][,nearest.dist:=Inf]
    }
    a
}
ss1=lapply(ss,f)    

ad1=rbindlist(ss1)
head(ad1)
ad1 <- ad1[nearest.dist<1e+5,]
ad1[,minp:=min(p.value),by=c("trait","nearest")]
ad1[,sfdr.proxy:=p.value==minp]
sum(ad1$sfdr.proxy,na.rm=TRUE)

ad1[pid=="1:117263868",]
ad1[nearest=="1:117263868",]

ad1[p.value==minp & nearest.dist>0,]

dim(ad1)



sum(ad1$sfdr,na.rm=TRUE)


alldata[,possig:=ifelse(sfdr==TRUE),pos,NA]
alldata[,mindistsig:=min(pos

int01 <- FDR[fdr0<0.01 & p.value>5e-8,] # 85 SNP/trait/PC trios
dim(int01)
int01[,.(pid,trait,pc,p.value,fdr,fdr0)]

unt01 <- unique(int01[,.(pid,trait,gene,p.value)])) # 50 unique SNP/trait pairs



FDR999[,pid:=paste0("P",pid)]
fwrite(FDR999[p.bonf < 0.05 | fdr0 < 0.01, ],file="basis-sig-drivers.csv",quote=TRUE)
