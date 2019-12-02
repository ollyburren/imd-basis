#module load r-3.6.1-gcc-5.4.0-zrytncq
#export R_LIBS_USER=/home/ob219/R/x86_64-pc-linux-gnu-library/3.6
#library(devtools)
#install_github("ollyburren/cupcake")

SUM_STATS_DIR <- '/home/ob219/share/as_basis/GWAS/for_fdr_13_traits_0919'
ROOT.DIR <- '~/share/as_basis/basis-creation'
SNP_MANIFEST_FILE <- file.path(ROOT.DIR,'support/13_snp_manifest.RDS')
TRAIT_MANIFEST_FILE <- file.path(ROOT.DIR,'support/13_trait_manifest.tab')
SHRINKAGE_FILE <- file.path(ROOT.DIR,'support/13_trait_shrinkage.RDS')
BASIS_FILE <- file.path(ROOT.DIR,'support/13_trait_basis.RDS')
GWAS_DATA_DIR <- file.path(ROOT.DIR,'/sum_stats/')



others <- 'SRC_FILE SSIMP_FILE  Label Category
pm_myogen gwas_pm_new_ssimp PM  Myositis
jdm_myogen  gwas_jdm_new_ssimp  JDM Myositis
dm_myogen gwas_dm_new_ssimp  DM Myositis
myositis_myogen  gwas_dmjdmpm_new_ssimp combined Myositis
renton_mg_early YoungOnset_ssimp  Early_onset Myasthenia_gravis
renton_mg_late  LateOnset_ssimp Late_onset  Myasthenia_gravis
renton_mg  Overall_ssimp  _combined Myasthenia_gravis
NMO_IgGPos  NMO_IgGPos_ssimp  IgGPos  NMO
NMO_IgGNeg  NMO_IgGNeg_ssimp  IgGNeg  NMO
NMO_combined NMO_combined_ssimp combined NMO
methotrexate  NA  Methotrexate  NA
hasnoot_uveitis_jia NA  Uveitis NA
birdshot_retinopathy  NA  Birdshot_chorioretinopathy  NA
jia_sys_19  NA  Systemic  JIA
jia_RFpos_19  NA  RFPos JIA
jia_RFneg_19  NA  RFNeg JIA
jia_PsA_19  NA  PsA JIA
jia_PO_19 NA  PO  JIA
jia_ERA_19  NA  ERA JIA
jia_EO_19 NA  EO  JIA
jia_case_19 NA  combined JIA
jia_undiff_19 NA Undifferentiated JIA
mpo_gwas1 NA  MPO+  Vasculitis
pr3_gwas1 NA  PR3+ Vasculitis
all.egpa_lmm  NA  _combined Vasculitis
anca.negative.egpa_lmm  NA ANCA-  EGPA
mpo.anca.positive.egpa_lmm  NA  MPO+  EGPA
span_psa  span_psa_ssimp  Spanish PSA
na_psa  na_psa_ssimp  N_American  PSA
bowes_psa NA  UK  PSA
li_ankspond NA  Turkish/Iranian Ank.Spond
ank_spond NA  International Ank.Spond
' %>% read.table(text=.,header=TRUE) %>% data.table

library(cupcake)




gwas.DT<-get_gwas_data(TRAIT_MANIFEST_FILE,SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE)
basis <- readRDS(BASIS_FILE)
pids <- rownames(basis$rot)
sDT <- readRDS(SHRINKAGE_FILE)


## compute perpendicular residuals from x=y
pres <- function(x,y){
  ## intercept of x=-y which is orthogonal to x=y
  c <- (x + y)/2
  ## pythagoras to compute perpendicular distance
  sqrt((x-c)^2 + (y-c)^2)
}

library(parallel)
all.res <- mclapply(1:nrow(others),function(i){
  trait <- sprintf("%s_%s",others$Category[i],others$Label[i]) %>% gsub("^NA\\_","",.)
  sprintf("Processing %s",trait) %>% message
  dat <- sprintf("%s_source.RDS",others$SRC_FILE[i]) %>% file.path(SUM_STATS_DIR,.) %>% readRDS
  dat <- dat[pid %in% pids,]
  miss.idx <- which(is.na(dat$or))
  miss.count <- length(miss.idx)
  if(miss.count==0)
    miss.idx<-1
  ## project on in the usual fashion.
  X<-cupcake::project_basis(dat[-miss.idx,],sDT,basis,trait)$proj
  ## next create a basis using only the SNPs present to see how this affects the basis
  filt.gwas.DT <- gwas.DT[pid %in% dat$pid[-miss.idx],]
  filt.sDT <-compute_shrinkage_metrics(filt.gwas.DT)
  filt.basis <- create_basis(filt.gwas.DT,filt.sDT)
  Y <- cupcake::project_basis(dat[-miss.idx,],filt.sDT,filt.basis,paste(trait,'_filt'))$proj
  ret <- rbind(X,Y) %>% reshape2::melt(.) %>% reshape2::dcast(Var2~Var1) %>% data.table
  setnames(ret,c('PC','ori','filt'))
  ret[,c('trait','missing'):=list(trait,miss.count)][]
},mc.cores=8) %>% rbindlist
others[,trait:=sprintf("%s_%s",Category,Label) %>% gsub("^NA\\_","",.)]
all.res[,tlabel:=NA]
#all.res[trait %in% others[!is.na(SSIMP_FILE),]$trait,tlabel:=trait]
## add labels to the ones that we impute
all.res[,residuals:=pres(ori,filt)]
all.res[,Z.res:=(residuals-mean(residuals))/sd(residuals),by=PC]
all.res[,label:='']
all.res[abs(Z.res)>1.96,tlabel:=trait]

library(cowplot)
library(ggplot2)
library(ggrepel)
theme_set(theme_cowplot())
pg <- ggplot(all.res[PC!='PC14',],aes(x=ori,y=filt,label=tlabel)) + geom_point() + facet_wrap(~PC,scales='free') +
geom_abline(col='red',lty=2) + geom_label_repel(box.padding=1.5,cex=3.5,label.size=0.01) + xlab("Original Projection Score") + ylab("Tailored Projection Score")
save_plot(file="~/tmp/missing.pdf",pg,base_height=9)
