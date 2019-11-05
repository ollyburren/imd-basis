library(devtools)
install_github('ollyburren/cupcake@v0.1.0.0')
library(cupcake)
ROOT.DIR <- '~/share/as_basis/basis-creation'
SNP_MANIFEST_FILE <- file.path(ROOT.DIR,'support/13_snp_manifest.RDS')
TRAIT_MANIFEST_FILE <- file.path(ROOT.DIR,'support/13_trait_manifest.tab')
SHRINKAGE_FILE <- file.path(ROOT.DIR,'support/13_trait_shrinkage.RDS')
GWAS_DATA_DIR <- file.path(ROOT.DIR,'/sum_stats/')
P.thresh <- 0.001
seb.thresh <- 0.2

gwas.DT<-get_gwas_data(TRAIT_MANIFEST_FILE,SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE)
shrink.DT<-compute_shrinkage_metrics(gwas.DT)

dt2mat <- function(DT,var){
  B <- dcast(DT,pid ~ trait,value.var=var)
  tmp.mat <- as.matrix(B[,-1]) %>% t()
  colnames(tmp.mat) <- B[,1]$pid
  tmp.mat
}

## comparison with degas approach
## for our weighted approach we need to remove the control vector prior to projection
## so that we have a valid comparison

M <- merge(gwas.DT,shrink.DT,by='pid')
M[,metric:=lor * shrinkage]
cupcake.beta <- dt2mat(M, var="lor") %>% prcomp(.,center=TRUE,scale=FALSE)
cupcake.shrinkage <- dt2mat(M, var="metric") %>% prcomp(.,center=TRUE,scale=FALSE)


gwas.DT[,lor:=log(or)][,z:=qnorm(p.value/2,lower=FALSE) * sign(lor)][,se:=lor/z]
filt.DT <- copy(gwas.DT)
## get pruned SNPs see ob-ldprune.R
prune.snps <- readRDS(SNP_LD_PRUNE_MANIFEST_FILE)[]$pid %>% unique
filt.DT <- gwas.DT[pid %in% prune.snps,]
filt.DT[p.value > P.thresh | se > seb.thresh, c("lor","z"):=list(0,0)]
## only get snps for which there is a non-zero value for at least one SNP
keep.snps <- filt.DT[filt.DT[,.I[max(abs(lor)!=0)],by='pid']$V1]$pid %>% unique
degas.beta <- dt2mat(filt.DT[pid %in% keep.snps,], var="lor") %>% prcomp(.,center=TRUE,scale=FALSE)
degas.z <- dt2mat(filt.DT[pid %in% keep.snps,], var="z") %>% prcomp(.,center=TRUE,scale=FALSE)
save(list=c('cupcake.beta','cupcake.shrinkage','degas.beta','degas.z'),file=DEGAS_COMPARISON_BASIS_FILE)
