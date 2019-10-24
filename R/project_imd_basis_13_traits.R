## project summary stats for all traits onto basis created in create_imd_basis_13_traits.R
## note this uses cupcake v0.1.0.1

library(devtools)
#load_all("~/git/cupcake")
library(parallel)
install_github('ollyburren/cupcake@v0.1.0.1')
ROOT.DIR <- '~/share/as_basis/basis-creation'
SNP_MANIFEST_FILE <- file.path(ROOT.DIR,'support/13_snp_manifest.RDS')
SHRINKAGE_FILE <- file.path(ROOT.DIR,'support/13_trait_shrinkage.RDS')
BASIS_FILE <- file.path(ROOT.DIR,'support/13_trait_basis.RDS')
NOWEIGHT_BASIS_FILE <- file.path(ROOT.DIR,'support/13_trait_basis-noweight.RDS')
DATA_DIR <- file.path(ROOT.DIR,'../basis-projection/project_13')
FULL_BASIS_PROJECTION_FILE <- file.path(ROOT.DIR,'../basis-projection/results/13_trait.RDS')
SPARSE_BASIS_PROJECTION_FILE <- file.path(ROOT.DIR,'../basis-projection/results/13_trait-sparse.RDS')

pc.emp <- readRDS(BASIS_FILE)
man.DT <- readRDS(SNP_MANIFEST_FILE)
sDT <- readRDS(SHRINKAGE_FILE)

all.filez <- list.files(path=DATA_DIR,pattern='*_source.RDS',full.names=TRUE)
proj <- mclapply(all.filez,function(f){
  message(f)
  trait <- basename(f) %>% gsub("\\_source.RDS","",.)
  sprintf("Processing %s",trait) %>% message
  dat <- readRDS(f)[,shrinkage:=NULL][]
  x<-cupcake::project_basis(dat,sDT,pc.emp,trait)
  ## optional save file with shrinkages for use with fdr approach for novel associations
  #saveRDS(x$data,file=file.path(OUT_DIR,sprintf("%s_source.RDS",trait)))
  ## optional save file with missing SNPs for comparison 
  #saveRDS(x$missing,file=file.path(OUT_DIR,sprintf("%s_missing.RDS",trait)))
  x$proj
},mc.cores=8)
full.proj <- do.call('rbind',proj)
saveRDS(full.proj,FULL_BASIS_PROJECTION_FILE)

## project sparse data

SPARSE_BASIS <- file.path(ROOT.DIR,'support/13_trait_basis-sparse.RData')
(load(SPARSE_BASIS))
library(parallel)
sparse.proj <- mclapply(all.filez,function(f){
  trait <- basename(f) %>% gsub("\\_source.RDS","",.)
  sprintf("Processing %s",trait) %>% message
  dat <- readRDS(f)[pid %in% rownames(rot.pca),]
  if(!any(names(dat)=='beta')){
     if(any(names(dat)=='or')){
        dat[,beta:=log(or)]
     }else{
        stop("GWAS input must have either an or or beta column")
     }
  }
  ## remove instances where beta is infinite or missing
  dat <- dat[!is.na(beta) | is.finite(beta),]
  if(!any(colnames(dat)=='seb'))
    dat[,seb:=abs(beta/qnorm(as.numeric(p.value)/2,lower.tail=FALSE))]
  cupcake::project_sparse(beta=dat$beta,seb=dat$seb,pids=dat$pid)[,trait:=trait][]
},mc.cores=8)

sparse.proj <- rbindlist(sparse.proj)
saveRDS(sparse.proj,SPARSE_BASIS_PROJECTION_FILE)

if(FALSE){
  ## compare projections for both methods - possible supp figure
  full.proj.delta <- (t(readRDS(FULL_BASIS_PROJECTION_FILE)) - pc.emp$x['control',]) %>% t
  full.proj.delta <- data.table(trait=rownames(full.proj.delta),full.proj.delta)
  full.proj.delta <- melt(full.proj.delta,id.vars='trait')
  sparse.proj <- readRDS(SPARSE_BASIS_PROJECTION_FILE)
  M <- merge(full.proj.delta[,.(trait,PC=variable,full.proj.delta=value)],sparse.proj[,.(trait,PC,sparse.proj=delta)],by=c('trait','PC'))
  par(mfrow=c(4,4))
  for(pc in paste('PC',1:13)){
    plot(M$full.proj.delta,M$sparse.proj,main=pc,xlab="Full basis delta",ylab="Sparse basis delta",cex=0.3,pch=16)
    abline(a=0,b=1,col='red',lty=2)
  }
}
