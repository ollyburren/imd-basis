## code to create and store basis and project all summary statistics

library(devtools)
#install_github('ollyburren/cupcake')
load_all('~/git/cupcake')
library(cupcake)
library(ggplot2)
library(cowplot)
library(ggrepel)


library(devtools)
install_github('ollyburren/cupcake@v0.1.0.0')
ROOT.DIR <- '~/share/as_basis/basis-creation'
SNP_MANIFEST_FILE <- file.path(ROOT.DIR,'support/6_trait_dense_snp_manifest.RDS')
TRAIT_MANIFEST_FILE <- file.path(ROOT.DIR,'support/6_trait_dense_manifest.tab')
SHRINKAGE_FILE <- file.path(ROOT.DIR,'support/6_trait_dense_shrinkage.RDS')
BASIS_FILE <- file.path(ROOT.DIR,'support/6_trait_dense_basis.RDS')
GWAS_DATA_DIR <- file.path(ROOT.DIR,'/sum_stats_dense/')

gwas.DT<-get_gwas_data(TRAIT_MANIFEST_FILE,SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE)
## make a list of sle variants with p=0 as shrinkage cannot cope with these.
rm.pid <- gwas.DT[p.value==0,]$pid %>% unique
gwas.DT <- gwas.DT[!pid %in% rm.pid,]
## are there any pallindromic SNPs that could be ambiguous or those with MAF<1%
shrink.DT<-compute_shrinkage_metrics(gwas.DT)
saveRDS(shrink.DT,file=SHRINKAGE_FILE)
## create shrunk basis
dense_basis <- create_basis(gwas.DT,shrink.DT)
saveRDS(basis,file=BASIS_FILE)

## make a basis using the same gwas data but sparse SNP map used in main paper analysis
SPARSE_SNP_MANIFEST_FILE <- file.path(ROOT.DIR,'support/13_snp_manifest.RDS')
SPARSE_SHRINKAGE_FILE <- file.path(ROOT.DIR,'support/6_trait_snp-sparse_shrinkage.RDS')
SPARSE_BASIS_FILE <- file.path(ROOT.DIR,'support/6_trait_snp-sparse_basis.RDS')
gwas.DT<-get_gwas_data(TRAIT_MANIFEST_FILE,SPARSE_SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE)
gwas.DT <- gwas.DT[!pid %in% rm.pid,]
shrink.DT<-compute_shrinkage_metrics(gwas.DT)
saveRDS(shrink.DT,file=SPARSE_SHRINKAGE_FILE)
## create shrunk basis
sparse_snp_basis <- create_basis(gwas.DT,shrink.DT)
saveRDS(basis,file=SPARSE_BASIS_FILE)


library("cowplot")
library("ggrepel")

## plot both scree and biplot

sbi <- function(pc.emp){
  vexp <- summary(pc.emp)[['importance']][2,]
  PC1.var<-signif(vexp["PC1"]*100,digits=3)
  PC2.var<-signif(vexp["PC2"]*100,digits=3)
  M <- cbind(as.data.table(pc.emp$x),trait=rownames(pc.emp$x))
  scp <- cbind(data.table(vexp),pcs=factor(names(vexp),levels=names(vexp)))
  scp[,cs:=cumsum(vexp)]
  scp$group=1
  text.size <- 20
  ## do a scree plot
  ppl <- ggplot(scp,aes(x=pcs,y=vexp,group=group)) + geom_point() + geom_line() + ylab("Variance Explained") + xlab("Principal Components") +
  theme(axis.text=element_text(size=text.size),axis.title=element_text(size=text.size),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  ppr<-ggplot(M,aes(x=PC1,y=PC2,label=trait)) + geom_point(size=3) + geom_text_repel(size=7) + # hjust = 0, nudge_x = 0.005)  +
  scale_color_discrete(guide=FALSE) + scale_alpha_discrete(guide=FALSE,range=c(0.3,1)) +
  xlab(sprintf("%s (%.1f%%)",'PC1',PC1.var)) + ylab(sprintf("%s (%.1f%%)",'PC2',PC2.var)) +  background_grid(major = "xy", minor = "none") +
  theme(axis.text=element_text(size=text.size),axis.title=element_text(size=text.size))
  list(ppl,ppr)
  plot_grid(ppl, ppr)
}

dense <- sbi(dense_basis)
sparse <- sbi(sparse_snp_basis)
plot_grid(dense,sparse,labels=c('a','b'),ncol=1)


## plot PC scores vs one another ?
small.DT <- data.table(trait=rownames(sparse_snp_basis$x),sparse_snp_basis$x)
small.DT <- melt(small.DT,id.vars="trait")
setnames(small.DT,'value','small')
big.DT <- data.table(trait=rownames(dense_basis$x),dense_basis$x)
big.DT <- melt(big.DT,id.vars="trait")
setnames(big.DT,'value','big')
m.DT <- merge(small.DT,big.DT,by=c('trait','variable'))
flip<-m.DT[,list(cor(small,big)<0),by='variable'][V1==TRUE,]$variable
m.DT[variable %in% flip, big:=big * -1]

pp2 <- ggplot(m.DT[variable!='PC7'],aes(x=small,y=big,color=trait)) + geom_point() +
geom_abline(intercept=0,slope=1,col='red',linetype=2,alpha=0.3) +
facet_wrap(~variable,nrow=3) + xlab("PC score sparse basis") + ylab("PC score dense basis") +
theme(legend.position="bottom") + labs(color = "Disease") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
save_plot("~/tmp/basis_comparison.pdf",pp2,base_height=5)
