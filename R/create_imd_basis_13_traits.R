## code to create and store basis note that this code uses cupcake v0.1.0.0

library(devtools)
install_github('ollyburren/cupcake@v0.1.0.0')
ROOT.DIR <- '~/share/as_basis/basis-creation'
SNP_MANIFEST_FILE <- file.path(ROOT.DIR,'support/13_snp_manifest.RDS')
TRAIT_MANIFEST_FILE <- file.path(ROOT.DIR,'support/13_trait_manifest.tab')
SHRINKAGE_FILE <- file.path(ROOT.DIR,'support/13_trait_shrinkage.RDS')
BASIS_FILE <- file.path(ROOT.DIR,'support/13_trait_basis.RDS')
NOWEIGHT_BASIS_FILE <- file.path(ROOT.DIR,'support/13_trait_basis-noweight.RDS')
GWAS_DATA_DIR <- file.path(ROOT.DIR,'/sum_stats/')

## get basis gwas data and annotate with ld blocks (HapMap) and maf from reference genotype set (1000 genomes EUR).
gwas.DT<-get_gwas_data(TRAIT_MANIFEST_FILE,SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE)
shrink.DT<-compute_shrinkage_metrics(gwas.DT)
saveRDS(shrink.DT,file=SHRINKAGE_FILE)
## create shrunk basis
basis <- create_basis(gwas.DT,shrink.DT)
saveRDS(basis,file=BASIS_FILE)
## create unweighted basis
noweight.basis <- create_basis(gwas.DT,shrink.DT,apply.shrinkage=FALSE)
saveRDS(noweight.basis,NOWEIGHT_BASIS_FILE)

## sanity checking plotting code
if(FALSE){
  library(cowplot)
  library(ggrepel)
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
    plot_grid(ppl, ppr, labels = "auto",label_size=text.size)
  }
  shrunk <- sbi(basis)
  unshrunk <- sbi(noweight.basis)
  plot_grid(shrunk,unshrunk)
}
