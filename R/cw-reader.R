source("R/cw-files.R")
OUT_DIR <- file.path(ROOT.DIR,'../basis-projection/project_13_new')

if(FALSE) {
    SPARSE_BASIS <- file.path(ROOT.DIR,'support/13_trait_basis-sparse.RData')
(load(SPARSE_BASIS))
library(parallel)

    all.filez <- list.files(path=DATA_DIR,pattern='*.RDS',full.names=TRUE)
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
  cupcake:::project_sparse(beta=dat$beta,seb=dat$seb,pids=dat$pid)[,trait:=sub("_source.RDS","",trait)]
},mc.cores=8)
    proj <- rbindlist(sparse.proj)
    saveRDS(proj,file="~/sparse.RDS")
}


source("R/cw-renamer.R")
traits.drop=c("mpo_gwas2","mpo_meta","pr3_gwas2","pr3_meta", # unpublished
"dm_myogen", "jdm_myogen", "myositis_myogen",  "NMO_combined", "NMO_IgGNeg",  "NMO_IgGPos", "pm_myogen", "renton_mg", "renton_mg_early", "renton_mg_late", # ssimp exists
"anca_Neg",  "egpa", "mpo_Pos", # lmm exists
"CD_prognosis"
##              ,"span_psa","na_psa" # imputed
##              ,"dm_myogen","jdm_myogen","pm_myogen" # imputed
##              ,"NMO_combined","NMO_IgGNeg","NMO_IgGPos" # imputed
              )
reader <- function(what=c("sparse","fullfat")) {
    what <- match.arg(what)
    proj <- readRDS(switch(what,
                           fullfat="~/share/as_basis/basis-projection/results/13_traits.RDS",
                           sparse="~/share/as_basis/basis-projection/results/13_trait-sparse_with_meta.RDS"
                                        #~/share/as_basis/basis-projection/results/13_trait-sparse.RDS"
                           ))
    ## attempt to add categories
    ## proj  %<>% rename.traits()  %<>% rename.cats()
    proj <- proj[!(category %in% c("DELETE","tachmazidou_osteo")) &
                 !(trait %in% traits.drop)]
    proj[,trait.label:=gsub(" ssimp|NMO|EGPA|AAV","",trait.label)]
    proj[,category.label:=ifelse(trait.label==category,"",category)]
    setnames(proj,"p","p.value",skip_absent=TRUE)
    proj[,fdrcat:=ifelse(category %in% c("Geneatlas","Geneatlas_Cancer","Cytokines","Geneatlas_ICD"),
                         category, "general")]
    proj[,newfdr:=p.adjust(p.value,method="BH"),by=c("fdrcat","PC")]
    proj[,fdr.overall:=p.adjust(p.overall,method="BH"),by=c("fdrcat","PC")]
    proj
}

## all.filez <- list.files(path=DATA_DIR,pattern='*.RDS',full.names=TRUE)
readraw <- function(trait,pids=NULL) {
    if(length(trait) > 1)
        return(rbindlist(lapply(trait, readraw, pids=pids), use.names=TRUE))
    cat(trait,"\t")
    f <- trait
    if(!grepl("_source.RDS",f))
        f  %<>%  paste0(.,"_source.RDS")
    f %<>% file.path(DATA_DIR,.)
    dat <- readRDS(f)
    if(is.character(dat$p.value))
        dat[,p.value:=as.numeric(p.value)]
    if(!is.null(pids))
        dat <- dat[pids,]
    if(!("beta" %in% names(dat)))
        dat[,beta:=log(or)]
    if(!("seb" %in% names(dat)))
        dat[,seb:=abs(beta/qnorm(p.value/2,lower.tail=FALSE))]
    dat
}

