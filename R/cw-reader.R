source("R/cw-files.R")
source("R/cw-utils.R")
source("R/cw-renamer.R")

## traits not used in this paper (unpublished/unimputed data)
.traits.drop=c("mpo_gwas2","mpo", "mpo_meta","pr3_gwas2","pr3_meta", # unpublished
              "dm_myogen", "jdm_myogen", "myositis_myogen",  "NMO_combined", "NMO_IgGNeg",  "NMO_IgGPos", "pm_myogen", "renton_mg", "renton_mg_early", "renton_mg_late", # ssimp exists
              "anca_Neg",  "egpa", "mpo_Pos", # lmm exists
              "CD_prognosis" # uncorrected
              )
## core traits from Atle paper
.traits.astle <- c('pdw','mpv','plt','irf','ret','rdw','hct','mch','mono','baso','eo','neut','lymph')  %>% paste0("ASTLE:",.)

reader <- function(what=c("sparse","fullfat")) {
    what <- match.arg(what)
    proj <- readRDS(switch(what,
                           fullfat="~/share/as_basis/basis-projection/results/13_traits.RDS",
                           sparse=#"~/share/as_basis/basis-projection/results/13_trait-sparse_with_meta.RDS"
                             "~/share/as_basis/basis-projection/results/13_trait-sparsev2.RDS"
                           ))
    ## attempt to add categories
    meta <- readRDS('~/share/as_basis/basis-projection/support/projection_meta.RDS')
    proj <- merge(proj,meta,by="trait",all.x=TRUE)
    ## mangle labels for consistency
    proj[grep("^roederer_",trait),c("category","category.label","trait.label"):=
                                    list("Roederer","Flow cytom",sub("roederer:","",trait))]
    proj[grep("^ASTLE:",trait),c("category","category.label","trait.label"):=
                                 list("ASTLE","Blood counts",sub("ASTLE:","",trait))]
    proj[,trait.label:=gsub(" ssimp|NMO|EGPA|AAV","",trait.label)]
    proj[trait.label=="systemic lupus erythematosis sle",trait.label:="SLE"]
    proj[trait=="li_ankspond",trait.label:="Turkish/Iranian"]
    proj[trait.label=="colitis not crohns or ulcerative colitis",trait.label:="colitis not Crohns or UC"]
    proj[category=="Vasculitis",trait.label:=sub("  1","",trait.label)]
    proj[,category.label:=ifelse(trait.label==category,"",category)]
    ## drop some things not for publication
    proj <- proj[category %notin% c("DELETE","tachmazidou_osteo","asthma") &
                 trait %notin% .traits.drop]
    # trim to Astle core traits
    proj <- proj[category!="Astle" | trait %in% .traits.astle]
    ## badly named by me
    setnames(proj,"p","p.value",skip_absent=TRUE)
    ## fdr on final dataset
    proj[,fdrcat:=ifelse(category %in% c("Geneatlas","Geneatlas_Cancer","Cytokines","Geneatlas_ICD","ASTLE","Roederer"),
                         category, "general")]
    proj[,newfdr:=p.adjust(p.value,method="BH"),by=c("fdrcat","PC")]
    proj[,fdr.overall:=p.adjust(p.overall,method="BH"),by=c("fdrcat","PC")]
    proj
}

## all.filez <- list.files(path=DATA_DIR,pattern='*.RDS',full.names=TRUE)
read_raw <- function(trait,pids=NULL) {
    if(length(trait) > 1)
        return(rbindlist(lapply(trait, read_raw, pids=pids), use.names=TRUE))
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

read_basis <- function(what=c("weighted","unweights"),delta=TRUE) {
    what <- match.arg(what)
    basis <- readRDS(switch(what,
                            weighted=BASIS_FILE,
                            unweighted=NOWEIGHT_BASIS_FILE))
    ## calc delta
    if(delta) {
        basis <- basis$x
        n <- nrow(basis)
        basis <- basis[1:(n-1),] - basis[rep(n,n-1),]
    }
    
    basis
}

