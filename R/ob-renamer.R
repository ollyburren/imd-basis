source('~/cw-renamer.R')
ROOT.DIR <- '~/share/as_basis/basis-creation'
FULL_BASIS_PROJECTION_FILE <- file.path(ROOT.DIR,'../basis-projection/results/13_trait.RDS')
full<-readRDS(FULL_BASIS_PROJECTION_FILE)
fdt<-data.table(trait=rownames(full))
rename.traits(fdt)
fdt <- fdt[-grep("^UKBB_NEALE:SRM|^roederer_|^ASTLE|^TIAN",trait),]
rename.cats(fdt)
RESULTS.FILE <- '~/share/as_basis/GWAS/RESULTS/14_10_13_traits_0919_summary_results_with_seb_variance.RDS'
res.DT <- readRDS(RESULTS.FILE)[,.(trait,old.cat=category),] %>% unique
M <- merge(fdt,res.DT,by='trait',all.x=TRUE)
M[is.na(category.label),category.label:=old.cat]
M[grep('^NMO',trait),category.label:='NMO']
M[trait %in% c('anca_Neg','egpa','mpo_Pos'),category.label:='EGPA']
M[grep('^gwas\\_.*\\_new_ssimp',trait),category.label:='Myositis']
M[trait == 'mpo',category.label:='Vasculitis']
M[trait == 'li_ankspond',category.label:='Ank.Spond']
M[trait == 'renton_mg',category.label:='Myasthenia gravis']
M<-M[,.(trait,trait.label,category=category.label)]
rename.cats(M)
M[!is.na(category.label),category:=category.label]
M[category=="",category:=trait]
M[is.na(category),category:='DELETE']
M[category=='geneatlas_icd',category:='Geneatlas_ICD']
M[category=='geneatlas_cancer',category:='Geneatlas_Cancer']
saveRDS(M[,.(trait,trait.label,category)],'~/share/as_basis/basis-projection/support/projection_meta.RDS')
