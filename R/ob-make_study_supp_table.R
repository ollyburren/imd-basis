library(xlsx)
ROOT.DIR <- '~/share/as_basis/basis-creation'
OUT_FILE <- '/home/ob219/share/as_basis/supp_tables/supplementary_tables_1_4_v4.xlsx'

## supp table 1 - basis traits
TRAIT_MANIFEST_FILE <- file.path(ROOT.DIR,'support/13_trait_manifest.tab')
trait.manifest <- fread(TRAIT_MANIFEST_FILE)
trait.manifest[,c('dtr','first_author'):=tstrsplit(disease,'_')]
trait.manifest <- trait.manifest[order(cases+controls),]
basis.out <- trait.manifest[,.(Trait=dtr,`First Author`=first_author,Reference=pmid,N0=controls,N1=cases)]
write.xlsx(basis.out,file=OUT_FILE,sheet='Basis',row.names=FALSE)

RESULTS.FILE <- '/home/ob219/share/as_basis/GWAS/RESULTS/30_08_19_0619_summary_results.RDS'
res.DT <- readRDS(RESULTS.FILE)[,.(trait,category,n0,n1,sdy)] %>% unique
keep.cat <- c('bb_disease','bb_cancer')
neale.out <- res.DT[category %in% keep.cat,]
neale.out[,c('pre','suf'):=tstrsplit(trait,':')]
neale.out[,trait:=suf]
neale.out <- neale.out[order(n1),.(Trait=trait,N0=n0,N1=n1),by=category]
neale.out[category=='bb_disease',category:='Self-reported disease']
neale.out[category=='bb_cancer',category:='Self-reported cancer']

imd <-
"asthma
type.1.diabetes
hypothyroidism.myxoedema
hyperthyroidism.thyrotoxicosis
adrenocortical.insufficiency.addison.s.disease
acute.infective.polyneuritis.guillain.barre.syndrome
multiple.sclerosis
ankylosing.spondylitis
myositis.myopathy
sarcoidosis
vasculitis
allergy.to.house.dust.mite
allergy.or.anaphylactic.reaction.to.food
allergy.hypersensitivity.anaphylaxis
allergy.or.anaphylactic.reaction.to.drug
systemic.lupus.erythematosis.sle
sjogren.s.syndrome.sicca.syndrome
scleroderma.systemic.sclerosis
hayfever.allergic.rhinitis
myasthenia.gravis
eczema.dermatitis
psoriasis
malabsorption.coeliac.disease
colitis.not.crohns.or.ulcerative.colitis
ulcerative.colitis
rheumatoid.arthritis
psoriatic.arthropathy
vitiligo
scleroderma.systemic.sclerosis
pernicious.anaemia
crohns.disease" %>% strsplit(.,"\n") %>% unlist
neale.out[,imd.trait:=0][Trait %in% imd,imd.trait:=1]
neale.out

write.xlsx(neale.out,file=OUT_FILE,sheet='UKBB Neale',row.names=FALSE,append=TRUE)


ga.out <- res.DT[category %in% c('geneatlas_srd'),]
ga.out[,c('pre','suf'):=tstrsplit(trait,':')]
ga.out[,trait:=suf]
ga.out <- ga.out[order(n1),.(Trait=trait,N0=n0,N1=n1),by='category']
ga.out[,category:='Self-reported disease']
ga.out[,imd.trait:=0][Trait %in% imd,imd.trait:=1]
write.xlsx(ga.out,file=OUT_FILE,sheet='UKBB GeneATLAS',row.names=FALSE,append=TRUE)

keep.cat <-
"bowes_jia_2019
psyc_consortium
myogen
lyons_egpa
estrada_NMO
wong_aav
bowes_psa
taylor_mtx
brown_as
kuiper_bs
li_as
psa_aterido
hasnoot_uveitis_jia" %>% strsplit(.,"\n") %>% unlist
rare.out <- res.DT[category %in% keep.cat,]
rare.out[category=='estrada_NMO',trait:=sprintf("%s_ssimp",trait)]
rare.out[trait=='na_psa',trait:='na_psa_ssimp']
rare.out[trait=='span_psa',trait:='span_psa_ssimp']
rare.out[trait=='li_as',trait:='li_ankspond']
rare.out[trait=='anca_Neg',trait:='anca.negative.egpa_lmm']
rare.out[trait=='mpo_Pos',trait:='mpo.anca.positive.egpa_lmm']
rare.out[trait=='egpa',trait:='all.egpa_lmm']



## show number of sparse basis SNPs missing
SPARSE_BASIS <- file.path(ROOT.DIR,'support/13_trait_basis-sparse.RData')
sparse.basis <- load(SPARSE_BASIS)
fs <- list.files(path='/home/ob219/share/as_basis/basis-projection/project_13',pattern='*missing.RDS',full.names=TRUE)
miss.dt <- lapply(fs,function(f){
  missing <- readRDS(f)
  sparse.missing <- sum(missing %in% rownames(rot.pca))
  message(sparse.missing)
  data.table(trait=basename(f) %>% gsub("\\_missing.RDS","",.),sparse.missing)
}) %>% rbindlist


setwd("~/git/imd-basis/")
source("./R/cw-reader.R")
other.traits <- reader()[fdrcat=="general" & category!="UKBB",.(category,trait,trait.label,category.label)]  %>% unique()
rare.out <- merge(other.traits,rare.out[,.(trait,n0,n1)],by='trait',all.x=TRUE)
rare.out <- merge(rare.out,miss.dt,by.x='trait',by.y='trait',all.x=TRUE)
## add missing mg data
rare.out[trait=='Overall_ssimp',c('n0','n1'):=list(1977,972)]
rare.out[trait=='LateOnset_ssimp',c('n0','n1'):=list(1977,737)]
rare.out[trait=='YoungOnset_ssimp',c('n0','n1'):=list(1977,235)]



rare.out[category=='ahola-olli_cytokine',trait:=gsub("CK:","",trait)]
rare.out[,c('fa','tr'):=tstrsplit(category,'_')]
rare.out <- rare.out[,.(Trait=trait,`First Author`=fa,Reference='unpublished',N0=n0,N1=n1,sdY=sdy,missing)]
rare.out[,Trait:=trait]

rare.out[grep('^jia',trait),Reference:='unpublished']
rare.out[Trait=='attention_deficit_hyperactivity_disorder',Reference:='29325848']
rare.out[Reference=='29325848',`First Author`:='martin']
rare.out[trait=='BIP',Reference:='29906448']
rare.out[trait=='SCZ',Reference:='29906448']
rare.out[Trait=='ADHD',Reference:='29325848']
rare.out[Reference=='29906448',`First Author`:='psyc_consortium']
rare.out[grep('NMO',Trait),Reference:='29769526']
rare.out[category=='EGPA',Reference:='https://doi.org/10.1101/491837']
rare.out[category=='Vasculitis',Reference:='22808956']
rare.out[Trait=='li_ankspond',Reference:='30946743']
rare.out[Trait=='na_psa_ssimp',Reference:='30552173']
rare.out[Trait=='span_psa_ssimp',Reference:='30552173']
rare.out[Trait=='hasnoot_uveitis_jia',Reference:='29513936']
#rare.out[`First Author`=='ahola-olli',Reference:='27989323']
rare.out[Trait=='birdshot_retinopathy',Reference:='24957906']
rare.out[Trait=='methotrexate',Reference:='29795407']
rare.out[category=='Myasthenia gravis',Reference:='25643325']
rare.out[Trait=='ank_spond',Reference:='21743469']
rare.out[category %in% c('PsA','Myositis'),Reference:='unpublished']
rare.out[,c('Trait','First Author'):=list(NULL,NULL)]
rare.out[,ssimp:=grepl('ssimp',trait)]
rare.out <- rare.out[,.(trait=trait.label,category=category.label,n0,n1,missing.basis.snps=sparse.missing,ssimp,reference=Reference)]
rare.out <- rare.out[order(n1),.(trait,category,n0,n1,missing.basis.snps,ssimp,reference),by=category]
write.xlsx(rare.out,file=OUT_FILE,sheet='Other',row.names=FALSE,append=TRUE)

## astle traits
astle <- res.DT[category=='astle_blood',.(trait,n=n0)]
astle <- astle[order(n),]
write.xlsx(res.DT[category=='astle_blood',.(trait,n=n0)],file=OUT_FILE,sheet='Astle Blood',row.names=FALSE,append=TRUE)

## roederer

roed <- res.DT[category=='roederer_immunophenotypes',.(trait=gsub("roederer\\_","",trait),n=n0)]
roed <- roed[order(n),]
write.xlsx(roed,file=OUT_FILE,sheet='Roederer Immunophenotype',row.names=FALSE,append=TRUE)

## ahola-olli_cytokine
ahola <- res.DT[category=='ahola-olli_cytokine',.(trait=gsub("CK:","",trait),n=n0)]
ahola <- ahola[order(n),]
write.xlsx(ahola,file=OUT_FILE,sheet='Ahola-olli Cytokines',row.names=FALSE,append=TRUE)
