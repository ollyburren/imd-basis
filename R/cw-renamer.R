rename.traits.list <- read.table(text='original shiny short
CK:CTACK CTACK CTACK
CK:IFNg IFNg IFNg
CK:IL18 IL18 IL18
CK:IL5 IL5 IL5
CK:MCP1 MCP1 MCP1
CK:MIP1b MIP1b MIP1b
CK:TNFa TNFa TNFa
CK:Eotaxin Eotaxin Eotaxin
CK:IL10 IL10 IL10
CK:IL1b IL1b IL1b
CK:IL6 IL6 IL6
CK:MCP3 MCP3 MCP3
CK:PDGFbb PDGFbb PDGFbb
CK:TNFb TNFb TNFb
CK:FGFBasic FGFBasic FGFBasic
CK:IL12p70 IL12p70 IL12p70
CK:IL1ra IL1ra IL1ra
CK:IL7 IL7 IL7
CK:MCSF MCSF MCSF
CK:RANTES RANTES RANTES
CK:TRAIL TRAIL TRAIL
CK:HGF HGF HGF
CK:IL17 IL17 IL17
CK:IL4 IL4 IL4
CK:IP10 IP10 IP10
CK:MIP1a MIP1a MIP1a
CK:SDF1a SDF1a SDF1a
CK:GCSF GCSF GCSF
CK:IL13 IL13 IL13
CK:IL2 IL2 IL2
CK:IL8 IL8 IL8
CK:MIF MIF MIF
CK:SCF SCF SCF
CK:VEGF VEGF VEGF
CK:GROa GROa GROa
CK:IL16 IL16 IL16
CK:IL2ra IL2ra IL2ra
CK:IL9 IL9 IL9
CK:MIG MIG MIG
CK:SCGFb SCGFb SCGFb
CK:bNGF bNGF bNGF
ASTLE:eo eo eosinophils
jia_ERA_19 enthesitis_related_jia ERA
jia_EO_19 extended_oligo_jia EO
jia_PO_19 persistent_oligo_jia PO
jia_PsA_19 psoriatic_jia PsA
jia_RFneg_19 polyoligo_rf-_jia RF-
jia_RFpos_19 polyoligo_rf+_jia RF+
jia_undiff_19 undifferentiated_jia undiff.
jia_sys_19 systemic_jia systemic
jia_case_19 combined_jia _combined
CD_prognosis crohns_disease_prognosis Crohn\'s prognosis
ADHD attention_deficit_hyperactivity_disorder ADHD
BIP bipolar_disorder bipolar
SCZ schizophrenia schizophrenia
gwas_dmjdmpm_new_ssimp combined_myositis _combined
gwas_jdm_new_ssimp juvenile_dermatomyositis JDM
gwas_dm_new_ssimp dermatomyositis DM
gwas_pm_new_ssimp polymyositis PM
jdm_myogen juvenile_dermatomyositis JDM
dm_myogen dermatomyositis DM
pm_myogen polymyositis PM
myositis_myogen combined_myositis _combined
NMO_combined combined_neuromyelitis_optica _combined
NMO_IgGNeg igg-_neuromyelitis_optica IgG-
NMO_IgGPos igg+_neuromyelitis_optica IgG+
NMO_combined_ssimp combined_neuromyelitis_optica _combined
NMO_IgNeg_ssimp igg-_neuromyelitis_optica IgG-
NMO_IgPos_ssimp igg+_neuromyelitis_optica IgG+
anca_Neg anca-_egpa ANCA-_EGPA
anca.negative.egpa_lmm anca-_egpa ANCA-_EGPA
egpa combined_egpa _combined
all.egpa_lmm combined_egpa _combined
mpo mpo+_aav MPO+_AAV
mpo_Pos mpo+_egpa MPO+_EGPA
mpo.anca.positive.egpa_lmm mpo+_egpa MPO+_EGPA
mpo_gwas1 mpo+_aav_1 MPO+_AAV_1
mpo_gwas2 mpo+_aav_2 MPO+_AAV_2
mpo_meta mpo+_aav_meta MPO+_AAV_meta
pr3_gwas1 pr3+_aav_1 PR3+_AAV_1
pr3_gwas2 pr3+_aav_2 PR3+_AAV_2
pr3_meta pr3+_aav_meta PR3+_AAV_meta
bowes_psa psoriatic_arthritis_unpublished UK
li_ankspond ankylosing_spondylitis_li Turkish/Iranian
ank_spond ankylosing_spondylitis_tasc International
mahajan_t2d type_2_diabetes type_2<-delete
na_psa_ssimp psoriatic_arthritis_north_america N.American
span_psa_ssimp psoriatic_arthritis_spanish Spanish
hasnoot_uveitis_jia uveitis_jia Uveitis
LateOnset_ssimp renton_mg_late late onset
YoungOnset_ssimp renton_mg_early early onset
Overall_ssimp renton_mg_combined _combined
renton_mg_late renton_mg_late late onset
renton_mg_early renton_mg_early early onset
renton_mg_combined renton_mg_combined _combined
SRD.allergy.hypersensitivity.anaphylaxis allergy anaphylaxis allergy anaphylaxis
SRD.colitis.not.crohns.or.ulcerative.colitis colitis colitis
SRD.sjogren.s.syndrome.sicca.syndrome sjogrens.sicca Sjogren\'s/sicca
SRD.systemic.lupus.erythematosis.sle SLE SLE
SRD.chronic.obstructive.airways.disease.copd COPD COPD'
,header=TRUE,stringsAsFactors=FALSE)

rename.cats.list=read.csv(text='original short
ahola-olli_cytokine Cytokines
astle_blood Blood counts
asthma Allergy
disease UKBB
bb_disease UKBB
bb_cancer UKBB
iga IgA Nephropathy
kiryluk_iga_neph IgA Nephropathy
kiryluk_neph IgA Nephropathy
jia JIA
bowes_jia_2019 JIA
egpa EGPA
lyons_egpa EGPA
lyons_egpa_lmm EGPA
lyons_vasculitis Vasculitis
wong_aav Vasculitis
ferreira_asthma asthma
psyc_consortium Psychiatric
psa PsA
bowes_psa PsA
lee_CD_prognosis
taylor_mtx
kuiper_bs
aterido PsA
psa_aterido PsA
arterido_psa_ssimp PsA
lada Diabetes
cousminer_lada Diabetes
uveitis Uveitis
cytokine Cytokines
NMO NMO
estrada_NMO NMO
estrada_nmo_ssimp NMO
myositis Myositis
myogen Myositis
myogen_myositis_ssimp Myositis
as Ank.Spond
brown_as Ank.Spond
li_as Ank.Spond
ahola−olli_cytokine Cytokines
hasnoot_uveitis_jia Uveitis
bowes_jia_2019 JIA
tian_infectious_disease Infections
renton_mg_ssimp Myasthenia gravis
renton_mg Myasthenia gravis
GA Geneatlas
CK Cytokines
ASTLE Blood counts
UKBB_NEALE UKBB
',header=TRUE,stringsAsFactors=FALSE,sep= )

message(dt=rename.traits(dt) to rename traits)
rename.traits <- function(dt) {
    dt[,trait:=as.character(trait)]
    dt[,trait.label:=gsub(^roederer_|^.*:|^bb_|^SRD.|^SRM.|^SRC.,,trait)]
    renames <- subset(rename.traits.list, rename.traits.list$original %in% dt$trait)
    if(nrow(renames))
        for(i in 1:nrow(renames))
            dt[trait==renames$original[i],trait.label:=renames$short[i]]
    dt[,trait.label:=gsub(., ,trait.label,fixed=TRUE)]
    dt[,trait.label:=gsub(., ,trait.label,fixed=TRUE)]
    dt[,trait.label:=gsub(_, ,trait.label,fixed=TRUE)]
    dt[,trait.label:=gsub(combined,_combined,trait.label,fixed=TRUE)]
    dt
}

message(dt=rename.cats(dt) to rename cats)
rename.cats <- function(dt) {
    dt[,trait:=as.character(trait)]
    if(!(category %in% names(dt))) {
        dt[,category:=GWAS]
        dt[grep(:,trait),category:=sub(:.*,,trait)]
        dt[grep(^GA:[A-Z][0-9][0-9].,trait),category:=Geneatlas_ICD]
        dt[grep(^UKBB_NEALE:SRM:,trait),category:=SRM]
        dt[grep(^roederer_,trait),category:=Flow cytom]
    }
    dt[,category:=as.character(category)]
    renames <- subset(rename.cats.list, rename.cats.list$original %in% dt$category)
    if(nrow(renames))
        for(i in 1:nrow(renames))
            dt[category==renames$original[i],category.label:=renames$short[i]]
    dt
}
