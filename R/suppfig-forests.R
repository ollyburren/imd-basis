library(data.table)
library(magrittr)
library(ggplot2)
source("R/cw-reader.R")

## proj <- res$proj[!(category %in% c("astle_blood","geneatlas_icd","geneatlas_srd","roederer_immunophenotypes","ahola-olli_cytokine"))]
proj <- reader() 
proj[,cat.orig:=category]

## calculate FDR
with(proj,table(category,fdrcat))
proj[,stars:=ifelse(newfdr<0.01,"*","")]
proj[,FDR.overall:=fdr.overall]
proj[,variance:=var.proj]

## rename
proj[is.na(trait.label),] # should be empty
proj[is.na(category.label),] # should be empty
head(proj[grepl("UKBB",category),])
proj[category.label=="Myasthenia gravis",category.label:="Myasthenia\ngravis"]

basis <- read_basis()

pc1.traits <- c(proj[newfdr<0.01 & PC=="PC1" & fdrcat=="general"]$category.label  %>% unique(),
                proj[newfdr<0.01 & PC=="PC1" & category.label=="UKBB"]$trait.label  %>% unique(),
                rownames(basis))  %>% unique()  %>%
  setdiff(., c("control","UKBB"))  %>% sort()
cat(pc1.traits,sep="\n")

## add in basis traits for comparison
basis.DT <- data.table(trait=rownames(basis),basis) %>%
  melt(.,id.vars='trait',value.name="delta",variable.name="PC")
basis.DT[,category:="basis"]
## rename to match projected data
nm <- c(asthma="asthma",
        CD="crohns disease",
        UC="ulcerative colitis",
        PSC="primary sclerosing cholangitis",
        CEL="malabsorption coeliac disease",
        MS="multiple sclerosis",
        RA="rheumatoid arthritis",
        PBC="primary billiary cholangitis",
        T1D="type 1 diabetes",
        VIT="vitiligo",
        SLE="SLE",
        "cousminer_lada"= "cousminer_lada",
        "kiryluk_neph" = "kiryluk_neph" )  
for(i in seq_along(nm))
    basis.DT[trait==names(nm)[i], trait:=nm[i]]
basis.DT[,category.label:=category]
basis.DT  %<>% rename.traits() #[,trait.label:=trait]
all.traits <- traits<-split(proj$trait,proj$category) %>% lapply(.,unique)

## plot

btraits <- unique(basis.DT$trait)
source("R/cw-forests.R")
theme_set(theme_cowplot(16))
theme_update(plot.title = element_text(hjust = 0))

bak=copy(proj)

sigtraits <- function(iPC) {
    sigcat <- unique(proj[PC==iPC # right PC
        & FDR.overall<0.01 # overall significant
        & newfdr<0.01 # this PC significant
        ]$category.label)
    ## but will do UKBB, characterising datasets separately
    sigcat  %<>%  setdiff(.,c("UKBB","",NA,"Blood counts","Cytokines"))
    traits <- unique(proj[category.label %in% sigcat # signif cat
                          | trait.label %in% nm # basis traits
                            | (PC==iPC & newfdr < 0.01 & FDR.overall<0.01) # signif trait
                            | (category.label %in% c("Blood counts","Cytokines") & PC==iPC & newfdr < 0.01) # signif characterising trait
                            ]$trait) 
    traits  %<>%  c(.,names(nm))  %>% unique()  # add basis traits
}

## all basic 
## drop GA
proj <- bak[fdrcat=="general"]
#!(category %in% c("astle_blood","geneatlas","geneatlas_icd","geneatlas_srd","roederer_immunophenotypes","ahola-olli_cytokine"))]
## setnames(proj,c("trait.label","category.label"),c("trait","category"))
proj[category.label=="Myasthenia gravis", category.label:="Myasthenia\ngravis"]

## pdf("~/share/as_basis/figures/suppfig-forest-everything.pdf",height=15,width=12)
for(i in 1:13) {
    ipc=paste0("pc",i)
    iPC=paste0("PC",i)
    traits <- sigtraits(iPC)
    ## to debug
    proj.dat=proj[trait %in% c(traits),]; basis.dat=basis.DT; pc=iPC; focal=traits; fdr_thresh=0.01
    p <- forest_labelled(proj[trait %in% c(traits),],
                     basis.DT,pc=iPC,
                     focal=traits,
                     order.within=TRUE)
    ggsave(paste0("figures/suppfig-forest-pc",i,".pdf"),plot=p,height=10,width=8,scale=1.5)
}
## dev.off()

if(interactive())
    system("evince ~/share/as_basis/figures/suppfig-forest-pc1.pdf &")

################################################################################

## specific PC plots

aab <- read.table(text="trait,Pathogenic auto-antibodies,Any autoantibodies,
ankylosing spondylitis,No,No,
asthma,No,No,
basal cell carcinoma,No,No,
birdshot retinopathy,No,No,
coeliac disease,No,Yes,
colitis,No,No,
Latent Autoimmune Diabetes in Adults,Yes,-,
crohns disease,No,No, (some anti-commensal antibodies but not truly auto and not pathogenic)
diabetes,No,No, (assuming mostly type 2)
eczema dermatitis,No,No,
gastroenteritis dysentry,No,No,
hayfever allergic rhinitis,No,No,
hyperthyroidism thyrotoxicosis,Yes,-,
hypothyroidism myxoedema,Yes,-,
juvenile idiopathic arthritis,No,No,
IgA nephropathy,No,No,no true autoantibodies to my knowledge (Ken will know better)
malignant melanoma,No,No,
MS,No,No,
Myasthenia gravis,Yes,-,
Myositis,Yes,-,
neuro myolitis optica,Yes in seropositive pts,Yes in seropositive pts,
IgGPos_ssimp,Yes,Yes,
IgGNeg_ssimp,No,No,
PBC,No,Yes,
pernicious anaemia,Yes,Yes, - but only if caused by autoimmune gastritis
PSC,No,No,
psoriatic arthritis,No,No,
rheumatoid arthritis,Yes,-,
Sjogrens/sicca,Yes,-,
SLE,Yes,-,
type 1 diabetes,Yes,-,
ulcerative colitis,No,No,
vitiligo,No,Yes,",
sep=",",header=TRUE)

## Figure 4



i=1
tmp <- bak[!(category %in% c("ASTLE","geneatlas_icd","geneatlas_srd","roederer_immunophenotypes","ahola-olli_cytokine"))]
table(aab$trait %in% tmp$trait.label)
setdiff(aab$trait, tmp$trait.label)
tmp <- merge(tmp,aab,by.x="trait.label",by.y="trait",all.x=TRUE)
basis.DT[trait=="IgA_NEPH",trait:="IgA_nephropathy"]
basis.DT[trait.label=="IgA NEPH",trait.label:="IgA nephropathy"]
tmp[trait %in% c("multiple sclerosis","ank_spond","li_as") |
    trait.label %in% c("multiple sclerosis","ank_spond","li_as","colitis not crohns or ulcerative colitis","IgA nephropathy") |
    category.label %in% c("PsA","JIA","asthma","Ank.Spond"),
    c("Pathogenic.auto.antibodies","Any.autoantibodies"):=list("No","No")]
tmp[trait.label==" _combined" & category.label=="JIA",
    c("Pathogenic.auto.antibodies","Any.autoantibodies"):=list(NA,NA)]
tmp[trait.label %in% c("primary billiary cholangitis","malabsorption coeliac disease","primary sclerosing cholangitis") |
    category.label=="",
    c("Pathogenic.auto.antibodies","Any.autoantibodies"):=list("No","Yes")]
tmp[trait.label %in% c("RF+","systemic lupus erythematosis sle","SLE") |
    category.label %in% c("Myositis","Myasthenia gravis","Myasthenia\ngravis"),
    c("Pathogenic.auto.antibodies","Any.autoantibodies"):=list("Yes","-")]
tmp[,aab:=ifelse(Pathogenic.auto.antibodies=="Yes","P",
                 ifelse(Any.autoantibodies=="Yes","NP","No"))]
tmp[category.label=="Myasthenia gravis", category.label:="Myasthenia\ngravis"]

tmp <- rbind(tmp,
             data.table(trait.label=c("primary billiary cholangitis","primary sclerosing cholangitis","IgA nephropathy","LADA"),
                        category.label=c("UKBB"),
                        aab=c("No","NP","No","P")),
             fill=TRUE)
tmp.basis <- merge(basis.DT,tmp[,.(trait.label,aab)],by="trait.label",all.x=TRUE)


ipc=paste0("pc",i)
iPC=paste0("PC",i)
traits <- sigtraits(iPC)
proj.dat=tmp[trait %in% c(traits) & trait!="ASTLE:eo",]
    pc1 <-  forest_labelled(tmp[trait %in% c(traits),],#,traits.i),],
                     tmp.basis,pc="PC1",
                     focal=traits)
pc1$data$aab %<>% factor(., levels=c("No","NP","P"))
levels(pc1$data$aab) <- c("None","Non-pathogenic","Pathogenic")

theme_set(theme_cowplot(10))
p.pc1 <- pc1 + 
  geom_point(aes(y=0.15,x=trait.label,fill=aab,shape=aab),
             size=4,
             data=subset(pc1$data,!is.na(pc1$data$aab))) +
  theme(legend.position = "top",
        legend.justification = "right",
        legend.text=element_text(size=rel(0.8)),
        legend.title=element_text(size=rel(1)),
        plot.title=element_blank())  +
  scale_shape_manual("AAB",breaks=levels(pc1$data$aab),
                     values=c(None=21,"Non-pathogenic"=22,"Pathogenic"=23,"TRUE"=23,"FALSE"=22)) +
  scale_fill_manual("AAB",breaks=levels(pc1$data$aab),
                    values=c(None="white","Non-pathogenic"="lightblue","Pathogenic"="royalblue","TRUE"="white","FALSE"="white")) +
  guides(linetype="none",shape="legend",fill="legend") 
p.pc1

ggsave("figures/fig4-pc1.pdf",height=5,width=4,scale=1.8)

save(p.pc1, file="~/basis-pc1-forest.RData")

################################################################################
## PC3

proj <- bak[fdrcat=="general" | category=="Cytokines"]
i=3
iPC=paste0("PC",i)
traits <- sigtraits(iPC)  
proj.dat=proj[trait %in% c(traits),]
proj.dat$category.label  %<>% factor(., c("Ank.Spond", "JIA", "Myasthenia\ngravis", "PsA", "Cytokines", "UKBB"))
p.pc3 <-  forest_labelled(proj.dat, #proj[trait %in% c(traits),],#,traits.i),],
                     basis.DT,pc="PC3",
                     focal=traits)
p.pc3

save(p.pc3, file="~/basis-pc3-forest.RData")
ggsave("figures/fig4-pc3.pdf",height=10,width=8)

## cytokines for urs
## proj <- bak[!(category %in% c("astle_blood","geneatlas_icd","geneatlas_srd","roederer_immunophenotypes"))]
## i=3
## ipc=paste0("pc",i)
## iPC=paste0("PC",i)
## sigcat <- unique(proj[(PC==iPC & FDR<0.01 & category!="Basis" & category!="bb_disease" & category!="bb_cancer")]$category.label)
## sigcat  %<>%  setdiff(.,c("Infections","ahola−olli_cytokine","","Psychiatric"))
## traits <- unique(proj[category.label %in% c(sigcat,"Basis") | (PC==iPC & FDR < 0.01)]$trait)
## traits  %<>%  c(.,nm)  %>% unique()
## proj.dat=proj[trait %in% c(traits),]
## p.urs <-  forest_labelled(proj[trait %in% c(traits),],#,traits.i),],
##                      basis.DT,pc="PC3",
##                      focal=traits,
##                      order.within=TRUE) + scale_y_reversed()
## p.urs
## proj[trait %in% c(traits) & PC==iPC & FDR < 0.1 & category.label=="cytokines",]
## ggsave("pc3-cytokines.pdf",height=10,width=8,scale=1.5)

proj <- bak[fdrcat=="general" | (category=="ASTLE" & trait=="ASTLE:eo")]
proj[category=="ASTLE",category.label:="Blood counts"]
proj[trait.label=="eo",trait.label:="eosinophils"]
traits  <- sigtraits("PC13")  %>% setdiff(., c("baso","pdw","rdw","neut","plt","hct","mpv","ret","lymph","mono","irf","mch")) # remove inconsistent blood count traits
proj.dat=proj[trait %in% c(traits),]
    p.pc13 <-  forest_labelled(proj[trait %in% c(traits),],#,traits.i),],
                     basis.DT,pc="PC13",
                     focal=traits)
p.pc13
save(p.pc13, file="~/basis-pc13-forest.RData")
ggsave("figures/fig4-pc13.pdf",height=10,width=8)


pdf("figures/fig4-pc1-pc13.pdf",height=8*1.4,width=10*1.4)
cowplot::plot_grid(p.pc1,NULL,p.pc13,nrow=1,
                   labels=c("a","","b"),
                   rel_widths=c(1,0.05,1))
dev.off()
## ggsave("figures/fig4-pc1-pc13.pdf",height=8,width=10,scale=1.4)

cowplot::plot_grid(p.pc1,p.pc3,p.pc13,nrow=1)
ggsave("figures/fig4-pc1-pc3-pc13.pdf",height=8,width=12)

if(!interactive())
    q("no")


