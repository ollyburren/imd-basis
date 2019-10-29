library(data.table)
library(magrittr)
library(ggplot2)
source("R/cw-reader.R")

## proj <- res$proj[!(category %in% c("astle_blood","geneatlas_icd","geneatlas_srd","roederer_immunophenotypes","ahola-olli_cytokine"))]
proj <- reader() 
proj[,cat.orig:=category]

## calculate FDR
with(proj,table(category,fdrcat))
proj[,variance:=var.proj]
proj[,stars:=ifelse(newfdr<0.01,"*","")]
proj[,FDR.overall:=p.adjust(p.overall,method="BH"),by=c("PC","fdrcat")]

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

pdf("~/share/as_basis/figures/suppfig-forest-everything.pdf",height=15,width=12)
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
    print(p)
}
dev.off()

if(interactive())
    system("evince ~/share/as_basis/figures/suppfig-forest-everything.pdf &")

## if(!interactive())
##     q("no")

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
p.pc1 <- pc1 + 
  geom_point(aes(y=0.15,x=trait.label,fill=aab,pch=aab),size=3,
             data=subset(pc1$data,!is.na(pc1$data$aab))) +
  ## theme(legend.position = c(0.8,0.8),legend.justification=c(1,0),
  theme(legend.position = "top",legend.justification = "left")  +
      guides(linetype="none") +
  scale_shape_manual("AAB",breaks=c("No", "NP", "P"),
                     values=c(No=21,NP=22,P=23,"TRUE"=23,"FALSE"=22)) +
  scale_fill_manual("AAB",breaks=c("No","NP","P"),
                      values=c(No="white",NP="lightblue",P="royalblue"))
p.pc1
save(p.pc1, file="~/basis-pc1-forest.RData")
ggsave("figures/fig4-pc1.pdf",height=10,width=8)


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


cowplot::plot_grid(p.pc1,NULL,p.pc3,nrow=1,
                   labels=c("a","","b"),
                   rel_widths=c(1,0.05,1))
ggsave("figures/fig4-pc1-pc3.pdf",height=8,width=10,scale=1.4)

cowplot::plot_grid(p.pc1,p.pc3,p.pc13,nrow=1)
ggsave("figures/fig4-pc1-pc3-pc13.pdf",height=8,width=12)


################################################################################

## junk below here

################################################################################

q()




## general heatmap of basis
library(pheatmap)
m=dcast(basis.DT,trait~PC,value.var="delta")
m2=as.matrix(m[,-1])
rownames(m2) <- m$trait
pheatmap(m2)

sigs=proj[FDR<0.01]$trait  %>% unique()
p=dcast(proj[PC %in% colnames(m2) & trait %in% sigs],trait~PC,value.var="delta")
p2=as.matrix(p[,-1])
rownames(p2) <- p$trait
pheatmap(p2)

t2=rbind(m2,p2)


pheatmap(t11,file="~/11.pdf",width=8,height=10,fontsize=8)
pheatmap(t13,file="~/13.pdf",width=8,height=10,fontsize=8)

ggplot(basis.DT,aes(y=trait,x=PC,fill=delta)) + geom_tile() +
    scale_fill_distiller(type="div")

## bigr <- function(p,size.rel = 1) {
##   p +
##     theme(
##       plot.title    = element_text(size = rel(2.0 * size.rel), hjust = 0),
##       plot.subtitle = element_text(size = rel(1.5 * size.rel), hjust = 0),
##       legend.title  = element_text(size = rel(1.8 * size.rel)),
##       legend.text   = element_text(size = rel(1.3 * size.rel)),
##       axis.title    = element_text(size = rel(1.5 * size.rel)),
##       axis.text     = element_text(size = rel(1.2 * size.rel)),
##       strip.text.x  = element_text(size = rel(1.8 * size.rel)),
##       strip.text.y  = element_text(size = rel(1.8 * size.rel)),
##       legend.key.height = unit(1.3 * size.rel, "line"),
##       legend.key.width  = unit(1.3 * size.rel, "line")
##     )
## }
################################################################################

## by cat group
sigcat.any <- unique(proj[(FDR<0.01 & category!="Basis" & category!="UKBB" & category!="cancer")]$category.label)  %>% setdiff(.,c("Infections","UKBB","Psychiatric"))
for(cat in sigcat.any) {
    pcs <- unique(proj[category.label==cat & FDR<0.01]$PC)  %>% as.character()
    traits <- unique(proj[category.label %in% c(cat,"Basis") |
                            (category.label=="UKBB" & PC %in% pcs & FDR < 0.01)]$trait)
    pdf(paste0("~/forest-",cat,".pdf"),height=15,width=12)
    forest_manypc(proj[trait %in% c(traits),],#,traits.i),],
                  basis.DT,pc=pcs,
                  focal=traits)
    print(plot_grid(plotlist=plots))
    dev.off()
}

proj[PC=="PC1" & category.label=="JIA",.(trait.label,n1)]

## not sig
proj[FDR<0.01]$trait.label->sigs
nosigs <- setdiff(proj[!(category.label %in% c("UKBB","Infections"))]$trait.label,sigs)
unique(proj[trait.label %in% nosigs,.(category.label,trait.label,n0,n1)])
unique(proj[category.label!="UKBB" & trait.label %in% sigs,.(category.label,trait.label,n1)][order(n1)])

for(i in 1:11) {
    ipc=paste0("pc",i)
    iPC=paste0("PC",i)
    sigcat2  %<>%  setdiff(.,c("Infections","ahola−olli_cytokine","","Psychiatric"))
    traits <- unique(proj[category.label %in% c(sigcat,"Basis") | (PC==iPC & FDR < 0.01)]$trait)
    traits  %<>%  c(.,nm)  %>% unique()
    ## AS custom
    if("AS (Brown)" %in% traits)
        traits  %<>% c("ankylosing spondylitis",.)  %>% unique
    if( "ankylosing spondylitis" %in% traits)
        traits  %<>% c("AS (Brown)",.)  %>% unique
    proj.dat=proj[trait %in% c(traits),]
    basis.dat=basis.DT
    pc=iPC
    focal=traits
    fdr_thresh=0.01
    ## p <- forest_everything(proj[trait %in% c(traits),],#,traits.i),],
    p <- forest_labelled(proj[trait %in% c(traits),],#,traits.i),],
                     basis.DT,pc=iPC,
                     focal=traits)
    print(p)
}
dev.off()
    


## individual customized
forest_plot(proj.dat=proj[category %in% c("UKBB/GWAS"),],
            basis.dat=basis.DT,
            pc="PC1",
            title="")  #%>% bigr()

## MYGEN, PC1
cats.pub=c("bb_disease","myogen","brown_as","cousminer_lada","ferreira_asthma","estrada_NMO","li_as","tian_infectious_disease","psyc_consortium")
forest_plot(proj.dat=proj[cat.orig %in% cats.pub],
                  basis.DT,pc="PC1",focal=all.traits[c('MYOGEN')])
ggsave("~/myogen-pc1.pdf",height=10,width=8)
forest_plot(proj.dat=proj[cat.orig %in% cats.pub],
                  basis.DT,pc="PC11",focal=all.traits[c('MYOGEN')])
ggsave("~/myogen-pc11.pdf",height=10,width=8)

## AS, PC1
forest_plot(proj.dat=proj[category %in% c("UKBB/GWAS","ank.spond")],
                  basis.DT,pc="PC1",focal=all.traits[c('ank.spond')])
## JIA, PC1, PC3
for(i in c(1,3,4)) {
forest_plot(proj.dat=proj[cat.orig %in% c(cats.pub,"bowes_jia_2019")],
            basis.DT,pc="PC3",focal=all.traits[c('JIA')])
forest_plot(proj.dat=proj[cat.orig %in% c(cats.pub,"bowes_jia_2019")],
            basis.DT,pc="PC3",focal=all.traits[c('JIA')])
}

################################################################################

## pairs
proj=copy(bak)

library(ggforce,quietly=TRUE)
t1="type 1 diabetes"
t2="cousminer lada"
pairplot <- function(t1,t2) {
    ## dat[is.na(variance),variance:=0]
    d=merge(proj[trait==t1],proj[trait==t2],by="PC",suffixes=c(".1",".2"))
    ggplot(d,aes(x=delta.1,y=delta.2,label=PC,col=PC)) +
      geom_abline(colour="grey") +
      geom_hline(yintercept = 0,colour="grey") +
      geom_vline(xintercept = 0,colour="grey") +
      geom_ellipse(aes(x0=delta.1,y0=delta.2,a=1.96*sqrt(variance.1),b=1.96*sqrt(variance.2),angle=0),lty=2) +
      geom_point() +
      geom_label() +
      labs(x=t1,y=t2)
}


library(coloc)
paircoloc <- function(t1,t2) {
   ## dat[is.na(variance),variance:=0]
    d=merge(proj[trait==t1],proj[trait==t2],by="PC",suffixes=c(".1",".2"))
    ct=coloc.test.summary(b1=d$delta.1,b2=d$delta.2,
                          V1=diag(d$variance.1),
                          V2=diag(d$variance.2))
    p=plot(ct) + geom_abline(slope=ct@result["eta.hat"]) + labs(x=t1,y=t2)
    print(ct)
    print(p)
}
     
grep("AS",proj$trait,value=TRUE)
pairplot("AS (Brown)","AS (Li)")
paircoloc("AS (Brown)","AS (Li)")

pairplot("type 1 diabetes","cousminer lada")
paircoloc("type 1 diabetes","cousminer lada")
paircoloc("type 1 diabetes","AS (Brown)")
paircoloc("na psa","span psa")
paircoloc("bowes psa","na psa")
paircoloc("bowes psa","span psa")
