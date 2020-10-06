library(data.table)
library(magrittr)
library(ggplot2)
library(cowplot)
source("R/cw-reader.R")

imd.trait <-
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
crohns.disease" %>% strsplit(.,"\n") %>% unlist()

source("R/cw-reader.R")
## reader()[fdrcat=="general" & category!="UKBB",.(category,trait)]  %>% unique()

proj <- reader()[category=="UKBB"]
proj[,imd:=sub(".*:","",trait) %in% imd.trait]
table(proj$imd)
length(unique(imd.trait)) * 13


fdr <- proj[PC=="PC1"][order(fdr.overall)]
dim(fdr)

library(gsheet)
N <- gsheet2tbl("https://docs.google.com/spreadsheets/d/1bej1XFq6WwvTlURxhRFmfad2X0BGWZXf7D6IlK6WN0Q/edit#gid=1742820075")  %>%
  as.data.table()

fdr[,Trait:=sub(".*:","",trait)]
fdr <- merge(fdr,N,by="Trait")
fdr <- fdr[!((trait=="UKBB_NEALE:SRD:unclassifiable" & N1 < 1000) |
             (trait=="UKBB_NEALE:SRC:unclassifiable" & N1 > 1000))] # distinguish two kinds of unclassifiable


library(ggrepel)

fdr <- fdr[order(N1)]
fdr[,trait.label:=factor(as.character(trait.label),levels=unique(fdr$trait.label))]
fdr[,sig:=fdr.overall<0.01]
mx <- max(fdr$N1)

fdr[,cat:=paste(ifelse(imd,"immune-related","not immune-related"),
                ifelse(sig,"sig","not sig"))]
                

ggplot(fdr,aes(x=trait.label,fill=imd,col=sig)) +
  ## geom_col(aes(y=mx),position="dodge",alpha=0.5) +
  geom_col(aes(y=N1),position="dodge") +
  geom_label_repel(aes(y=N1,label=trait.label), data=fdr[imd==TRUE | sig==TRUE]) +
  scale_fill_viridis_d(option="A") +
  scale_x_discrete(breaks=fdr[imd==TRUE]$trait.label) +
  scale_y_log10(limits=c(1,1e+5)) +
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))

fdr <- fdr[order(fdr.overall)]
fdr[,trait.label:=factor(as.character(trait.label),levels=unique(fdr$trait.label))]
set.seed(42)
fdr[,ry:=runif(1:nrow(fdr))]
ggplot(fdr,aes(x=trait.label,fill=imd)) +
  geom_col(aes(y=1),col="grey",position="dodge",alpha=0.2) +
  geom_col(aes(y=fdr.overall),col="grey",position="dodge") +
  geom_label_repel(aes(y=ry,label=trait.label), data=fdr[imd==TRUE]) +
  scale_fill_viridis_d(option="A") +
  ## scale_x_discrete(breaks=fdr[imd==TRUE]$trait.label) +
  scale_y_continuous("FDR") + 
  theme(axis.text.x=element_blank())
  ## theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))

if(is.logical(fdr$imd))
    fdr[,imd:=ifelse(imd,"IMD","Other")]
fdr[,imd2:=ifelse(imd=="IMD","IMD",paste0(imd,"2"))]

theme_set(theme_cowplot(12))
p <- ggplot(fdr,aes(x=pmin(-log10(fdr.overall),6), y=N1,fill=imd)) +
  ## geom_rect(xmin=0.1,ymin=1100,xmax=0.01,ymax=mx,colour="red",fill="grey") +
  geom_point(pch=21,size=3,data=fdr[imd=="Other"],fill="DarkSlateBlue",alpha=0.5) +
  geom_point(pch=21,size=3,data=fdr[imd=="IMD"],fill="#fcfdbf") +
  geom_vline(xintercept=2) +
  geom_label_repel(aes(label=trait.label,fill=imd), data=fdr[sig==TRUE | imd=="IMD"],
                   force=2) +
  ## scale_fill_viridis_d("IMD",option="A",begin=0.9,direction=-1) +
  ## scale_colour_viridis_d("IMD",option="A",begin=0.9,direction=-1) +
  scale_colour_manual(values=c("IMD"="#fcfdbf","Other"="grey70")) +
  scale_fill_manual(values=c("IMD"="#fcfdbf","Other"="grey90")) +
  scale_y_log10("Number of cases",breaks=10^(c(2,3,4,5)),labels=c("100","1,000","10,000","100,000")) +
  scale_x_continuous("-log10 FDR",limits=c(0,6.5)) +
  background_grid()+
  guides(fill = guide_legend(title = "Disease",
                             override.aes = aes(label = ""))) +
  theme(legend.position="none")
## c(0.9,0.1),
##         legend.box.margin=margin(12, 12, 12, 12, "pt"),
##         ## legend.spacing=unit(12,"cm"),
##         legend.box.background=element_rect(size=0.5))
p

library( 'gridExtra' )

tt <- with(fdr, table("IMD"=imd,
                      "FDR<0.01"=ifelse(sig, "FDR<0.01","no")))
tt <- cbind(tt,"Total"=rowSums(tt))
tt <- cbind(tt,"(%)"=round(100 * tt[,1]/tt[,3],2))
tt <- tt[,c("Total","FDR<0.01","(%)")]
tt

tt3 <- ttheme_minimal(
  core=list(bg_params = list(fill = c("#fcfdbf", "grey90"), col=NA),
            fg_params=list(fontface=3)),
  colhead=list(bg_params = list(fill = "white", col=NA)))


x <- tableGrob(tt,theme=tt3)
ggdraw() +
  draw_plot( plot = p, x = 0, y = 0, width = 1, height = 1 ) +  # z
  draw_plot( plot = x, x = 0, y = 0.75, width = .4, height = .25 )


ggsave("figures/suppfig-ukbb-sig-by-imd.pdf",height=6,width=8,scale=1.6)

