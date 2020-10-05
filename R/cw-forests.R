library(cowplot)
source("R/cw-palette.R")
cols <- hrbrthemes::ipsum_pal()(7); c3 <- cols[3]
order.traits <- function(dat) {
    ## ukbb always left to right
    sdat <- split(dat,dat$newcat)
    strait <- lapply(sdat, function(idat) {
        with(idat,
             trait.label[order(delta,decreasing=TRUE)])  %>% unique()
    })
    traits <- unlist(strait)
    combi <- grep("_combined",traits)
    if(length(combi))
        traits <- c(" _combined",traits[-combi])
    traits
}

forest_labelled <- function(proj.dat,basis.dat=basis.DT,pc,focal=NULL,title="",fdr_thresh=0.01,samples=FALSE,order.within=TRUE){
    dat <- proj.dat[PC==pc,]
    dat[,trait.label:=as.character(trait.label)]
    dat[,category.label:=as.character(category.label)]
    ## remove non-sig UKBB
    ## dat[category.label %in% c("bb_disease","bb_cancer"),category.label:="UKBB"]
    dat <- dat[category.label!="UKBB" | newfdr < fdr_thresh | trait.label %in% basis.dat$trait.label]
    ## add CI
    dat[,ci:=1.96 * sqrt(variance)]
    dat[,c('lower','upper'):=list(delta-ci,delta+ci)]
    ## add basis
    dat <- rbind(basis.dat[PC==pc & trait.label!='control',],dat,fill=TRUE)

    ## order categories left to right, UKBB at bottom
    dat[,cat.delta:=mean(delta),by="category.label"]
    ocat <- c("basis","UKBB")
    ncat=unique(dat[order(cat.delta)]$category.label)  %>% sort()
    ## if cytokines or blood counts, put last before UKBB
    if(length(w <- which(ncat %in% c("Cytokines","Blood counts"))))
        ncat <- c(ncat[w],ncat[-w])
    ocat=c(setdiff(ncat,ocat),ocat)
    dat[,category.label:=factor(category.label,levels=ocat)]
    dat[,newcat:=category.label]
    dat[category.label %in% c("Basis","basis"),newcat:="UKBB"]
    dat[,full.label:=paste(category.label,trait.label,sep="/")]

    ## order traits left to right within UKBB, and alphabetically otherwise
    otraits <- order.traits(dat)
    dat[,full.label:=factor(trait.label,levels=otraits)]
    
    forests <- ggplot(dat,aes(x=full.label,y=delta,#colour=category.label,
                              linetype=newfdr>fdr_thresh)) +
      geom_linerange(aes(ymin=lower,ymax=upper),col=c3) +
      geom_point(aes(pch=grepl("_combined",trait.label)),
                 size=4,col="grey80",
                 data=dat[!(category.label %in% c("basis","Basis"))]) + #,
                                        #bg="black") +
      geom_point(data=dat[category.label %in% c("basis","Basis")],
                 col="red",pch=15) +
      geom_hline(yintercept=0,col='red',linetype=3) +
      ggtitle(paste("Component",sub("PC","",pc))) +
      facet_grid(newcat ~ .,scales="free_y",space="free_y",switch="y") +
      scale_y_continuous("Difference from control",breaks=0) +
      ## scale_x_discrete(breaks=seq_along(otraits), labels=olabels) +
      #scale_shape_manual(values=c("TRUE"=23,"FALSE"=22),guide=FALSE) +
      scale_shape_manual(values=c("TRUE"=18,"FALSE"=15),guide=FALSE) +
      scale_linetype(guide=FALSE) +
      ## scale_colour_gradient2(low=darken(myblue),high=darken(myred),mid="grey95",guide="none") +
      coord_flip() +
      background_grid(major="y") +
      theme(plot.title=element_text(hjust = 0),
            strip.placement = "outside",
            ## strip.background=element_rect(fill="grey90"),
strip.background=element_blank(),
strip.text = element_text(hjust = 1,vjust=1,face="bold"),
panel.spacing = grid::unit(2, "lines"),
strip.text.y.left = element_text(angle=0),
## axis.ticks.x=element_blank(),
            ## axis.text.x=element_blank(),
            ## legend.position="none",
            axis.line=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title=element_blank(),
            strip.text.y = element_text(angle=180))
    
if(samples) {
    p.samples <- ggplot(dat[],aes(x=trait.label,ymin=0,ymax=n1,y=n1)) +
        geom_pointrange() +
         coord_flip() + 
        facet_grid(newcat ~ .,scales="free_y",space="free_y",switch="y") +
      background_grid() +
      theme(plot.title=element_text(hjust = 0),
            strip.background=element_blank(),
            axis.line.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title.y=element_blank(),
            strip.text.y = element_blank(),
            legend.position="none") +
    scale_y_log10(breaks=10^c(2,3,4),labels=c("100","1000","10000")) +
      xlab("Trait") + ylab("Number of cases (log scale)")

    plot_grid(forests,p.samples,rel_widths=c(0.8,0.2),align="h",axis="b")
} else {
    forests
    }
    
}
