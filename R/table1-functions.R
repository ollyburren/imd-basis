
orp2bv <- function(or,p) {
    beta <- log(or)
    z <- qnorm(p/2,lower.tail=FALSE)
    se <- abs(beta/z)
    list(beta=beta,vbeta=se^2)
}


lapply_with_names <- function(x,...) {
    ret <- lapply(x, ...)
    names(ret) <- x
    ret
}


lapply_copy_names <- function(x,...) {
    ret <- lapply(x, ...)
    names(ret) <- names(x)
    ret
}

geta <- function(target) {
    mx <- target$pos + 1e+6
    mn <- target$pos - 1e+6
    chr <- target$chr
    cat(chr,":",mn,"-",mx,"\n")
    comm=paste0("zcat ~/share/Data/reference/1000GP_Phase3/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.0.01MAF.txt.gz | awk '$1==",target$chr," && $2 > ",mn," && $2 < ",mx,"'")
    kg <- fread(cmd=comm,sep="\t")[,1:5]
    setnames(kg,c("chr","pos","rsid","ref","alt"))
    kg[,pid:=paste(chr,pos,sep=":")]
    copy(kg)
}

limit_chr_pos <- function(chr,position,w=1e+6) {
    keep <- rep(FALSE,length(chr))
    position <- as.numeric(position)
    for(target in targets)
        keep <- keep | ( chr==target$chr & position > target$pos- w & position < target$pos+ w)
    keep
}

library(annotSnpStats)

## fix alleles
switch_indicator <- function(pid,a1,a2,verbose=FALSE) {
    m <- match(pid,alleles$pid)
    tmp <- data.table(pid=pid,a1=a1,a2=a2,
                      ref=alleles$ref[m],alt=alleles$alt[m],
                      order=seq_along(pid))
    cl <- with(tmp, g.class(paste(a1,a2,sep="/"),
                            paste(ref,alt,sep="/")))
    print(table(cl))
    if(any(cl %in% c("rev","revcomp"))) {
        if(verbose)
            print(tmp[cl %in% c("rev","revcomp")])
              classes <- c(ambig=NA,
                     impossible=NA,
                     nochange=FALSE,
                     rev=TRUE,
                     comp=FALSE,
                     revcomp=FALSE)
        ret <- classes[cl]
        if(any(cl=="rev") && any(cl=="revcomp")) { # can't tell what to do with ambig
            return(ret)
        } else if (any(cl=="rev")) { # no comp
            ident <- with(tmp, paste(a1,a2,sep="/") == paste(ref,alt,sep="/"))
            rev <- with(tmp, paste(a1,a2,sep="/") == paste(alt,ref,sep="/"))
            ret[ident] <- FALSE
            ret[rev] <- TRUE
            return(ret)
        }
        stop("revcomp not implemented")
    } else {
        return(rep(TRUE,length(pid)))
    }
}

## coloc
make.data <- function(target,t) {
    message(t)
    x=target$data[[t]]
    x <- x[!(chr==5 & position < 110500000)]
    d1 <- list(beta=x$beta,
               varbeta=x$vbeta,
               snp=x$pid,
               type="cc",
               N=with(N[[t]], sum(N1+N0)),
               s=with(N[[t]], sum(N1)/sum(N0+N1)))
}

library(coloc)
fcoloc <- function(target) {
    ret <- lapply_with_names(target$test, function(t1) {
        d1 <- make.data(target,t1)
        lapply_with_names(target$cond, function(t2) {
            d2 <- make.data(target,t2)
            tmp <- coloc.abf(d1,d2,p12=5e-6)
            wh <-  which.max(tmp$results$SNP.PP.H4)
            best <- tmp$results$snp[wh]
            best.pp <- tmp$results$SNP.PP.H4[wh]
            pp4 <- tmp$summary[["PP.H4.abf"]]
            list(PP4=pp4, best.snp=best,best.pp=best.pp)
        })
    })
    target$coloc <- ret
    target
}

## functions to add LD
library(simGWAS)
gethap <- function(cmd) {
    y=fread(cmd=cmd)
    ha <- simGWAS:::vcf2haps(as.matrix(y[,-c(1:9)]))
    rownames(ha) <- paste(gsub("chr","",y$V1),y$V2,sep=':')
    t(ha)
}
cor2 <- function (x) {
    1/(NROW(x) - 1) * crossprod( scale(x, TRUE, TRUE) )
}
BCF_TOOLS='/usr/local/Cluster-Apps/bcftools/1.2/bin/bcftools'
BCF_FILE='~/share/Data/reference/UK10K/chr16.bcf.gz'
bcf_maf=0.01
getld <- function(target) {
    mx <- target$pos + 1e+6
    mn <- target$pos - 1e+6
    chr <- target$chr
    allpid <- lapply(target$data, "[[", "pid")  %>% unlist()  %>% unique()
    cat(chr,":",mn,"-",mx,"\n")
    ld.cmd <- sprintf("%s view -H /home/ob219/rds/rds-cew54-wallace-share/Data/reference/UK10K/BCF/chr%s.bcf.gz --min-af %f:minor --max-alleles 2 --min-alleles 2 -r chr%s:%s-%s -Ov",
                  BCF_TOOLS,chr,bcf_maf,chr,mn,mx)
    h <- gethap(ld.cmd)
    use <- apply(h,2,var)>0 & colnames(h) %in% allpid
    h <- h[,use,drop=FALSE]
    list(ld=cor2(h),maf=colMeans(h))
}


maxN <- function(x, N=2){
    len <- length(x)
    if(any(N==0))
        N[N==0] <- 1
    if(any(N>len))
    ##     warning('N greater than length(x).  Setting N=length(x)')
        N <- pmin(N,length(x))
    ## }
    sort(x,partial=len-N+1)[len-N+1]
}

cbind_minlength <- function(x,y) {
    n <- min(length(x),length(y))
    cbind(x[1:n],y[1:n])
}

splitter <- function(m,thr=0.1) {
    if(any(m<0)) # r, not rsq
        m <- m^2
    n <- ncol(m)
    ## mx <- sapply(2:n, function(i) {
    ##     j <- i-1
    ##     maxN(m[1:j,i:n],round(((n-i+1) * j) * c(0.001,0.005,0.01)))
    ## })  %>% t()
    mx <- sapply(2:(n-1), function(i) {
        j <- (i-1):1
        k <- (i+1):n
        vec <- m[cbind_minlength(j, k)]
        max(vec)
    })  %>% t()

    ## iterate.
    ## 1. find min, record
    ## 2. set ,mx=1 at all snps in neighbourhood
    ## 3. stop if remaining min > the
    splits <- numeric(100)
    i <- 0
    mx0 <- mx
    while(min(mx0) < thr) {
        i <- i+1
        splits[i] <- which.min(mx0)
        idx <- (splits[i]-100):(splits[i]+100)
        idx <- idx[ idx >=1 & idx<=(n-2) ]
        mx0[ idx ] <- 1
        if(i>100)
            break
    }
    splits <- splits[1:i]
    
    ## plot(2:(n-1),mx,type="l")
    ## abline(h=0.1)
    ## abline(v=splits,col="red",lty=2)
    
    ## corrplot(m,method="color")
    ## abline(h=n-(splits + 1))
    ## abline(v=splits + 1)
  
    splitat <- structure(splits + 1, names=colnames(m)[ splits + 1 ])
    splitpids=colnames(m)[ splits + 1 ]
    grp=cumsum( colnames(m) %in% splitpids)
    grppids <- split(colnames(m), grp)
}
    
plotmann <- function(nm,pp4=TRUE,addsnps=TRUE,ti="") {
    target <- targets[[nm]]
    if(ti=="")
        ti <- nm
    z <- rbindlist(target$data,use.names=TRUE,fill=TRUE)
    z[,n:=.N,by="pid"]
    z <- z[n==length(target$data) | pid %in% c("13:43056036","5:110401872")]
    p <- ggplot(z,aes(x=position,y=-log10(p.value),colour=trait)) +
      geom_vline(xintercept = target$pos) +
      geom_point() +
      facet_grid(trait ~ ., scales="free") +
      labs(x=paste0("Chr ",target$chr)) +
      ggtitle(ti) +
      scale_color_brewer(palette="Paired") +
      theme(plot.title=element_text(hjust=0.5),
            legend.position="none",
            strip.text.y = element_text(angle = 0))
    if(addsnps) {
    if("13:43056036" %in% z$pid)
        p <- p + geom_vline(xintercept=43056036,col=col1,lty="dotted") +
          geom_text(x=43056036,y=4,label="rs34132030",hjust="outward",col=col1,angle=-90)
    if("5:110401872" %in% z$pid)
        p <- p + geom_vline(xintercept=110401872,col=col1,lty="dotted") +
          geom_text(x=110401872,y=4,label="rs1837253",hjust="outward",col=col1,angle=-90)
    }

    cdata <- lapply(target$coloc, function(x) {
        y <- do.call("rbind",lapply(x, as.data.frame))
        y$cond <- rownames(y)
        y})
    for(nm in names(cdata))
        cdata[[nm]]$trait <- nm
    cdata %<>% do.call("rbind",.)
    cdata$position <- min(z$position)
    rownames(cdata) <- NULL

    ctab <- with(cdata, matrix(round(PP4,2),length(unique(cond)),length(unique(trait)),
                               dimnames=list(basis=unique(cond),unique(trait))))  %>%
      t()  %>%
      as.table()
    print(ctab)

    if(pp4) {
        title <- textGrob("PP4",gp=gpar(fontface=2))
        t1 <- tableGrob(ctab,theme=tt2)
        padding <- unit(5,"mm")
        
        table <- gtable_add_rows(
     t1, 
     heights = grobHeight(title) + padding,
     pos = 0)
table <- gtable_add_grob(
    table, 
    title, 
    1, 1, 1, ncol(table))

    ggdraw() +
      draw_plot(p) +
      draw_plot( plot = table,
                   x = 0.1, y = 0.8, width = 0.3, height = .2 ) +
      background_grid()
    } else {
          return(p)
        }
}

## coloc with LD
make.data_ld <- function(target,t,usepids) {
    message(t)
    x=copy(target$data[[t]])
    ## x <- x[!(chr==5 & position < 110500000)]
    x <- x[pid %in% usepids]
    ## x[sw==TRUE,beta:=-beta]
    d1 <- list(beta=x$beta,
               varbeta=x$vbeta,
               snp=x$pid,
               type="cc",
               ## MAF=ld$maf[x$pid],
               ## LD=ld$ld[x$pid,x$pid],
               N=with(N[[t]], sum(N1+N0)),
               s=with(N[[t]], sum(N1)/sum(N0+N1)))
}


devtools::load_all("~/RP/coloc")
fcoloc_ld <- function(nm) {
    target <- targets[[nm]]
    ## ld <- LD[[nm]]
    targetpid <- paste(target$chr,target$pos,sep=":")
    use <- sapply(splits[[nm]], function(x) targetpid %in% x)  %>% which()
    usepids <- splits[[nm]][[use]]
    ret <- lapply_with_names(target$test, function(t1) {
        d1 <- make.data_ld(target,t1,usepids)
        lapply_with_names(target$cond, function(t2) {
            d2 <- make.data_ld(target,t2,usepids)
            tmp <- coloc.abf(d1,d2,p12=5e-6)
            ## tmp <- coloc.signals(d1,d2,p12=5e-5,method="mask")
            ## cbind(data.table(test=t1,cond=t2),
            ##       tmp$summary[,.(hit.test=hit1,hit.cond=hit2,PP.H0.abf,PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf)] )
            wh <-  which.max(tmp$results$SNP.PP.H4)
            best <- tmp$results$snp[wh]
            best.pp <- tmp$results$SNP.PP.H4[wh]
            pp4 <- tmp$summary[["PP.H4.abf"]]
            list(PP4=pp4, best.snp=best,best.pp=best.pp)
        })
    })  #%>% rbindlist()
}

plotmann_ld <- function(nm,pp4=TRUE,addsnps=TRUE,ti="") {
    target <- targets[[nm]]
    if(ti=="")
        ti <- nm
    z <- rbindlist(target$data,use.names=TRUE,fill=TRUE)
    z[,n:=.N,by="pid"]
    z <- z[n==length(target$data) | pid %in% c("13:43056036","5:110401872")]
    targetpid <- paste(target$chr,target$pos,sep=":")
    use <- sapply(splits[[nm]], function(x) targetpid %in% x)  %>% which()
    usepids <- splits[[nm]][[use]]
    zuse <- z[pid %in% usepids,]
    p <- ggplot(z,aes(x=position,y=-log10(p.value),colour=trait)) +
      geom_vline(xintercept=c(min(zuse$position),max(zuse$position)),col="grey",lty=2) +
      geom_vline(xintercept = target$pos) +
      geom_point() +
      facet_grid(trait ~ ., scales="free") +
      labs(x=paste0("Chr ",target$chr)) +
      ggtitle(ti) +
      scale_color_brewer(palette="Paired") +
      theme(plot.title=element_text(hjust=0.5),
            legend.position="none",
            strip.text.y = element_text(angle = 0))
    if(addsnps) {
    if("13:43056036" %in% z$pid)
        p <- p + geom_vline(xintercept=43056036,col=col1,lty="dotted") +
          geom_text(x=43056036,y=4,label="rs34132030",hjust="outward",col=col1,angle=-90)
    if("5:110401872" %in% z$pid)
        p <- p + geom_vline(xintercept=110401872,col=col1,lty="dotted") +
          geom_text(x=110401872,y=4,label="rs1837253",hjust="outward",col=col1,angle=-90)
    }

    cdata <- lapply(target$scoloc, function(x) {
        y <- do.call("rbind",lapply(x, as.data.frame))
        y$cond <- rownames(y)
        y})
    for(nm in names(cdata))
        cdata[[nm]]$trait <- nm
    cdata %<>% do.call("rbind",.)
    cdata$position <- min(z$position)
    rownames(cdata) <- NULL

    ctab <- with(cdata, matrix(round(PP4,2),length(unique(cond)),length(unique(trait)),
                               dimnames=list(basis=unique(cond),unique(trait))))  %>%
      t()  %>%
      as.table()
    print(ctab)

    if(pp4) {
        title <- textGrob("PP4",gp=gpar(fontface=2))
        t1 <- tableGrob(ctab,theme=tt2)
        padding <- unit(5,"mm")
        
        table <- gtable_add_rows(
     t1, 
     heights = grobHeight(title) + padding,
     pos = 0)
table <- gtable_add_grob(
    table, 
    title, 
    1, 1, 1, ncol(table))

    ggdraw() +
      draw_plot(p) +
      draw_plot( plot = table,
                   x = 0.1, y = 0.8, width = 0.3, height = .2 ) +
      background_grid()
    } else {
          return(p)
        }
}
