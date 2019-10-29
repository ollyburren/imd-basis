#!/usr/bin/env Rscript

## basis p values
## remotes::install_github('ollyburren/cupcake')
library(data.table)
library(magrittr)
library(ggplot2)
library(cupcake)
library(parallel)
library(pheatmap)
source("R/cw-files.R") # sets various file locations
source("R/cw-utils.R")

basis.DT<-get_gwas_data(TRAIT_MANIFEST_FILE,SNP_MANIFEST_FILE,GWAS_DATA_DIR,filter_snps_by_manifest=TRUE) 

shrink.DT<-readRDS(SHRINKAGE_FILE)
basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,"shrinkage")
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))

pc.emp <- readRDS(BASIS_FILE)

dt2mat <- function(dt,...) {
    tmp <- dcast(dt,...)
    rn <- tmp[[1]]
    m <- as.matrix(tmp[,-1])
    rownames(m) <- rn
    m
}

## check centering and rotating recreates pc.emp$x
zm <- scale(basis.mat.emp,center=TRUE,scale=FALSE)
zm.centre <- attr(zm,"scaled:center")
zz <- zm %*% pc.emp$rotation
cor(as.vector(zz),as.vector(pc.emp$x))
sum(abs(as.vector(zz)-as.vector(pc.emp$x)))

################################################################################

## how sparse are the rotations?
sX <- pc.emp$rotation 
m=reshape2::melt(sX)  %>% as.data.table()
m <- m[Var2!="PC14"]
head(m)

library(cowplot); theme_set(theme_cowplot())
ggplot(m,aes(x=value)) + geom_histogram(binwidth=0.01) + facet_wrap(~Var2) +
  scale_y_sqrt("Count (sqrt scale)",
               breaks=c(1,1e+3,1e+4,1e+5,#2e+5,
                        3e+5),labels=c("0","1,000","10,000","100,000",#"200,000",
                                       "300,000"),
               limits=c(0,3e+5)) + background_grid()
ggsave("~/share/as_basis/figures/suppfig-sparsebasis-coefficients.pdf",
       height=8,width=8)

## for diff quantiles of rotation vectors, what is correlation?
qthr <- seq(0.0001,0.005,by=0.0001)
cr0 <- cor(as.vector(zz),as.vector(pc.emp$x))
cr <- mclapply(qthr, function(q) {
    quants <- apply(-abs(sX),2,quantile,q)
    X <- pc.emp$rotation
    use <- matrix(TRUE,nrow(X),ncol(X),dimnames=dimnames(X))
    for(j in seq_along(quants)) {
        wh <- which(-abs(X[,j])>=quants[j])
        X[wh,j] <- 0
        use[wh,j] <- FALSE
    }
    message(q)
    user <- apply(use[,-14],1,any)
    beta0 = zm %*% X
    diag(cor(beta0,pc.emp$x))
})  %>% do.call("rbind",.)
head(cbind(qthr,cr))
tail(cbind(qthr,cr))

## what quantile gives minimum correlation?
CTHR=0.999 # minimum correlation
qlim <- qthr[apply(cr>CTHR,2,function(x) which(x)[1])][1:13]
quants <- sapply(1:13, function(i) quantile(-abs(sX[,i]), qlim[i]))
nquants <- sapply(1:13, function(i) sum(-abs(sX[,i]) <= quants[i]))
message(CTHR)
stats <- cbind(qlim,quants,nquants)
print(stats)

## 0.001 is good
X <- pc.emp$rotation[,1:13]
use <- matrix(TRUE,nrow(X),ncol(X),dimnames=dimnames(X))
for(j in seq_along(quants)) {
    wh <- which(-abs(X[,j])>quants[j])
    use[wh,j] <- FALSE
    X[wh,j] <- 0
}
beta0 = zm %*% X
diag(cor(beta0,pc.emp$x))

user <- apply(use[,-14],1,any)
table(user)
X <- pc.emp$rotation[user,-14]
Z <- zm[,which(user)]
use <- use[user,-14]
use.pca <- use
user.pca <- user
rot.pca <- X
rot.pca[!use.pca]=0

## par(mfrow=c(1,1))
hist(rowSums(use),main=paste("SNP use in sparse PCA",CTHR))

apply(use,2,mean) # 16-60% in each col
pids.use <- rownames(use)
sub(":.*","",pids.use)  %>% table()  ## all chromosomes represented
pids.bychr <- split(pids.use, sub(":.*","",pids.use))

man.DT <- readRDS(SNP_MANIFEST_FILE)

beta.centers <- structure(attr(zm,"scaled:center"),names=colnames(zm))[rownames(rot.pca)]

ishrink <- shrink.DT[match(rownames(rot.pca),pid),]
shrinkage <- structure(ishrink$shrinkage,names=ishrink$pid)

## LD
man.DT[,chr:=sub(":.*","",pid)]
man.DT <- man.DT[pid %in% pids.use]
SNP.manifest <- copy(man.DT)[match(pids.use,pid)]
s.DT <- split(man.DT, man.DT$chr)

LD <- lapply(names(pids.bychr), function(chr) {
        ss.file <- file.path(REF_GT_DIR, sprintf("%s.RDS", chr))
        sm <- readRDS(ss.file)
        pids <- colnames(sm)
        dup.idx <- which(duplicated(pids))
        if (length(dup.idx) > 0) {
            sm <- sm[, -dup.idx]
            pids <- pids[-dup.idx]
        }
        sm.map <- match(s.DT[[chr]]$pid, pids)
        if (any(is.na(sm.map))) {
            message("SNPs in manifest that don't have genotypes")
        }
        r <- ld(sm[, sm.map], sm[, sm.map], stats = "R")
        r[is.na(r)] <- 0
        r
})  %>% bdiag_with_dimnames(.)

LD <- LD[rownames(rot.pca),rownames(rot.pca)]

README <- c(rot.pca = "sparse rotation matrix (Q in manuscript)",
            use.pca = "hard thresholded rot.pca!=0",
            beta.centers = "vector of column means of 13 GWAS + control at sparse basis snps",
            shrinkage = "vector of w/sigma_maf at sparse basis snps",
            LD = "matrix of genotype correlations in 1000 Genomes EUR at sparse basis snps",
            SNP.manifest = "data.table of sparse basis snps + alleles, maf, ld.block",
            project.sparse="function to project new beta + seb into sparse basis")

project.sparse <- function(beta,seb,pids) {
    ##' assumes basis-sparse-13-0.999.RData has been loaded and that LD,
    ##' rot.pca, beta.centers, shrinkage are defined in current environment
    ##' beta = new beta, seb=se(beta), pids=snp ids in order of beta
    require(Matrix)
    if(length(beta)!=length(seb) || length(beta)!=length(pids) || !length(beta))
        stop("arguments must be equal length vectors > 0")
    if(!all(pids %in% SNP.manifest$pid))
        stop("all pids must be members of sparse basis (SNP.manifest$pid)")
    if(length(pids) < 0.95 * nrow(rot.pca))
        warning("more than 5% sparse basis snps missing")
    b <- beta * shrinkage[pids] - beta.centers[pids] # * shrinkage[pids]
    proj <- b %*% rot.pca[pids,]
    v <- seb * shrinkage[pids] * rot.pca[pids,]
    var.proj  <- t(v) %*% LD[pids,pids] %*% v  %>% diag()
    ctl <-  (-beta.centers[pids])  %*% rot.pca[pids,]
    ret <- data.table(PC=colnames(proj),
                      proj=proj[1,],
                      var.proj=var.proj,
                      delta=(proj - ctl)[1,])
    ret$z=ret$delta/sqrt(ret$var.proj)
    ret$p=pnorm(abs(ret$z),lower.tail=FALSE) * 2
    copy(ret)
}


save(rot.pca,use.pca,beta.centers,shrinkage,LD,SNP.manifest,
     project.sparse,README,stats,
     file="~/share/as_basis/sparse-basis/basis-sparse-0.999.RData",
     version=2)
