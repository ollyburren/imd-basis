##' creates list of genes from ensembl, and a function `annot` to
##' annotate any vector of pid by nearby genes within a window, w

library(GenomicRanges)
if(!file.exists("~/genes.RData")) {
    library(randomFunctions)
    library(biomaRt)
    mart <- useMart('ENSEMBL_MART_ENSEMBL',host="grch37.ensembl.org")
    ensembl.gene <- useDataset("hsapiens_gene_ensembl",mart=mart)
    genes<-getBM(
        ## filters= c("chromosome_name","start","end"),
        attributes= c('ensembl_gene_id','chromosome_name','start_position','end_position','external_gene_name','strand','gene_biotype','entrezgene'),
        ## values= ss,
        mart= ensembl.gene)
    genes <- subset(genes, !grepl("HG|HS|GL",chromosome_name))
    save(genes,file="~/genes.RData")
} else {
    load("~/genes.RData")
}

library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
load("~/genes.RData")
genes <- subset(genes, !grepl("HG|HS|GL",chromosome_name))
genes <- subset(genes, gene_biotype=="protein_coding")
ggr <- GRanges(Rle(genes$chromosome_name),
               IRanges(start=genes$start_position,
                       end=genes$end_position))
mcols(ggr) <- genes[,c("ensembl_gene_id","external_gene_name","entrezgene")]

annot <- function(pid,snp=TRUE,w=5e+5) {
    ss <- strsplit(pid,":") 
    chr <- sapply(ss,"[[",1)
    bp <- sapply(ss,"[[",2)  %>% as.numeric()
    left <- bp-w
    right <- bp + w
    gr0 <- GRanges(seqnames=Rle(chr,rep(1,length(chr))),
                   IRanges(start=bp,end=bp))
    mcols(gr0) <- DataFrame(pid=pid)
    ## add snps
    if(snp) {
        df <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP144.GRCh37, gr0)
        df$pid <- paste(seqnames(df),start(df),sep=":")
        m <- match(gr0$pid,df$pid)
        mcols(gr0) <- cbind(mcols(gr0),mcols(df)[m,])
    }
    ## add genes
    start(gr0) <- left
    end(gr0) <- right
    ## ggr <- subsetByOverlaps(ggr,gr0)
    hits <- findOverlaps(ggr, gr0)
    ret <- ggr[ queryHits(hits) ]
    mcols(ret)  %<>%  cbind(., mcols(gr0[ subjectHits(hits) ]))
    ret %<>% as.data.frame()
    ret
}

