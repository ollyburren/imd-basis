## code to perform LD pruning as implemented in DeGAS using plink
ROOT.DIR <- '~/share/as_basis/basis-creation'
SNP_MANIFEST_FILE <- file.path(ROOT.DIR,'support/13_snp_manifest.RDS')
BCFOUTDIR <- '/home/ob219/share/as_basis/GWAS/1KG_Nov_19'
man.DT <-readRDS(SNP_MANIFEST_FILE)[,c('chr','pos'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)][order(chr,pos),]
OUTDIR <- '/home/ob219/share/as_basis/GWAS/1KG_Nov_19'
options('scipen'=999)
lapply(split(man.DT,man.DT$chr),function(dt){
  chrname<-dt[1,]$chr
  fname <- file.path(OUTDIR,sprintf("%s.txt",chrname))
  dt <- dt[order(pos),]
  write.table(dt[,.(chr,pos)],file=fname,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
})
options('scipen'=0)

## prepare sample file put in BCFOUTDIR
## wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.xlsx
library(xlsx)
samples <- read.xlsx(file.path(BCFOUTDIR,'20130606_sample_info.xlsx'),1) %>% data.table
EUR <- samples[Population %in% c('CEU','FIN','GBR','IBS','TSI'),.(Sample,Population)]
write(EUR$Sample,file.path(BCFOUTDIR,'eur_samples.txt'))


UK1KDIR <- '/home/ob219/rds/hpc-work/reference_panels/1000genomes'
BCFOUTDIR <- '/home/ob219/share/as_basis/GWAS/1KG_Nov_19'
for(chrname in 1:22){
   message(chrname)
   ukf <- file.path(UK1KDIR,sprintf("ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",chrname))
   rfile <- file.path(OUTDIR,sprintf("%s.txt",chrname))
   ofile <- file.path(BCFOUTDIR,sprintf("%s.bcf",chrname))
   cmd <- sprintf("/usr/local/software/spack/spack-0.11.2/opt/spack/linux-rhel7-x86_64/gcc-5.4.0/bcftools-1.6-6p3lqyearbsrree33gn7blgmuxjiybc5/bin/bcftools view --force-samples -R %s -S %s -Ob %s > %s",rfile,file.path(BCFOUTDIR,'eur_samples.txt'),ukf,ofile)
   write(cmd,file="~/tmp/qstuff/filter.txt",append=TRUE)
   #system(cmd)
   ## run these outside of R so as to preserve session from watchdog
}

## do in shell
if(FALSE){
  for i in `\ls *.bcf`; do
  echo "Processing $i"
  cp $i $i.back
  bcftools sort $i | bcftools view -Ob -e 'ALT[*]~"CN"' > $i.sorted
  mv $i.sorted $i
  done

  ## this creates a list of all SNPs in 1K genomes that are also in SNP manifest.
  for i in `\ls *.bcf`; do
  echo "Processing $i"
  bcftools view -H $i -Ov | cut -f1-5 >> all.snps.txt
  done

  cat *.prune.in > all.prune.in

  ## if thing go wrong this gets the backup
  for i in `\ls *.back`;do
  outfile=$(basename $i .back)
  echo "mv $i $outfile"
  mv $i $outfile
  done

  for i in `\ls *.bcf`; do
  prefix=$(basename $i .bcf)
  plink --bcf $i --indep 50 5 2 --out $prefix --allow-extra-chr
  done

}

gen.snps <- fread(file.path(BCFOUTDIR,'all.snps.txt'))
sel.snps <- scan(file.path(BCFOUTDIR,'all.prune.in'),'character()')

gen.snps <- gen.snps[V3 %in% sel.snps,][,pid:=paste(V1,V2,sep=':')]

pruned.snp.manifest <- man.DT[pid %in% gen.snps$pid,]
saveRDS(pruned.snp.manifest,file=file.path(ROOT.DIR,'support/13_trait_snp_manifest-degas-ld-pruned.RDS'))
