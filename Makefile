FPATH=figures # path to figures dir

## file to store sparse basis driver snps
SPARSEDRIVERS=/home/cew54/share/as_basis/sparse-basis/basis-sparse-0.999.RData
## change this, and everything changes, because new results
.DEPENDS=R/cw-reader.R R/cw-renamer.R R/cw-files.R R/cw-colours.R
# ## obsdata is the file of significant results
# OBSDATA=$(HOME)/basis13-delta-fdr0.01.csv

all: $(SPARSEDIVERS) figures/suppfig-consistency.pdf 

# $(OBSDATA) figures/figure2-hclust-shrinkage.png figures/figure3-big-cluster.pdf figures/suppfig-forest-everything.pdf 
## this is slow, but its a dependency, so don't wish to queue it
$(SPARSEDRIVERS) : R/make-sparse-basis-driver-snps-13.R
	Rscript $< > $(<)out 2>&1

## this is slow, and nothing depends on it
figures/suppfig-consistency.pdf : R/consistency.R $(SPARSEDRIVERS)
	qR.rb -r $< 

figures/figure2-hclust-shrinkage.png: R/figure2-cluster-with-without-shrinkage.R $(.DEPENDS)
	Rscript $< > $(<)out 2>&1



# $(OBSDATA) : extract-data-for-paulk.R $(.DEPENDS)
# 	Rscript $< > $(<)out 2>&1

# numbers.txt: numbers.R $(.DEPENDS)
# 	Rscript $< > numbers.txt 2>&1

# /home/cew54/sparse-projections.RDS: project.R
# 	qR.rb -r -c 2 $< > $(<)out 2>&1

# figures/figure3-big-cluster.pdf : figure3-cluster-everything-sparse.R $(.DEPENDS) $(OBSDATA)
# 	Rscript $< > $(<)out 2>&1

# figures/suppfig-forest-everything.pdf : forests-v3.R
# 	Rscript $< > $(<)out 2>&1

# figures/suppfig-ancestry.pdf : ancestry.R
# 	Rscript $< > $(<)out 2>&1

# figures/figures.pdf: figures/figures.tex figures/figure2-hclust-shrinkage.png
# 	cd figures && pdflatex figures && cd ..

# ## TODO: bb-icd comparison, code in comments in consistency.R
