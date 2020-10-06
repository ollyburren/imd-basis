FPATH=figures # path to figures dir

## file to store sparse basis driver snps
SPARSEDRIVERS=/home/cew54/share/as_basis/sparse-basis/basis-sparse-0.999.RData
## change this, and everything changes, because new results
DEPENDS= $(wildcard R/cw-*.R) #R/cw-reader.R R/cw-renamer.R R/cw-files.R R/cw-colours.R
# ## obsdata is the file of significant results
# OBSDATA=$(HOME)/basis13-delta-fdr0.01.csv

#! all : main targets - figures.pdf and suppfigures.pdf
all: figures/figures.pdf figures/suppfigures.pdf

#! clean : remove every file in figures
clean:
	rm -r figures/*

.PHONY : help
help : Makefile
	@sed -n 's/^#!//p' $<

figures/figures.pdf : tex/figures.tex figures/figure2-hclust-rivas.png figures/figure3-big-cluster.pdf figures/fig4-pc1.pdf figures/fig4-pc3.pdf
	cp tex/figures.tex tex/overview_trimmed.pdf figures/
	cd figures && pdflatex figures && cd ..

figures/suppfigures.pdf : tex/suppfigures.tex figures/suppfig-forest-pc13.pdf figures/suppfig-ukbb-sig-by-imd.pdf figures/suppfig-consistency.pdf figures/suppfig-sparsesig-qqplots.pdf
	cp tex/suppfigures.tex figures/
	cd figures && pdflatex suppfigures && cd ..

################################################################################

## sub targets

## SPARSEDRIVERS : this is slow, but its a dependency, so don't wish to queue it
$(SPARSEDRIVERS) : R/make-sparse-basis-driver-snps-13.R $(DEPENDS)
	Rscript $< > $(<)out 2>&1

## these are slow, nothing else depends on them
figures/suppfig-consistency.pdf : R/consistency.R $(SPARSEDRIVERS)
	qR.rb -j consistency -r $< 

## these are quick enough to run interactively
figures/figure2-hclust-rivas.png: R/figure2+rivas.R $(DEPENDS)
	Rscript $< > $(<)out 2>&1

figures/figure3-big-cluster.pdf : R/figure3-cluster.R $(DEPENDS) 
	Rscript $< > $(<)out 2>&1

#! figures/fig4-pc1.pdf also makes suppfig forests
figures/suppfig-forest-pc13.pdf figures/fig4-pc1.pdf figures/fig4-pc3.pdf :: R/suppfig-forests.R
	Rscript $< > $(<)out 2>&1

figures/suppfig-ukbb-sig-by-imd.pdf : R/ukbb-by-imd.R $(DEPENDS)
	Rscript $< > $(<)out 2>&1

stats-sparse-drivers-significance.txt figures/suppfig-sparsesig-qqplots.pdf : R/table1.R
	Rscript $< > $(<)out 2>&1

figures/suppfig-proportionality.pdf : R/proportionality-v2.R

# $(OBSDATA) : extract-data-for-paulk.R $(DEPENDS)
# 	Rscript $< > $(<)out 2>&1

# numbers.txt: numbers.R $(DEPENDS)
# 	Rscript $< > numbers.txt 2>&1

# /home/cew54/sparse-projections.RDS: project.R
# 	qR.rb -r -c 2 $< > $(<)out 2>&1


# figures/suppfig-ancestry.pdf : ancestry.R
# 	Rscript $< > $(<)out 2>&1

# figures/figures.pdf: figures/figures.tex figures/figure2-hclust-shrinkage.png
# 	cd figures && pdflatex figures && cd ..

