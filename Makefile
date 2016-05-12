all:

	Rscript -e "require(knitr); knit('main.rnw'); purl('main.rnw')" && texi2pdf main.tex
