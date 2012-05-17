# Makefile to use knitr for package vignettes

# put all PDF targets here, separated by spaces
PDFS= Lab2.pdf Lab3.pdf Lab4.pdf Lab5.pdf Lab6.pdf Lab7.pdf Lab8.pdf

all: $(PDFS) 

clean:
	rm -f *.tex *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary

%.pdf: %.Rnw
	$(R_HOME)/bin/Rscript -e "library(knitr); knit2pdf('$*.Rnw')"
