RMD := $(wildcard *.rmd)
HTML := $(RMD:.rmd=.html)
PDF := $(RMD:.rmd=.pdf)
DOCX := $(RMD:.rmd=.docx)

COMMON := nature.csl gtex.bib

all: $(HTML)

docx: $(DOCX)

pdf: $(PDF)

%.html: %.rmd $(COMMON)
	Rscript -e "library(rmarkdown); render('$<');"

%.pdf: %.rmd $(COMMON)
	Rscript -e "library(rmarkdown); render('$<', 'pdf_document');"

%.docx: %.rmd $(COMMON)
	Rscript -e "library(rmarkdown); render('$<', 'word_document');"

