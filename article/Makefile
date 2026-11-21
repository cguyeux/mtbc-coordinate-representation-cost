main.pdf: main.tex references.bib
	pdflatex -interaction=nonstopmode main.tex
	bibtex main
	pdflatex -interaction=nonstopmode main.tex
	pdflatex -interaction=nonstopmode main.tex

clean:
	rm -f main.aux main.bbl main.blg main.log main.out

.PHONY: clean

supplementary.pdf: supplementary.tex references.bib supplementary_generated/tableS1.tex supplementary_generated/tableS2.tex
	pdflatex -interaction=nonstopmode supplementary.tex
	bibtex supplementary
	pdflatex -interaction=nonstopmode supplementary.tex
	pdflatex -interaction=nonstopmode supplementary.tex
