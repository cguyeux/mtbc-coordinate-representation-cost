figures/fig_gyra94.pdf: ../analyses/p717_figure_gyra94.py
	cd .. && python3 analyses/p717_figure_gyra94.py
	cd figures && pdflatex -interaction=nonstopmode fig_gyra94.tex

main.pdf: main.tex references.bib figures/fig_gyra94.pdf
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
