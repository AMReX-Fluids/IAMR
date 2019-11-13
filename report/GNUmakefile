EPStoPDF = epstopdf

ALL: report.pdf 

eps_source = $(wildcard *.eps) $(wildcard */*.eps)

pdf_source = $(eps_source:.eps=.pdf)

section_sources =

report.dvi: report.tex $(eps_source) $(section_sources) $(eps_source)
	latex report.tex < /dev/null
#	bibtex report
#	latex report.tex < /dev/null
#	latex report.tex < /dev/null

report.pdf: report.tex $(pdf_source) $(section_sources) $(pdf_source)
	pdflatex report.tex < /dev/null
#	bibtex report
#	pdflatex report.tex < /dev/null
#	pdflatex report.tex < /dev/null

pdf:	report.pdf 

%.ps: %.dvi
	dvips -t letter -o $@ $<

%.pdf: %.eps
	$(EPStoPDF) $<

clean:
	$(RM) $(pdf_source) report.dvi
	$(RM) report.blg report.log report.toc
	$(RM) report.aux report.ps report.bbl report.pdf
	$(RM) report.lot report.lof report.out
	$(RM) *~

.PHONY: clean
