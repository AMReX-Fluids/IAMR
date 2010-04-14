EPStoPDF = epstopdf

ALL: UsersGuide.ps UsersGuide.pdf 

eps_source = $(wildcard *.eps) \
             $(wildcard ConvertCheckpoint/*.eps) \
             $(wildcard Scaling/*.eps) \
             $(wildcard Verification/*.eps)

pdf_source = $(eps_source:.eps=.pdf)

section_sources = Introduction/Introduction.tex \
                  GettingStarted/GettingStarted.tex \
	          Inputs/Inputs.tex \
                  Equations/Equations.tex \
	          FlowChart/FlowChart.tex \
                  LevelSet/LevelSet.tex \
                  AMR/AMR.tex \
                  Visualization/Visualization.tex \
                  Software/Software.tex \
                  Suggestions/Suggestions.tex

UsersGuide.dvi: UsersGuide.tex $(eps_source) $(section_sources) $(eps_source)
	latex UsersGuide.tex < /dev/null
	bibtex UsersGuide
	latex UsersGuide.tex < /dev/null
	latex UsersGuide.tex < /dev/null
	latex UsersGuide.tex < /dev/null

UsersGuide.pdf: UsersGuide.tex $(pdf_source) $(section_sources) $(pdf_source)
	pdflatex UsersGuide.tex < /dev/null
	bibtex UsersGuide
	pdflatex UsersGuide.tex < /dev/null
	pdflatex UsersGuide.tex < /dev/null
	pdflatex UsersGuide.tex < /dev/null

pdf:	UsersGuide.pdf 

%.ps: %.dvi
	dvips -t letter -o $@ $<

%.pdf: %.eps
	$(EPStoPDF) $<

clean:
	$(RM) $(pdf_source) UsersGuide.dvi
	$(RM) UsersGuide.blg UsersGuide.log UsersGuide.toc
	$(RM) UsersGuide.aux UsersGuide.ps UsersGuide.bbl UsersGuide.pdf
	$(RM) *~

.PHONY: clean
