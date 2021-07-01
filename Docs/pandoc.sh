#!/bin/sh

for i in `find ../UsersGuide/* -name "*.tex"`
do
    echo "\\input{textsymbols}" > _temp.tex
    cat $i >> _temp.tex
    pandoc _temp.tex --mathjax --wrap=preserve -o docs_converted/`basename $i .tex`.rst
    # Appending this to the command will make an in-page bib without links
    #--bibliography sphinx_documentation/source/refs.bib
    rm -f _temp.tex
done

