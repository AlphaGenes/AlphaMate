#!/bin/sh

rst2latex.py --stylesheet preamble.sty AlphaMate.txt AlphaMate.tex

pdflatex AlphaMate.tex
bibtex AlphaMate.tex
makeindex AlphaMate.tex
pdflatex AlphaMate.tex
pdflatex AlphaMate.tex

open AlphaMate.pdf
