#!/bin/sh

rst2latex.py --stylesheet test.sty AlphaMate.rst AlphaMate.tex

pdflatex AlphaMate.tex
bibtex AlphaMate.tex
makeindex AlphaMate.tex
pdflatex AlphaMate.tex
pdflatex AlphaMate.tex

open AlphaMate.pdf
