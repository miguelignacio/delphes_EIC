#!/bin/bash

#jupyter-nbconvert --to pdf --template hidecode.tplx PHYS\ 1303\ -\ Spring\ 2018\ -\ Problems\ and\ Solutions.ipynb
# Get the root filename
fileroot=$(echo $1 | sed -e 's/\(.*\)\.ipynb/\1/')
jupyter-nbconvert --to latex --template hidecode.tplx "$1"
pdflatex ${fileroot}.tex

