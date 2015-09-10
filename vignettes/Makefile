all: meanShiftR.pdf 


meanShiftR.pdf:  meanShiftR.tex 
	latexmk -pdf -pdflatex="pdflatex -interaction=nonstopmode" -use-make meanShiftR.tex

clean: 
	latexmk -CA
