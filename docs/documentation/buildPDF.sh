# Compile latex files
python3 -m sphinx -b latex -D language=en -d _build/doctrees . _build/latex

# Switch to latex source files
cd _build/latex

# Compile pdf file(s)
latexmk -r latexmkrc -pdf -f -dvi- -ps- -jobname=piclas -interaction=nonstopmode
