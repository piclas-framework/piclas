# Documentation

## Building documentation

The user and developer guides are built automatically via Read The Docs.
When changing the guides, build the html and pdf files locally before committing changes to the repository.

**Prerequisites**

The following shows how to create the html and pdf files for the user/developer guide.
First, install the required prerequisites. Install *python3* and make sure that pip (third-party Python packages) is installed.
If you have an old version of, e.g., Ubuntu visit [this website](https://phoenixnap.com/kb/how-to-install-python-3-ubuntu).

    sudo apt install python3-pip

Navigate to the documentation folder from the PICLas top level directory

    cd docs/documentation

Run pip to install the required extensions and packages for compiling the user guide (only once)

    python3 -m pip install --exists-action=w --no-cache-dir -r requirements.txt

Make sure that *latexmk* is installed on the system for compiling the PDF version of the user guide. For Ubuntu, follow
[this link](https://zoomadmin.com/HowToInstall/UbuntuPackage/latexmk) for installation.

    sudo apt-get install latexmk


**HTML Version**

Compile the html version of the user guide via

    python3 -m sphinx -T -E -b html -d _build/doctrees -D language=en . _build/html

Check that no errors occur during compilation and then navigate to the created html files

    cd _build/html

Open index.html to see if everything has worked out correctly (e.g. with your favourite browser).
Note that you can simply run the script *buildHTML.sh* in the *documentation* directory for this task.


**PDF Version**

Next, create the pdf output.

    python3 -m sphinx -b latex -D language=en -d _build/doctrees . _build/latex

and switch into the output directory

    cd _build/latex

Finally, compile the pdf file

    latexmk -r latexmkrc -pdf -f -dvi- -ps- -jobname=piclas -interaction=nonstopmode

and check if the pdf exists

    ls _build/latex/piclas.pdf

Note that you can simply run the script *buildPDF.sh* in the *documentation* directory for this task.

## Writing documentation

### Figures

Graphics are supported as long as their size is not above 1 MB, recommended size is below 100 KB. The conversion of PDF plots/graphics produced with pgfplots/tikz can be performed via terminal using the `libvips` package available for most distributions

    sudo apt install libvips-dev
    vips copy example.pdf[dpi=150] example.jpg[Q=90,strip]

Modify `dpi=150` to scale the PDF and `Q=90` (between 0 and 100) to change the quality of the JPEG. Vector graphics might give better quality using SVG as the end format by converting from PDF

    pdf2svg example.pdf example.svg