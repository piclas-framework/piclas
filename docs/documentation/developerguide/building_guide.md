# Building the Documentation

The user and developer guides are built automatically via Read The Docs.
When changing the guides, build the html and pdf files locally before committing changes to the repository.

**Prerequisites**

First, install the required prerequisites. Install *python3* and make sure that pip (third-party Python packages) is installed.
If you have an old version of, e.g., Ubuntu visit [this website](https://phoenixnap.com/kb/how-to-install-python-3-ubuntu).

    sudo apt install python3-pip

The following examples shows how to create the html and pdf files for the user guide.
Navigate to the user guide directory from the PICLas top level directory

    cd docs/userguide

Run pip to install the required extensions and packages for compiling the user guide (only once)

    python3 -m pip install --exists-action=w --no-cache-dir -r ../requirements.txt

Make sure that *latexmk* is installed on the system for compiling the PDF version of the user guide. For Ubuntu, follow
[this link](https://zoomadmin.com/HowToInstall/UbuntuPackage/latexmk) for installation.

    sudo apt-get install latexmk


**HTML Version**

Compile the html version of the user guide via

    python3 -m sphinx -T -E -b html -d _build/doctrees -D language=en . _build/html

Check that no errors occur during compilation and then navigate to the created html files

    cd _build/html

Open index.html to see if everything has worked out correctly (e.g. with your favourite browser).
Note that you can simply run the script *buildHTML.sh* in the *userguide* directory for this task.


**PDF Version**

Next, create the pdf output.

    python3 -m sphinx -b latex -D language=en -d _build/doctrees . _build/latex

and switch into the output directory

    cd _build/latex

Finally, compile the pdf file

    latexmk -r latexmkrc -pdf -f -dvi- -ps- -jobname=piclas -interaction=nonstopmode

and check if the pdf exists

    ls _build/latex/piclas.pdf

Note that you can simply run the script *buildPDF.sh* in the *userguide* directory for this task.
