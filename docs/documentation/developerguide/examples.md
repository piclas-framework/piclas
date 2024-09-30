# Markdown Examples
## hyperlinks
[one with a title](http://fsf.org "click here for a good time!"). Unclear how to enforce new window.

## Code environment
Either use fenced style (tildes) 

~~~~~~~
if (a > 3) {
  moveShip(5 * gravity, DOWN);
}
~~~~~~~

or indented style (4 whitespaces)

    if (a > 3) {
      moveShip(5 * gravity, DOWN);
    }

Both works with pandoc and wordpress. Also see [pandoc verbatim code](http://pandoc.org/README.html#verbatim-code-blocks "pandoc verbatim code").

## Equations
(@gleichung1) $$a=b*c$$
As (@gleichung1) shows, blabla.

## Bibtex, cite
Hindenlang [@Hindenlang2015]. Only works with pandoc!

[bibshow file=references.bib]

Hindenlang [bibcite key=Hindenlang2015], Gassner [bibcite key=gassner2011disp]


## section references
## Figures, caption

```{figure} https://github.com/piclas-framework/piclas/blob/master/docs/logo.png?raw=true
---
name: fig:mylabel
width: 400px
align: center
---

This is an example caption.
```
See {numref}`fig:mylabel` for an image from the web embedded in this documentation.

```{figure} figures/mpi_shared_mesh/dev_mpi_shared_mesh.png
---
name: fig:example
width: 200px
align: center
---

This is an example caption.
```
See {numref}`fig:example` for embedding a local file.

## tables
## unnumbered section headings
  just add

    {-}

 after the heading

## Code blocks for various languages

```{code-block} C

int a = 32;
int a = 32;

```
