# Markdown Examples

## Linking

### Hyperlinks
[Here is a hyperlink with a title](http://fsf.org "click here for a good time!") that re-directs to the URL "http://fsf.org".

### Referencing other sections of the documentation
Referencing parts of the documentation is always useful when something that has already been documented is required in another
context. To create a link to a different part of the documentation, use referencing via


    {ref}`label`

where `label` refers to a labelled part of the documentation, either automatically created or set by hand.
The two methods are described in the following.

#### Automatic labelling
The automatic label for this sub-section is created as

    developerguide/examples.html:automatic labelling

and can be used via

    {ref}`developerguide/examples:automatic labelling`

which when used in text gives {ref}`developerguide/examples:automatic labelling` and should point to this section.
Note that duplicate labels must be prevented when using this method.
Therefore, it is possible to set labels by hand.

(sec:my-own-label)=
#### Setting labels by hand
To set a label for a section, write the following code style directly in the line before the section is defined

    (sec:my-own-label)=

and then reference this label in the text as

    {ref}`sec:my-own-label`

which when used in text gives {ref}`sec:my-own-label` and should point to this section.

## Code environment
Either use fenced style (tildes) 

~~~~~~~
if (a > 3) {
  moveShip(5 * gravity, DOWN);
}
~~~~~~~

or indented style (4 white spaces)

    if (a > 3) {
      moveShip(5 * gravity, DOWN);
    }

Both works with pandoc and wordpress. Also see [pandoc verbatim code](http://pandoc.org/README.html#verbatim-code-blocks "pandoc verbatim code").

## Equations
(@gleichung1) $$a=b*c$$
As (@gleichung1) shows, blabla.

## Bibtex, cite
Cications are stored in

    ./piclas/docs/documentation/references.bib

using bibtex format, e.g., 

    @article{mathiaud2022bgk,
      title   = {An ES-BGK model for diatomic gases with correct relaxation rates for internal energies},
      author  = {Mathiaud, Julien and Mieussens, Luc and Pfeiffer, Marcel},
      journal = {arXiv preprint arXiv:2202.10906},
      year    = {2022}
    }

and is accessed via

    {cite}`mathiaud2022bgk`

which results in {cite}`mathiaud2022bgk` when used in text.

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
