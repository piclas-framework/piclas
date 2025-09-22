(sec:padaption)=
# Variable polynomial degree (p-adaption)

The field solver of PICLas for Maxwell's and Poisson's equations supports a variable, element-local polynomial degree. This functionality
can be enabled with

    pAdaptionType = non-periodic-BC

The currently available options are:

| Option                | Description                                                                                                                                                                |
| --------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `none` (0)            | Default for setting all elements to $N$                                                                                                                                    |
| `random` (1)          | Elements get random polynomial degree between $N_{\mathrm{min}}$ and $N_{\mathrm{max}}$                                                                                    |
| `non-periodic-BC` (2) | Refinement of elements with non-periodic boundary conditions, see Section {ref}`sec:padaption-BC`                                                                          |
| `half-half` (3)       | Elements in the lower half domain in x-direction are set to $N_{\mathrm{min}}$ and the upper half are set to $N_{\mathrm{max}}$. The domain must be centered around $x=0$. |

An option to adapt the polynomial degree based on the Debye length is currently in development.

(sec:padaption-BC)=
## Refining boundary elements (pAdaptionType = non-periodic-BC)

To refine elements at a non-periodic boundary condition, which can be enabled by

    pAdaptionBCLevel = 1st-and-2nd-NMin+1

Several options are available:

| Option                           | Description                                                                                                                      |
| -------------------------------- | -------------------------------------------------------------------------------------------------------------------------------- |
| `1st-and-2nd-NMin+1` (-2)        | Elements with non-periodic boundary conditions receive $N_{\mathrm{max}}$, 2$^{\mathrm{nd}}$ layer receives $N_{\mathrm{min}}+1$ |
| `directly-connected-NMin+1` (-1) | Elements with non-periodic boundary conditions receive $N_{\mathrm{min}}+1$                                                      |
| `directly-connected-NMax` (1)    | Elements with non-periodic boundary conditions receive $N_{\mathrm{max}}$                                                        |
| `1st-and-2nd-NMax`  (2)          | First two layers of elements with non-periodic boundary conditions receive $N_{\mathrm{max}}$                                    |

An example is available in the regression tests, e.g. in `regressioncheck/NIG_PIC_poisson_Leapfrog_single_node/box_VDL_and_linPhi`