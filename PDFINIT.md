## Initialization of An Initial Particle Distribution Function (PDF)

Currently, the initial particle distributions are simple functions which do not consider great changes in any
coordinate. Here is a fast HOW-TO of creating any initial particle distribution.

Following information are given analytically: n=n(x),v=v(x). The initial distribution which are non-constant
are evaluated in each element using the DSMC cell_local sampling. The required values are evaluated at the
element-origin. This requires a fine sampling mesh with a sufficient high resolution to sample the PDF. These information
are then stored in the InitPart_DSMCHOState.h5 and used for a restart with the MacroRestartValues, resulting in the
InitPDF_State.h5. From this file, the PartData is copied to the initial file, generated with the coarse computational grid


    h5copy -i InitPDF_State.h5 -s PartData -o Coarse_State.h5 -d PartData


Finally, the last entry of PartInt of Coarse_State.h5 has to be set to the number of all particles in InitPDF_State.h5. Note,
that PICLas is executed on a single core.
