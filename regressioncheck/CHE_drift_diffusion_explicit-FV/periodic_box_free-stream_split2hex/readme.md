# FV Drift-diffusion Free Stream
- most basic test with constant initial conditions, which must not change
- periodic box
- split2hex settings in hopr.ini with 2x2x2 elements that are split into 768 hex elements

    elemtype       = 104                         ! element type (108: Hexahedral)
    meshTemplate   = 3
    SplitToHex     = T
    nFineHexa      = 1
