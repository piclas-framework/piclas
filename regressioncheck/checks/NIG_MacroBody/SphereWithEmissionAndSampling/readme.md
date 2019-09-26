# Macro sphere with insertion/emission
*  Macro sphere is positioned on surface and partially inside mesh
*  2 cubic meses are used: 10x10x10 and 1x1x1
    *  2 Boundaries: 1 open an 5 reflective
*  Inserted DSMC particles via surfaceflux, cell_local emission and cuboid emission

*  Tests:
    *  particle insertion with macro part in domain
    *  no DSMC particle inserted inside Sphere (num. particles is analyzed and compared with ref)
    *  tests DSMC macro write out because of MacroSphere volumereduction
