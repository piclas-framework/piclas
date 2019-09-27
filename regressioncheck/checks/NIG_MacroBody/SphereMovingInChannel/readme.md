# Macro sphere in channel analyze
*  Macro sphere moves through channel
*  Mesh with Bilinear sides is used
    *  10x1x1 and volume portion data is compared to reference files
*  Additionally DSMC particles are inserted and tracked (collismode=0)

*  Tests:
    *  Volume portions are calculated corrrectly
    *  MacroPartInElem check assigns elements that contain macro sphere for tracking correctly
    *  Particles are tracked without moving out of channel or inside of spheres
    *  Volume portion calculation, which is be called each timestep because of moving Sphere

*  Note:
   If reggie fails due to elemdata (volumeportion) analyze a reduction of the mesh size to 8x1x1 might help
   The reason are the RANDOM points used for Monte Carlo volume portion calculation

