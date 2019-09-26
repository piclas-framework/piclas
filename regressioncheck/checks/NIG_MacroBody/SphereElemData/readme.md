# Macro sphere Elemdata analyze
*  Macro sphere moves through channel
*  2 meshes with Bilinear sides are used
   *  100x10x10 and 10x1x1 and volume data is compared to reference files
*  Additionally DSMC particles are inserted and tracked
*  Tests:
  *  Volume portions are calculated corrrectly
  *  MacroPartInElem check assigns elements that contain macro sphere for tracking correctly
  *  Particles are tracked without moving out of channel or inside of spheres

