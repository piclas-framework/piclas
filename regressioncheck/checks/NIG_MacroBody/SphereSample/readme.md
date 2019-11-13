# Macro sphere sampling analyze
*  Macro sphere fixed in channel
*  Meshes with Bilinear and planar sides is used (10x3x3 Elems)
*  Additionally, DSMC particles are inserted which do not move

*  sampling data is compared to reference files
*  Tests:
    *  DSMC macro write out because of MacroSphere volumereduction
    *  Volume portions are calculated corrrectly?

The comparison of DSMCHO-states is problematic because components of velocity can be almost 0
of order +/-1 and therefore the relative comparison is greater then maximum allowed difference.
Consequently, the comparison needs to be deterministic

