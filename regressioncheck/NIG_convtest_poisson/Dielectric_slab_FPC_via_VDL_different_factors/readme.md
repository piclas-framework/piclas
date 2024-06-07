# Dielectric + FPC modelled via VDL
- This test does NOT require PETSC for FPC because the FPC is not used and instead the VDL model is applied (set PARTICLES=ON and PICLAS_CODE_ANALYZE=ON to successfully run)
- The actual dielectric region from -10e-9m to 0 with epsR=10 is not considered but replaced by a virtual dielectric layer (VDL)
  model
  - it shifts impacting charged particles from the real interface between the dielectric and the vacuum region to a different position
  - the position is determined by keeping the ratio between the thickness and the permittivity epsR constant
  - the resulting position that yields a new permittivity of eprR=1 (vacuum) is used, therefore, eliminating the normal dielectric
    model that is used in piclas
- comparison with the analytical solution in 1D
- Convergence test with N=1 for analyze_Convtest_h_cells=1,2,4,8,16,32,64 yield an order of convergence of O(0.5) because the jump
  between the dielectric and the vacuum region is always within an element and not exactly at the interface
- Different stretching factors are tested with different meshes
  - factor 10: 1D_mortonZ_48_90e-9_mesh.h5
  - factor 100: 1D_mortonZ_48_990e-9_mesh.h5
  - factor 1000: 1D_mortonZ_48_9990e-9_mesh.h5
