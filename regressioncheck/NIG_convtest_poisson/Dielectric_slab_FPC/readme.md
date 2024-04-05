# Dielectric + FPC
- This test requires PETSC (for FPC), PARTICLES=ON and PICLAS_CODE_ANALYZE=ON to successfully run
- Dielectric region from -10e-9m to 0 with epsR=10
- FPC at x=0 to model surface charge on a dielectric interface to vacuum
- comparison with the analytical solution in 1D
- Convergence test with N=1 for analyze_Convtest_h_cells=1,2,4,8,16,32,64 yield and order of convergence of O(2)
- Higher orders than 2 makes the numerical solution match the analytical solution (up to machine precision)
