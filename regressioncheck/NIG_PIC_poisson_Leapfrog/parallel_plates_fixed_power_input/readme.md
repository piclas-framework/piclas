# parallel plate testcase
* Single electron is accelerated according to electric field in plate capacitor that is adjusted to yield a fixed input power of
  1e-10 Watt by setting the following parameter

        BoundaryName = BC_right ! Adjust the right BC automatically depending on the actually input power
        BoundaryType = (/2,2/)  ! 2: Dirichlet with automatically adjusted electric potential

        CoupledPowerPotential = (/10. , 1000. , 2000./) ! lower, starting and maximum values for the electric potential at all BoundaryType = (/2,2/) BCs
        CoupledPowerTarget    = 1e-10 ! target power of 1e-10 Watt

* Coupled Power output and adjusted electric potential are compared with reference solution
* parameter.ini for activating coupled Power output

          CalcCoupledPower = T
