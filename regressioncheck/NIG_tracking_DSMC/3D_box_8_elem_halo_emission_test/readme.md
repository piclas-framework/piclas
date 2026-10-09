# 3D_box_8_elem_halo_emission_test
- particle emission test, which makes sure that a particles is emitted only once when the coordinates of it are exactly on the MPI
interface between two adjacent MPI processes
- 1 particle is emitted in the corner of 8 connected elements
- the test is successful if only one particle is created (without the check in the halo emission routines, this test produces 8
particles when it is executed with 8 processes)