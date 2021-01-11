# plasma ball example
- fully 3D periodic domain for testing different shape function methods
- one element in x-direction and 4 in y- and z-direction
- A ball of 100000 particles is initialized and deposited on the grid
    in order to check the deposition method "shape_function" in different
    directions via

    PIC-shapefunction-direction = 2,1

    and depositing the charge one- and two-dimensionally via

    PIC-shapefunction-dimension = 1,2
- 1D and y-direction means variation in y and constant in all other directions
- 2D and x-direction means constant in x-direction and variation in y- and z-direction
