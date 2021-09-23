# plasma ball example
- fully 3D periodic domain for testing different shape function methods
- one element in y-direction and 4 in x- and z-direction
- A ball of 100000 particles is initialized and deposited on the grid
    in order to check the deposition method "shape_function" in different
    directions via

    PIC-shapefunction-direction = 3,2

    and depositing the charge one- and two-dimensionally via

    PIC-shapefunction-dimension = 1,2
- 1D and z-direction means variation in z and constant in all other directions
- 2D and y-direction means constant in y-direction and variation in x- and z-direction
