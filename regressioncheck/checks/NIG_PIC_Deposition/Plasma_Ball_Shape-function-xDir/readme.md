# plasma ball example
- fully 3D periodic domain for testing different shape function methods
- one element in z-direction and 4 in x- and y-direction
- A ball of 100000 particles is initialized and deposited on the grid
    in order to check the deposition method "shape_function" in different
    directions via

    PIC-shapefunction-direction = 1,3

    and depositing the charge one- and two-dimensionally via

    PIC-shapefunction-dimension = 1,2
- 1D and x-direction means variation in x and constant in all other directions
- 2D and z-direction means constant in z-direction and variation in x- and y-direction
