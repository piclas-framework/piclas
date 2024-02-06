# Particle-In-Cell

(sec:PIC-deposition)=
## Charge and Current Deposition

Charge and current deposition can be performed using different methods, among others, shape
functions, B-splines or locally volume-weighted approaches.

|  **PIC-Deposition-Type**  |                             **Description**                            |
|      :--------------:     |                  :-----------------------------------:                 |
|   *cell_volweight_mean*   | Linear distribution in each element (continuous on element boundaries) |
|      *shape_function*     |                standard shape function with fixed radius               |
|    *shape_function_cc*    |            charge corrected shape function with fixed radius           |
| *shape_function_adaptive* |      charge corrected shape function with element-dependent radius     |

### Linear Distribution Over Cell Interfaces
A linear deposition method that also considers neighbouring elements can be selected by

    PIC-Deposition-Type = cell_volweight_mean

and is referred to as the CVWM method.
This method also considers the corner nodes of each element to which all neighbouring elements
contribute, hence, resulting in a non-local deposition scheme.
Note that the CVWM method allows switching of charge deposition on Dirichlet boundaries via

    PIC-DoDirichletDeposition = F

which simply nullifies the deposited charge on wall boundary nodes for Dirichlet sides to account for mirror charges.
The default value for this parameter is true and it is currently only available for the CVWM method in combination with the HDG
method.

### Shape Function

High-order field solvers require deposition methods that reduce the noise, e.g., shape functions {cite}`Jacobs2006`. The standard 3D shape function is selected by

    PIC-Deposition-Type = shape_function

or

    PIC-Deposition-Type = shape_function_cc

where `shape_function_cc` is a numerically charge-conserving method that adjusts the deposited charge by
comparing its integral value to the total charge given by the particles.

The shape function sphere might be truncated at walls or open boundaries, which is compensated when using `shape_function_cc` by
increasing the deposited charge of truncated particles.

Additionally, an element-local shape function radius can be used, which is determined for each element separately depending on the
size of the element and its direct neighbours by setting

    PIC-Deposition-Type = shape_function_adaptive

The shape function radius in this case is limited by the size of the surrounding elements and may not reach past its direct
neighbours.

The direct influence of only the neibouring elements can be extended further by activating

    PIC-shapefunction-adaptive-smoothing = T

which increases the radius of influence and therefore takes more elements into account for the calculation of the shape function
radius in each element, hence, leading to a smoother transition in regions, where the element sizes rapidly change.

This shape function method also is numerically charge conserving by integrating each particle's deposited charge and
adjusting to this value. Depending on the polynomial degree N, the number of DOF that are within the shape function radius can be
changed via

    PIC-shapefunction-adaptive-DOF = 33

The default values (maximum allowed for each polynomial degree $N$) depend on the dimensionality of the deposition kernel,
1D: $2(N+1)$, 2D: $\pi(N+1)^2$, 3D: $(4/3)\pi(N+1)^3$.

The following polynomial isotropic shape functions are all designed to be used in three dimensions, where reductions to 2D and 1D
are possible.

#### Shape Function 1D
A one-dimensional shape function in $x$-direction is given by

$$
S_{1D}(r,R,\alpha)=\frac{\Gamma(\alpha+3/2)}{\sqrt{\pi}R\Gamma(\alpha+1)\Delta y \Delta z}\left( 1-\left( \frac{r}{R} \right)^{2} \right)^{\alpha}~,
$$

which is normalized to give $\int_{z_{1}}^{z_{2}}\int_{y_{1}}^{y_{2}}\int_{-R}^{R}S_{1D}(r,R,\alpha)dxdydz=1$,
where the radius ${r=|\boldsymbol{x}-\boldsymbol{x}_{n}|=|x-x_{n}|}$ is the distance between the position of the
grid point at position $\boldsymbol{x}$ and the $n$-th particle at position $\boldsymbol{x}_{n}$,
$R$ is the cut-off radius, $\Delta y=y_{2}-y_{1}$ and $\Delta z=z_{2}-z_{1}$ are the domain lengths in $y$- and $z$-direction,
respectively, and $\Gamma(z)$ is the gamma function given by

$$
  \Gamma(z)=\int_{0}^{\infty}x^{z-1}\exp(-x)dx~.
$$

The direction in which deposition is performed is chosen via

    PIC-shapefunction-direction = 1 ! for x-direction
                                  2 ! for y-direction
                                  3 ! for z-direction

and the dimensionality of the shape function is controlled by

    PIC-shapefunction-dimension = 1 ! for 1D
                                  2 ! for 2D
                                  3 ! for 3D

which has to be set to 1 in the case of 1D deposition.

#### Shape Function 2D
A two-dimensional shape function in $x$-$y$-direction is given by

$$
S_{2D}(r,R,\alpha)=\frac{\alpha+1}{\pi R^{2} \Delta z}\left( 1-\left( \frac{r}{R} \right)^{2} \right)^{\alpha}~,
$$

which is normalized to give $\int_{z_{1}}^{z_{2}}\int_{0}^{2\pi}\int_{0}^{R}S_{2D}(r,R,\alpha)rdr d\phi d\theta=1$,
where the radius ${r=|\boldsymbol{x}-\boldsymbol{x}_{n}|}$ is the distance between the position of the
grid point at position $\boldsymbol{x}$ and the $n$-th particle at position $\boldsymbol{x}_{n}$,
$R$ is the cut-off radius and $\Delta z=z_{2}-z_{1}$ is the domain length in $z$-direction.
The perpendicular direction to the two axes, in which deposition is performed is chosen via

    PIC-shapefunction-direction = 1 ! for const. depo in x-direction
                                  2 ! for const. depo in y-direction
                                  3 ! for const. depo in z-direction

when the charge is to be deposited const. along the $x$- or $y$- or $z$-direction.
If the charge is to be deposited over the area instead of the volume, the flag

    PIC-shapefunction-3D-deposition=F

must be set, which simply sets $\Delta z=1$ for the example described above.
Again, the dimensionality of the shape function is controlled by

    PIC-shapefunction-dimension = 1 ! for 1D
                                  2 ! for 2D
                                  3 ! for 3D

which has to be set to 2 in the case of 2D deposition.

#### Shape Function 3D
A three-dimensional shape function in $x$-$y$-direction is given by {cite}`Stock2012`

$$
S_{3D}(r,R,\alpha)=\frac{\Gamma(\alpha+5/2)}{\pi^{3/2}R^{3}\Gamma(\alpha+1)}\left( 1-\left( \frac{r}{R} \right)^{2} \right)^{\alpha}~,
$$

which is normalized to give $\int_{0}^{\pi}\int_{0}^{2\pi}\int_{0}^{R}S_{2D}(r,R,\alpha)r^{2}\sin(\phi)dr d\phi d\theta=1$,
where the radius ${r=|\boldsymbol{x}-\boldsymbol{x}_{n}|}$ is the distance between the position of the
grid point at position $\boldsymbol{x}$ and the $n$-th particle at position $\boldsymbol{x}_{n}$ and
$R$ is the cut-off radius.

