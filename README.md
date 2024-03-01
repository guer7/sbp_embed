# sbp_embed
Repository for code using summation-by-parts (SBP) embedding operators. 

## The two-dimensional Maxwell equations
The directory 'Maxwell' contains code for the 4-block 2D Maxwell simulation with interface conditions imposed using the SBP embedding method, described in the paper 'Operator-based analysis of summation-by-parts methods'.

## Laplace operator
The directory 'Laplace' contains code producing the numerical results in 'Efficient discretization of the Laplacian on complex geometries'. Included in the directory is sbplib, https://sourceforge.net/projects/sbplib/, and code to construct glue-grid interpolation operators presented in 'J.E. Kozdon, L.C. Wilcox, Stable coupling of nonconforming, high-order finite difference methods, SIAM J. Sci. Comput. (2016) A923–A952, https://doi.org/10.1137/15M1022823.'
