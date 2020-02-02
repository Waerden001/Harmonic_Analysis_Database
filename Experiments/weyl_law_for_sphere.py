from scipy.special import gamma, factorial, binom as gamma, factorial, binom

def spherical_dim(degree, dimension):
    """Dimension of degree d homogeneous polynomials with dimension n variables restricted to S^{n-1}.

    """
    return binom(degree+dimension-1, dimension -1)-binom(degree+dimension-3, dimension -1)

def unit_ball_volume(dimension):
    """This only computes the constant coefficient of the volume formula.
       e.g. 2, 1 (\piR^2), 4/3 (\pi R^3), 1/2 (\pi^2 R^4) 
       The reason is that the \pi and R factors are not important in the actually computation
       we only need to know that the missing part is given by \pi^{floor(dimension)}R^{dimension}
    """
    if dimension % 2 == 0:
        return 1/factorial(dimension/2)
    if dimension % 2 == 1:
        count = dimension/2
        val = 1
        while count >0:
            val *= count
            count -= 1
        return 1/val
    

def unit_sphere_area(dimension):
    """Similarly, this only computes the coefficient. The missing part is given by \pi^{floor(dimension+1)}R^{dimension}.

    """
    return unit_ball_volume(dimension+1)*(dimension+1)

def weyl_law(bound, dim):
    """Note that the dim parameter of the spherical_dim function is the dimension of the ambient space 
       while the dim parameter of weyl_law is the dimension of our sphere.

       Weyl law says that N(bound)/bound^{dimension/2} is given by (2\pi)^{dim}*Vol(unit_ball)*Vol(Manifold) asymptotically. 
       On spheres, eigenvalues of  -Laplacian are degree(degree+dim+1-2), given by all degree d homogeneous polynomials in dim+1
       variables restricted to the sphere. 

       This description works for S^2 or higher dimensions.
    """
    assert dim >= 2, "Only works for S^2 or higher, do it yourself for S^1!"
    degree = 0
    count = 0
    while(degree*(degree+dim+1-2)<bound):
        count += spherical_dim(degree, dim + 1)/bound**(dim/2)
        degree += 1
    return count*(2**dim)/(unit_ball_volume(dim)*unit_sphere_area(dim))

"""------------------------------------------------------TEST--------------------------------------------------------------------------"""
for dim in range(2,10):
    print(weyl_law(100000000, dim))