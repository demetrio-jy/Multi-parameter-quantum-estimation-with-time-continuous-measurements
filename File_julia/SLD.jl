# This file contains a function 
#
# L(ρ, dρ; abstol)
#
# that returns the symmetric logarithmic derivative operator in the computational basis.
# The arguments of the function are:  
# ρ	> density matrix
# dρ	> derivative of the density matrix wrt a parameter 


using ZChop
using LinearAlgebra


function L(ρ, dρ; abstol = 1e-5)
    L = zero(ρ)
    eigval, eigvec = eigen(Matrix(zchop(ρ,1e-10)))
    eigval = real(eigval)
    dim = length(eigval)
    for n=1:dim, m=1:dim
        L[n,m] = 2 * ((eigval[n] + eigval[m] > abstol) ? (1. / (eigval[n] + eigval[m])) * 
                    (eigvec[:,n]' * dρ * eigvec[:,m]) : 0.) 
    end
    return eigvec * L * inv(eigvec)
end
