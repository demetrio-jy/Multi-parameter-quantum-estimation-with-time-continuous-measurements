# This file contains a function 
#
# QFI(ρ, d1ρ, d2ρ; abstol)
#
# that returns the quantum Fisher information for a given state.
# The arguments of the function are:  
# ρ	> density matrix
# d1ρ	> derivative of the density matrix wrt a first parameter 
# d2ρ	> derivative of the density matrix wrt a second parameter 
#
# We enforce the Hermiticity of ρ so that the algorithm is more efficient and returns real values.


using ZChop
using LinearAlgebra


function QFI(ρ, d1ρ, d2ρ; abstol = 1e-5)
    eigval, eigvec = eigen(Matrix(zchop(ρ,1e-10)))
    eigval = real(eigval)
    dim = length(eigval)
    return (real(sum( [( (eigval[n] + eigval[m] > abstol) ? (1. / (eigval[n] + eigval[m])) * 
                    ((eigvec[:,n]' * d2ρ * eigvec[:,m]) * (eigvec[:,m]' * d1ρ * eigvec[:,n]) + 
                        (eigvec[:,n]' * d1ρ * eigvec[:,m]) * (eigvec[:,m]' * d2ρ * eigvec[:,n])) : 0.) for n=1:dim, m=1:dim])))
end