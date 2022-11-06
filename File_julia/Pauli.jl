# This file contains the function 
#
# σ(direction) 	
#
# that returns one of the Pauli matrices dipending on the direction given as input.
# The symbols accepted as entry are x, y and z.		
# Pauli matrices are used in the code to build noise operators for the dechoerence dynamics.


function σ(direction::Symbol)
    @assert direction ∈ (:x, :y, :z) "Direction must be :x, :y, :z"

    sigma = Dict(:x => ([0im 1.; 1. 0.] ),
        :y => ([0. -1im; 1im 0.]),
        :z => ([1. 0im; 0. -1.]))
    
    return sigma[direction]
end
