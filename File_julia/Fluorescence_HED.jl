# This file runs an experiment for the continuous measurement of a qubit through heterodyne detection.
# The FI for the continuous measurement, the QFI and the FI for a strong measurement on the conditional state are calculated in order to provide a local estimation of two parameters.
# The function returns a tuple containing time vectors, effective FI and QFI, the Uhlmann curvature and its modulus.


using ZChop 
using Distributed


function Eff_QFI_HED(
    Ntraj::Int64,                       # Number of trajectories
    Tfinal::Float64,                    # Final time
    dt::Float64;                        # Time step
    κ = 1.,                             # Noise coupling
    ωx = 0.,                            # First parameter
    ωy = 0.,                            # Second parameter
    η = 1.)                             # Noise parameter 

    Ntime = Int(floor(Tfinal/dt)) # Number of timesteps

    # Define the possible noise channels 
    c1 = sqrt(κ/2) * (σ(:x) - 1im * σ(:y))/2
    c2 = -1im * sqrt(κ/2) * (σ(:x) - 1im * σ(:y))/2
    
    # ωz and ωy are the parameters that we want to estimate
    H = ωx * σ(:x) / 2 + ωy * σ(:y) / 2   # Hamiltonian of the spin system
    dxH = σ(:x) / 2                 # Derivative of H wrt ωz
    dyH = σ(:y) / 2                 # Derivative of H wrt ωy
    
    
    # Kraus-like operator, trajectory-independent part
    M0 = I - 1im * H * dt - 0.5 * dt * c1' * c1 - 0.5 * dt * c2' * c2

    # Derivative of the Kraus-like operator wrt to the parameters
    dxM = -1im * dxH * dt          
    dyM = -1im * dyH * dt          
    
    # Initial state of the qubit: set in |1⟩ 
    ψ0 = zeros(Complex{Float64}, 2)
    ψ0[1] = 1.
    ψ0[2] = 0. 
    ρ0 = ψ0 * ψ0'
    
    #States for the POVM on the conditional ρ
    ψ1 = zeros(Complex{Float64}, 2)
    ψ2 = zeros(Complex{Float64}, 2)
    ψ3 = zeros(Complex{Float64}, 2)
    ψ4 = zeros(Complex{Float64}, 2)
    ψ1[1] = 1. / sqrt(2.)
    ψ1[2] = 1. / sqrt(2.)
    ψ2[1] = -1. / sqrt(2.)
    ψ2[2] = 1. / sqrt(2.)
    ψ3[1] = 1im / sqrt(2.)
    ψ3[2] = 1. / sqrt(2.)
    ψ4[1] = -1im / sqrt(2.)
    ψ4[2] = 1. / sqrt(2.)

    t = (1 : Ntime) * dt

    # Run evolution for each trajectory, and build up the average for FI and final strong measurement QFI
    result = @distributed (+) for ktraj = 1 : Ntraj
        ρ = ρ0 # Assign initial state to each trajectory
        dxρ = zero(ρ) # Derivative of ρ wrt the parameters, the initial state does not depend on the paramter
        dyρ = zero(ρ)
        τx = dxρ
        τy = dyρ
        dW1() = sqrt(dt) * randn()              #Wiener increment on the σx measure
        dW2() = sqrt(dt) * randn()              #Wiener increment on the σy measure
        
        # Matrixes for the FI and QFI for each time
        FisherT = [([0. 0.;0. 0.]) for i=1:Ntime]
        FisherCondT = [([0. 0.;0. 0.]) for i=1:Ntime]
        QFisherT = [([0. 0.;0. 0.]) for i=1:Ntime]
	Dxy = [0. for i=1:Ntime]
        absDxy  = [0. for i=1:Ntime]

        for jt = 1 : Ntime
            
            dy1 = dt * sqrt(η) * tr(ρ * (c1+c1')) + dW1()     # Current for σx 
            dy2 = dt * sqrt(η) * tr(ρ * (c2+c2')) + dW2()     # Current for σy 
            
            # Euler-Milstein correction
            C = c1 * c1 * (dy1^2 - dt) + (c1 * c2 + c2 *c1) * dy1 * dy2 + c2 * c2 * (dy2^2 - dt)

            # Kraus operator 
            M = M0 + sqrt(η) * c1 * dy1 + sqrt(η) * c2 * dy2 + η/2 * C

            # Evolve the density operator
            new_ρ = M * ρ * M' + (1 - η) * dt * c1 * ρ * c1' + (1 - η) * dt * c2 * ρ * c2' 

            zchop!(new_ρ) # Round off elements smaller than 1e-14

            tr_ρ = tr(new_ρ)
            
            # Evolve the unnormalized derivative wrt ωz and ωy
            τx = (M * (τx * M' + ρ * dxM') + dxM * ρ * M' + (1 - η)* dt * c1 * τx * c1' + (1 - η)* dt * c2 * τx * c2')/ tr_ρ
            τy = (M * (τy * M' + ρ * dyM') + dyM * ρ * M' + (1 - η)* dt * c1 * τy * c1' + (1 - η)* dt * c2 * τy * c2')/ tr_ρ

            zchop!(τx) # Round off elements smaller than 1e-14
            zchop!(τy)

            tr_τx = tr(τx)
            tr_τy = tr(τy)
            
            # Now we can renormalize ρ and its derivative wrt ω1 and ω2
            ρ = new_ρ / tr_ρ
            dxρ = τx - tr_τx * ρ
            dyρ = τy - tr_τy * ρ

            # We evaluate the classical FI for the continuous measurement
            FisherT[jt] = ([real(tr_τx * tr_τx)  real(tr_τx * tr_τy); real(tr_τx * tr_τy)  real(tr_τy * tr_τy)])

            # We evaluate the classical FI for a POVM on the conditional state
            Lx = L(ρ,dxρ)
            Ly = L(ρ,dyρ)
            P1 = real(ψ1' * ρ * ψ1)/2
            dxP1 = real(tr(ρ * ψ1 * ψ1' * Lx))/2
            dyP1 = real(tr(ρ * ψ1 * ψ1' * Ly))/2
            P2 = real(ψ2' * ρ * ψ2)/2
            dxP2 = real(tr(ρ * ψ2 * ψ2' * Lx))/2
            dyP2 = real(tr(ρ * ψ2 * ψ2' * Ly))/2
            P3 = real(ψ3' * ρ * ψ3)/2
            dxP3 = real(tr(ρ * ψ3 * ψ3' * Lx))/2
            dyP3 = real(tr(ρ * ψ3 * ψ3' * Ly))/2
            P4 = real(ψ4' * ρ * ψ4)/2
            dxP4 = real(tr(ρ * ψ4 * ψ4' * Lx))/2
            dyP4 = real(tr(ρ * ψ4 * ψ4' * Ly))/2
            FisherCondT[jt] = ([dxP1^2 / P1 + dxP2^2 / P2 + dxP3^2 / P3 + dxP4^2 / P4    dxP1 * dyP1 / P1 + dxP2 * dyP2 / P2 + dxP3 * dyP3 / P3 + dxP4 * dyP4 / P4;
                                dxP1 * dyP1 / P1 + dxP2 * dyP2 / P2 + dxP3 * dyP3 / P3 + dxP4 * dyP4 / P4    dyP1^2 / P1 + dyP2^2 / P2 + dyP3^2 / P3 + dyP4^2 / P4])
            
            # We evaluate the QFI for a final strong measurement done at time t
            Qxx = QFI(ρ, dxρ, dxρ)
            Qxy = QFI(ρ, dxρ, dyρ)
            Qyy = QFI(ρ, dyρ, dyρ)
            QFisherT[jt] = ([Qxx Qxy; Qxy Qyy])
            
            # We evaluate the Uhlmann curvature for the non-diagonal element
            A = Lx * Ly - Ly * Lx
            Dxy[jt] = real(- 1im / 2 * tr(ρ * A))
            absDxy[jt] = abs(Dxy[jt])
            
        end
        
        hcat(FisherT, FisherCondT, QFisherT, absDxy) 
        
    end

    return (t, result[:,1] / Ntraj, result[:,2] / Ntraj, result[:,3] / Ntraj, result[:,4] / Ntraj)
    
end