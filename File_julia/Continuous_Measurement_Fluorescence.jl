# Fisher information for a continuously monitored qubit in a spontaneous emission context.
# Different detecting schemes are tested, namely:
# Fluorescence_HD	> homodyne detection
# Fluorescence_PD	> photodetection
# Fluorescence_HED	> heterodyne detection

	using Distributed
	
	addprocs(4)

	@everywhere include("Pauli.jl")
	@everywhere include("SLD.jl")
	@everywhere include("QFI.jl")
	@everywhere include("Fluorescence_HD.jl")
	@everywhere include("Fluorescence_PD.jl")
	@everywhere include("Fluorescence_HED.jl")

	function Eff_QFI(
   			Ntraj::Int64,		# Number of trajectories
    			Tfinal::Float64,	# Final time
    			dt::Float64;		# Time step
              		)

		hd = Eff_QFI_HD(Ntraj, Tfinal, dt)
		pd = Eff_QFI_PD(Ntraj, Tfinal, dt)
		he = Eff_QFI_HED(Ntraj, Tfinal, dt)

		Ntime = Int(floor(Tfinal/dt))

		output_hd = open("Output_Fluorescence_HD.dat","w")
		for jt = 1 : Ntime
			show(output_hd, hd[1][jt])
			write(output_hd, "\t")
			show(output_hd, hd[2][jt][1])
			write(output_hd, "\t")
			show(output_hd, hd[2][jt][2])
			write(output_hd, "\t")
			show(output_hd, hd[2][jt][4])
			write(output_hd, "\t")
			show(output_hd, hd[3][jt][1])
			write(output_hd, "\t")
			show(output_hd, hd[3][jt][2])
			write(output_hd, "\t")
			show(output_hd, hd[3][jt][4])
			write(output_hd, "\t")
			show(output_hd, hd[4][jt][1])
			write(output_hd, "\t")
			show(output_hd, hd[4][jt][2])
			write(output_hd, "\t")
			show(output_hd, hd[4][jt][4])
			write(output_hd, "\t")
			show(output_hd, hd[5][jt])
			write(output_hd, "\n")
		end
		close(output_hd)

		output_pd = open("Output_Fluorescence_PD.dat","w")
		for jt = 1 : Ntime
			show(output_pd, pd[1][jt])
			write(output_pd, "\t")
			show(output_pd, pd[2][jt][1])
			write(output_pd, "\t")
			show(output_pd, pd[2][jt][2])
			write(output_pd, "\t")
			show(output_pd, pd[2][jt][4])
			write(output_pd, "\t")
			show(output_pd, pd[3][jt][1])
			write(output_pd, "\t")
			show(output_pd, pd[3][jt][2])
			write(output_pd, "\t")
			show(output_pd, pd[3][jt][4])
			write(output_pd, "\t")
			show(output_pd, pd[4][jt][1])
			write(output_pd, "\t")
			show(output_pd, pd[4][jt][2])
			write(output_pd, "\t")
			show(output_pd, pd[4][jt][4])
			write(output_pd, "\t")
			show(output_pd, pd[5][jt])
			write(output_pd, "\n")
		end
		close(output_pd)

		output_he = open("Output_Fluorescence_HED.dat","w")
		for jt = 1 : Ntime
			show(output_he, he[1][jt])
			write(output_he, "\t")
			show(output_he, he[2][jt][1])
			write(output_he, "\t")
			show(output_he, he[2][jt][2])
			write(output_he, "\t")
			show(output_he, he[2][jt][4])
			write(output_he, "\t")
			show(output_he, he[3][jt][1])
			write(output_he, "\t")
			show(output_he, he[3][jt][2])
			write(output_he, "\t")
			show(output_he, he[3][jt][4])
			write(output_he, "\t")
			show(output_he, he[4][jt][1])
			write(output_he, "\t")
			show(output_he, he[4][jt][2])
			write(output_he, "\t")
			show(output_he, he[4][jt][4])
			write(output_he, "\t")
			show(output_he, he[5][jt])
			write(output_he, "\n")
		end
		close(output_he)
	
	end

	@time Eff_QFI(200, 1., 0.001)


