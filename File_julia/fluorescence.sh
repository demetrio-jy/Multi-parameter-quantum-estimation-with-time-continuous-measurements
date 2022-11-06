    #!/bin/bash
    #SBATCH -A piscia
    #SBATCH -c 1
    #SBATCH --job-name=fluorescence
    #SBATCH --mail-user=giacomo.piscia@studenti.unimi.it
    #SBATCH --mail-type=END

    julia Continuous_Measurement_Fluorescence.jl