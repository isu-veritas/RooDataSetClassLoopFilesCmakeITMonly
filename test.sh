#!/bin/bash
#Submit this script with: sbatch thefilename
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-user=achrmy@iastate.edu
#SBATCH --mail-type=END

bin/MLRooDataStore /work/LAS/amandajw-lab/users/achrmy/simulations/V5/ATM21/stage5v257fix/stage5.list grisudetATM21ITM00fix > grisudetATM21ITM00fix.log 2>&1
#bin/MLRooDataStore /work/LAS/amandajw-lab/users/achrmy/simulations/V5/ATM21/stage5v257fix/stage5.list grisudetATM21ITM00v257_0p2 > grisudetATM21ITM00v257_0p2.log 2>&1
#bin/MLRooDataStore stage5zen40.list grisudetATM21ITM40 > grisudetATM21ITM40.log 2>&1
#bin/MLRooDataStore stage5zen30.list grisudetATM21ITM30 > grisudetATM21ITM30.log 2>&1
#bin/MLRooDataStore stage5zen35.list grisudetATM21ITM35 > grisudetATM21ITM35.log 2>&1
