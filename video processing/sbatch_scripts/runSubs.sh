#!/bin/bash
#SBATCH -N 1
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --mem=64000
#SBATCH --array=1-3

export MATLABPATH=$'/MouseMotionMapper/'
cd $MATLABPATH

listOfFiles=""
SAVEPATH=/MouseMotionMapper/trainings/

rowNumber=$SLURM_ARRAY_TASK_ID
fileName=$(sed "${rowNumber}q;d" $listOfFiles)

module load matlab/R2013a
matlab -nosplash -nodesktop -nodisplay -singleCompThread -r "addpath(genpath('$MATLABPATH')); runSubsamples('$fileName',20,.25,20,'$SAVEPATH'); exit;"

