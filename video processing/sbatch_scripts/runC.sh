#!/bin/bash
#SBATCH -N 1
#SBATCH --time=54:00:00
#SBATCH --array=1
#SBATCH --mem=64000

export MATLABPATH=$'/MouseMotionMapper/'
cd $MATLABPATH

TRAINDATA=/MouseMotionMapper/demo/trainingSet_new10.mat
SAVEPATH=/MouseMotionMapper/Kmeans/
module load matlab/R2013a
matlab -nosplash -nodesktop -nodisplay -singleCompThread -r "addpath(genpath('$MATLABPATH')); runCluster('101','$SAVEPATH','$TRAINDATA'); exit;"

