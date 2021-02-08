#!/bin/bash
#SBATCH -N 1
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32000
#SBATCH --array=1-3

export MATLABPATH=$'/MouseMotionMapper/'
cd $MATLABPATH

## The inputs to makeHDK_proj are the following
########################
# PCA projection file list 
listOfFiles=""

# K-means clusters output
Kmeans=''

# Training set
TrainingSet='/MouseMotionMapper/demo/trainingSet_new10.mat'

# Save Path
SavePath='MouseMotionMapper/Reembedding/'

# The end of output-file names
EndName='_p10_HD_101.mat'

# Number of PCA used (20 in our case)
NumPCA=10
#################

rowNumber=$SLURM_ARRAY_TASK_ID
fileName=$(sed "${rowNumber}q;d" $listOfFiles)

echo $listOfFiles
echo $rowNumber
echo $fileName


module load matlab/R2013a
matlab -nosplash -nodesktop -nodisplay -singleCompThread -r "addpath(genpath('$MATLABPATH'));\
	makeHDK_proj('$fileName',$Kmeans,$TrainingSet,$SavePath,$EndName,$NumPCA); exit;"