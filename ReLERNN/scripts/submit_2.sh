#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --time=10:00:00
#SBATCH --exclusive
#SBATCH --partition gpu
#SBATCH --nodelist node94
#SBATCH -o rj_relernn.%N.%j.out # STDOUT
#SBATCH -e rj_relernn.%N.%j.err # STDERR


conda deactivate

module load python/3.10.14
module load cuda/11.8
module load cudnn/8.9.7.29_cuda11.x

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scratch/bioconsult/Maddy_Bursell/TensorRT-8.6.1.6/lib

python3 -c "import tensorflow as tf; print(tf.config.list_physical_devices('GPU'))"

./redjg_full.sh
