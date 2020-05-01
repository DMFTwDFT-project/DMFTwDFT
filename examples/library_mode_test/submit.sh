#!/bin/bash
#SBATCH --job-name=UO2
#SBATCH -N 4
#SBATCH -p bdwall
##SBATCH -A Meso 
##SBATCH -A dmft_for_oxides    
#SBATCH -A nickelates         
#SBATCH --ntasks-per-node=36
#SBATCH -t 40:00:00
#SBATCH -D ./ 

export I_MPI_FABRICS=shm:tmi

export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_F90=ifort

module load intel/17.0.4-74uvhji #intel/16.0.4-nzcw5zc
module load gsl/2.4
module load anaconda/4.4.0 
module load cmake/3.8.1-orygmpj
export PYTHONPATH=/home/hyowon/Codes/WANDMFT/bin:$PYTHONPATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/soft/bebop/gsl/2.4/lib

mkdir $SLURM_JOBID
cd $SLURM_JOBID
cp ../* .
ulimit -s unlimited
mpirun -n $SLURM_NTASKS ./test.x > out 2> err 
#./dmft.x > out 2> err 

