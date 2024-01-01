#!/bin/bash -l


#-----The '#' followed by SBATCH are SLURM prefixes. They are not comments-----#
#-----It is worth taking note of them when you write your own batch files------#


#SBATCH -D ./
#SBATCH --export=ALL
#SBATCH -p course
#SBATCH -t 5


SRC1=$1
SRC2=$2
SRC3=$3
#ADD INSTRUCTIONS TO LOAD THE MODULES HERE
module load compilers/intel/2019u5
module load mpi/intel-mpi/2019u5/bin


#ADD COMPILER INSTRUCTION HERE.
mpicc -fopenmp -O3 main-mpi-only.c coordReader.c ompnAddition.c -lm


#SLURM_NTASKS is given by the -n flag when using sbatch.
procs=${SLURM_NTASKS:-1}


#RUN THE PROGRAM HERE USING $procs WITH THE -np FLAG.
mpirun -np 2 ./a.out 9_coords.coord
