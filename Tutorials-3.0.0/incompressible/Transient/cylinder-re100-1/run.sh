#PBS -l walltime=240:30:00,nodes=1:ppn=4

#cd PATH/TO/YOUR/DIR

rm -rf log

decomposePar -force | tee -a log

mpirun -np 4 -machinefile $PBS_NODEFILE pisoCentralFoam -parallel | tee -a log




