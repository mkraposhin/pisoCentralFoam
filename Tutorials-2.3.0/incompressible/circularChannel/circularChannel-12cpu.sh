#PBS -l walltime=240:30:00,nodes=1:ppn=12

cd PATH/TO/YOUR/DIR

rm -rf log

blockMesh

decomposePar | tee -a log
mpirun -np 12 -machinefile $PBS_NODEFILE pisoCentralFoam -parallel | tee -a log
reconstructPar | tee -a log
sample | tee -a log

