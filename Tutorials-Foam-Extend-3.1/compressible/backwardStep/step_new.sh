#PBS -l walltime=04:30:00,nodes=1:ppn=12

cd PATH/TO/YOUR/DIR

rm -rf log

blockMesh | tee -a log
decomposePar -force | tee -a log
mpirun -np 12 -machinefile $PBS_NODEFILE pisoCentralFoam -parallel | tee -a log
reconstructPar -latestTime
sample

