#!/bin/bash

rm -rf log

blockMesh | tee -a log
decomposePar | tee -a log
mpirun -np 4 pisoCentralFoam -parallel | tee -a log

