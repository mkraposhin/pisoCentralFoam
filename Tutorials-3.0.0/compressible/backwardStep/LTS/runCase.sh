#!/bin/bash

blockMesh

rm -rf log

cat system/controlDict.firstIterations > system/controlDict; pisoCentralFoam | tee -a log

cat system/controlDict.otherIterations > system/controlDict; pisoCentralFoam | tee -a log

