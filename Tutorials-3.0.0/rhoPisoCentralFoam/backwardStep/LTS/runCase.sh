#!/bin/bash

blockMesh

rm -rf log

cat system/controlDict.firstIterations > system/controlDict; rhoPisoCentralFoam | tee -a log

cat system/controlDict.otherIterations > system/controlDict; rhoPisoCentralFoam | tee -a log

