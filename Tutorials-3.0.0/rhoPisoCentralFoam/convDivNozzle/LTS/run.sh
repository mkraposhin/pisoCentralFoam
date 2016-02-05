#!/bin/bash

rm -rf log

blockMesh | tee -a log
rhoPisoCentralFoam | tee -a log
sample | tee -a log
