#!/bin/bash

SCRIPT_PATH=`pwd`

if [ $# -eq 2 ]
then
    SCRIPT_PATH=$2
    cd $1
fi

CFILES=`ls -a`

for CFILE in $CFILES
do
    if [ -d $CFILE ]
    then
        if [ $CFILE != "." -a $CFILE != ".." ]
        then
            $SCRIPT_PATH/cleanCases.sh $CFILE $SCRIPT_PATH
        fi
    else
        if [ $CFILE == "cleanCase.sh" ]
        then
            echo "cleaningCase $PWD"
            ./cleanCase.sh
        fi
    fi
done

if [ $# -eq 2 ]
then
    cd ../
fi

#
#END OF FILE
#


