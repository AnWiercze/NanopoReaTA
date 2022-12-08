#!/bin/bash
echo "Nextflow started!"
echo $1
echo $2
echo $3
nextflow run $1 -params-file $2 -w $3
