#!/bin/bash
echo "STARTED"
python $1 --bamList $2 --gene $3 --converted_gtf $4 --output_dir $5 