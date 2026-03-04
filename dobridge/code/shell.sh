#!/bin/bash
#SBATCH --job-name=dobridge
#SBATCH --output=dobridge_%j.out
#SBATCH --error=dobridge_%j.err
#SBATCH --mail-user=logan_casey@berkeley.edu

stata -b do clean_data.do

if [ $? -ne 0 ]; then
    echo "Stata job failed"
    exit 1
else
    echo "Stata job completed successfully"
fi