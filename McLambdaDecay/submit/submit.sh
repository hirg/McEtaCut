#!/bin/bash
date

# . ./submit.sh

if [ $# -eq 0 ]
then
  sbatch --account star --ntasks=1 --array=5-20 --mem-per-cpu=2000 --time=02:00:00 McLambdaEta.slr
else
  echo "Wrong number of parameters!!"
fi

