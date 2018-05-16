#! /bin/bash
#PBS -N tt8
#PBS -q gpu
#PBS -l select=1:ngpus=1

cd $PBS_O_WORKDIR
/hpctmp/a0132576/apps/mmatlab2015b/bin/matlab -nodisplay -r main> output_file.txt

