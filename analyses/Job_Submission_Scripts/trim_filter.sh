#!/bin/bash

#PBS -N trim_filMiSeq
#PBS -S /bin/bash
#PBS -l nodes=1:ppn=4
#PBS -o trimfil.out
#PBS -e trimfil.err
#PBS -V

cd ~/CommonG/MiSeq
for ffile in *.fastq; 
	do TruncateFastq.pl $ffile 1 36 ${ffile/.fastq/_trun.fastq};
done  
