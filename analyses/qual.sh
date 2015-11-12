#!/bin/bash

#PBS -N qualMiSeq
#PBS -S /bin/bash
#PBS -l nodes=1:ppn=4
#PBS -o qual.out
#PBS -e qual.err
#PBS -V

cd ~/CommonG/MiSeq
for ffile in *_trun.fastq; 
	do QualFilterFastq.pl $ffile 20 18 ${ffile/_trun.fastq/_qual.fastq};
done  

