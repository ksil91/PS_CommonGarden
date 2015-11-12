#!/bin/bash

#PBS -N adpMiSeq
#PBS -S /bin/bash
#PBS -l nodes=1:ppn=4
#PBS -o adap.out
#PBS -e adap.err
#PBS -V

cd ~/CommonG/MiSeq
for ffile in *_qual.fastq; 
	do AdaptorFilterFastq.pl $ffile adapt1.fasta 13 ${ffile/_qual.fastq/_fil.fastq};
done  

