#!/bin/bash

for i in `seq 2 12`;
do
    echo "/home/jkvis/kmer/preproc.sh -k $i /exports/lgtc/projects/2016-07-06_exome-genome/genome/8503_R1.fastq.gz | /home/jkvis/kmer/kmer $i > /home/jkvis/kmer/8503_R1_$i.bin" | qsub -pe make 3 -l h_vmem=1G -l hostname="catshark|chimerashark|elephantshark|iridescentshark|leopardshark|sleepershark|swellshark"
    echo "/home/jkvis/kmer/preproc.sh -k $i /exports/lgtc/projects/2016-07-06_exome-genome/genome/8503_R2.fastq.gz | /home/jkvis/kmer/kmer $i > /home/jkvis/kmer/8503_R2_$i.bin" | qsub -pe make 3 -l h_vmem=1G -l hostname="catshark|chimerashark|elephantshark|iridescentshark|leopardshark|sleepershark|swellshark"
    echo "/home/jkvis/kmer/preproc.sh -k $i /exports/lgtc/projects/2016-07-06_exome-genome/exome/PGLHGT1_1.fastq.gz | /home/jkvis/kmer/kmer $i > /home/jkvis/kmer/PGLHGT1_1_$i.bin" | qsub -pe make 3 -l h_vmem=1G -l hostname="catshark|chimerashark|elephantshark|iridescentshark|leopardshark|sleepershark|swellshark"
    echo "/home/jkvis/kmer/preproc.sh -k $i /exports/lgtc/projects/2016-07-06_exome-genome/exome/PGLHGT1_2.fastq.gz | /home/jkvis/kmer/kmer $i > /home/jkvis/kmer/PGLHGT1_2_$i.bin" | qsub -pe make 3 -l h_vmem=1G -l hostname="catshark|chimerashark|elephantshark|iridescentshark|leopardshark|sleepershark|swellshark"
done

