# K-mer tool

Calculate k-mer count tables from fastq(.gz) files.

## Preprocessing

Extract all sequences from a fastq(.gz) file and split on Ns (actually
all other characters from `A`, `C`, `G`, `T`).

```
./preproc.sh -k 4 SRR043348_1.fastq.gz
```

## Building a model

- `make` with the correct `K` (line 13; `kmer.cc`);
- NB use `release` build;
- use streaming processing: `./preproc.sh -k 4 SRR043348_1.fastq.gz | ./kmer 4 > model.bin`.

## Comparing models

All functions are available in `kmer.cc`.

```
zcat /exports/lgtc/projects/2016-07-06_exome-genome/genome/8503_R1.fastq.gz | awk '{if(++count%4==1) printf "%s|", $0; else if(count%4==2) print $1;}' | ./kmer 3 counts/genome_3.bin counts/exome_3.bin | gzip > 8503_R1_select3.fasta.gz
```

