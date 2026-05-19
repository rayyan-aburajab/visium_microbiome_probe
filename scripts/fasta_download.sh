#!/bin/bash

while read acc; do
    efetch -db nucleotide -id "$acc" -format fasta > ../data/FASTA_Bacteria_NCBI/"$acc".fa
done < ../inputs/accessions.txt
