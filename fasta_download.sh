#!/bin/bash

while read acc; do
    echo "Downloading $acc..."
    efetch -db nucleotide -id "$acc" -format fasta > ../FASTA_NCBI_Ref2/"$acc".fa
done < accessions.txt