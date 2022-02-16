#!/bin/bash

if ! command -v blastn &> /dev/null
then
    printf "\nAborting: blastn could not be found in your path!\n\n"
    exit
fi

if ! command -v makeblastdb &> /dev/null
then
    printf "\nAborting: makeblastdb could not be found in your path!\n\n"
    exit
fi

if ! command -v bedtools &> /dev/null
then
    printf "\nAborting: bedtools could not be found in your path!\n\n"
    exit
fi

printf "Programs found! Proceeding...\n\n"

for f in ./16S_subunit_summaries/*_16S.tsv; do cat $f | cut -f1,4,5 > $(basename $f _16S.tsv)_16S.bed; done

for f in ./*.bed; do bedtools getfasta -fi $(basename $f _16S.bed).fasta -bed $f -fo $(basename $f _16S.bed)_16S.fasta; done

cat *_16S.fasta > all_genomes_16S.fasta

makeblastdb -in all_genomes_16S.fasta -dbtype nucl -parse_seqids


