#!/bin/bash

if [ $# -eq 0 ]
  then
    printf "\nUsage: <script> <full_path_to_blastdb> <full_path_to_query_fasta>\n\n"
    exit
fi


if [ -z "$1" ]
  then
    printf "\nAborted: No database supplied\n\n"
    exit
fi

if [ -z "$2" ]
  then
    printf "\nAborted: No query supplied\n\n"
    exit
fi

printf "\nRequirements satisfied, proceeding...\n\n"

printf "Running blast search\n\n"

blastn -db $1 -query $2 -outfmt 6 -out ASV_to_MAG.blastout -evalue 0.01

printf "Sorting blast output. Make sure to check column 4 alignment length!\n\n"

cat ASV_to_MAG.blastout | sort -nrk3,3 > ASV_to_MAG_sorted.blastout

rm ASV_to_MAG.blastout
