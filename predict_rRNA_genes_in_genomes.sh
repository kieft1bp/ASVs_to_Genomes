#!/bin/bash

## NOTE: All genome nucleotide fasta files should be in the cwd

if ! command -v blastn &> /dev/null
then
    printf "\nAborting: blastn could not be found in your path!\n\n"
    exit
fi

if test -f "get_fasta_seqs.py"; then
    printf "\nFound scripts, proceeding...\n"
else
    printf "\nAborted: You need to download the get_fasta_seqs.py script from the github repo\n\n"
    exit
fi

if test -d "5-16-23S_subunit_summaries"; then
    printf "\nRemove or rename old 5-16-23S_subunit_summaries/ directory, it already exists\n"
    printf "\nE.g., you can remove with <rm -rf 5-16-23S_subunit_summaries/>\n\n"
    exit
fi

if test -d "16S_subunit_summaries"; then
    printf "\nRemove or rename old 16S_subunit_summaries/ directory, it already exists\n"
    printf "\nE.g., you can remove with <rm -rf 16S_subunit_summaries/>\n\n"
    exit
fi

if test -f "./blastdb/16S_ribosomal_RNA.bti";then
    printf "\n16S blastn database found!\n"
else
    printf "\nAborted: You need to download the blastdb/ directory from the github repo\n\n"
    exit
fi

chmod +x ./get_fasta_seqs.py

printf "\nPart 1 of 3: Predicting ribosomal RNA genes\n\n"

touch removed_rRNA_summary.tsv ## make this to subvert the error for the rm command below

mkdir 5-16-23S_subunit_summaries
mkdir 16S_subunit_summaries

rm *rRNA_summary.tsv ## remove existing summary .tsv files since we're writing iteratively to later files

for f in ./*.fasta; do echo "Predicting ribosomal RNAs from" $(basename $f); echo $(basename $f .fasta) >> ./5-16-23S_subunit_summaries/$(basename $f .fasta)_rRNA_summary.tsv; barrnap --quiet $f >> ./5-16-23S_subunit_summaries/$(basename $f .fasta)_rRNA_summary.tsv; done

for f in ./5-16-23S_subunit_summaries/*rRNA_summary.tsv; do cat $f | grep "16S_" | grep -v "#" > ./16S_subunit_summaries/$(basename $f _rRNA_summary.tsv)_16S.tsv; done

printf "\nPart 2 of 3: Creating summary files of sequence and taxonomy from 16S subunits\n"

for f in ./16S_subunit_summaries/*16S.tsv; do cat $f | cut -f1 > list; ./get_fasta_seqs.py ./$(basename $f _16S.tsv).fasta list > seqs; blastn -query seqs -evalue 0.01 -db ./blastdb/16S_ribosomal_RNA -max_target_seqs 1 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames" 2>/dev/null > blastn_results.tsv; cat seqs | tr '\n' ' ' | sed 's/>/\n/g' | cut -d" " -f2- | sed 's/ //g' | sed '/^$/d' > seqs_clean; paste $f seqs_clean blastn_results.tsv > ./16S_subunit_summaries/$(basename $f .tsv)_clean.tsv; rm list seqs*; done

for f in ./16S_subunit_summaries/*16S_clean.tsv; do cat $f | cut -f1,13,17-20,23 | sort -nrk2,2 | sed '1iquery\tpercent_id\tqstart\tqend\tsstart\tsend\ttaxon' > ./16S_subunit_summaries/$(basename $f .tsv)_summarized.tsv; done

printf "\nPart 3 of 3: Counting the instances of 5/16/23S subunits (for bin "purity" estimate)\n"

touch rRNA_subunit_summary.tsv

rm rRNA_subunit_summary.tsv; for f in ./5-16-23S_subunit_summaries/*RNA_summary.tsv; do samplename=`echo $(basename $f _rRNA_summary.tsv)`; a16S_count=`cat $f | grep "16S" | wc -l`; a5S_count=`cat $f | grep "5S" | wc -l`; a23S_count=`cat $f | grep "23S" | wc -l`; cat $f | tr '\n' '\t' | grep "16S" | grep "5S" | grep "23S" | sed '1iheader' | awk -v var="$samplename" -v var1="$a16S_count" -v var2="$a5S_count" -v var3="$a23S_count" ' END { if (NR > 1) print var"\tyes\t"var1"\t"var2"\t"var3; if (NR < 1) print var"\tno\t"var1"\t"var2"\t"var3;}' >> rRNA_subunit_summary.tsv; done; sed -i '1iGenome\thas_three_subunits\t16S_count\t5S_count\t23S_count' rRNA_subunit_summary.tsv

printf "\nGood luck sifting through the mapping info!\n\n"
