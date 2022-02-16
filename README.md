# ASV to Genome mapping

1. Create ASVs from your paired-end .fastq files. This is accomplished with the dada2_workflow.R script.
2. Identify rRNA genes (and fragments) inside your genomes (e.g., MAGs, SAGs). This is accomplished with the predict_rRNA_genes_in_genomes.sh script.
3. Create a blastn database from identified 16S rRNA genes in your genomes. This is accomplished with the create_16S_blastn_ref_database.sh script.
4. Run homology search of ASV 16S gene fragments against your genome database. This is accomplished with the blast_ASV_16S_against_Genome_16S.sh script.
5. Interpret your results! You'll end up with a blast tabular file (outfmt 6) describing ASV sequence identity and alignment length to a genome 16S gene. Depending on your amplicon length, I would only consider ASV-Genome "pairs" when there is 0 or 1 polymorphisms between sequences. For example, if you amplify the V4 region (~250bp), you should have a sequence identity cutoff of ~99.5%, allowing for 1 nucleotide difference between the ASV fragment and the genome gene. Ideally, your entire ASV amplicon (in the example previous, all 250 bp) should align. 
