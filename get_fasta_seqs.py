#!/usr/bin/env python

from Bio import SeqIO
import sys

syntax = '''
------------------------------------------------------------------------------------
Syntax:        python extract_sequence_by_name_list.py *file1.fasta *file2.txt
*Sequences in fasta format 
**List of sequences to extract; must have the same name as in fasta file without '>'
------------------------------------------------------------------------------------
'''
if len(sys.argv) != 3:
        print(syntax)
        sys.exit()

from Bio import SeqIO                                                               
import sys                                                                          

wanted = [line.strip() for line in open(sys.argv[2])]                               
seqiter = SeqIO.parse(open(sys.argv[1]), 'fasta')                                    
SeqIO.write((seq for seq in seqiter if seq.id in wanted), sys.stdout, "fasta")
