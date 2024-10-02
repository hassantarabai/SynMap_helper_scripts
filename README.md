# SynMap_helper_scripts
Tools to facilitate Synteny analysis using SynMap

## **Inversion of nucleotide fasta and gff files**
### **invert_fasta.py & invert_gff.py **
The script * invert_fasta.py * is used to invert nueclotide sequences and takes a multi-contig fasta file as an input.
The * invert_gff.py * script can be used to preapre a gff file with the inverted coordinates of an input gff file. The script takes as input a gff file and its associated nucleotide fasta file. 

The scripts require [Biopython](https://biopython.org/) for thier function function.

## Installation

with conda:
```
conda install -c conda-forge biopython
```
