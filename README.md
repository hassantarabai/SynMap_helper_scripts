# SynMap_helper_scripts
Tools to facilitate Synteny analysis using SynMap

## **Inversion of nucleotide FASTA and GFF files**

### **invert_fasta.py & invert_gff.py**

- The script **invert_fasta.py:** This script inverts nucleotide sequences in a multi-contig FASTA file by generating their reverse complement.
- The **invert_gff.py:** This script prepares a GFF file with inverted coordinates. It requires an input GFF file along with the associated nucleotide FASTA file to correctly compute the inverted positions and strands.
 
Both scripts rely on [Biopython](https://biopython.org/) for functionality.

## Installation

with conda:
```
conda install -c conda-forge biopython
```
# Output
- invert_fasta.py: Outputs a new FASTA file with reverse-complemented sequences.
- invert_gff.py: Outputs a new GFF file with inverted coordinates and strand information.

# Note
You can use the **agat_convert_sp_gxf2gxf.pl** command from [agat package](https://anaconda.org/bioconda/agat) to modify the generated inverted gff files to be SynMap compatible:

# Use
```
agat_convert_sp_gxf2gxf.pl --gff inverted_gff -o Synmap_compatible_inverted_gff
```
