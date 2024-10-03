# SynMap_helper_scripts
Tools to facilitate Synteny analysis using SynMap

## **Inversion of nucleotide FASTA and GFF files**

### **invert_fasta.py & invert_gff.py**

- **invert_fasta.py:** This script inverts nucleotide sequences in a multi-contig FASTA file by generating their reverse complement.
- **invert_gff.py:** This script prepares a GFF file with inverted coordinates. It requires an input GFF3 file along with the associated nucleotide FASTA file to correctly compute the inverted positions and strands.

### Output
- invert_fasta.py: Outputs a new FASTA file with reverse-complemented sequences.
- invert_gff.py: Outputs a new GFF3 file with inverted coordinates and strand information.

### Note
You can use the **agat_convert_sp_gxf2gxf.pl** command from [agat package](https://anaconda.org/bioconda/agat) to modify the generated inverted gff files for SynMap compatibility:

#### Use
```
agat_convert_sp_gxf2gxf.pl --gff inverted_gff -o Synmap_compatible_inverted_gff
```

### **get_CDS_AA.py**

This script extract cds based on an input of GFF3 and corresponding nucelotide FASTA file, then translate the cds to amino acids (AA).

### Output
-Nucelotide FASTA file representing identified cds.
-AA FASTA file representing the translated cds.

 
## The scripts rely on [Biopython](https://biopython.org/) for functionality.

### Installation

with conda:
```
conda install -c conda-forge biopython
```
