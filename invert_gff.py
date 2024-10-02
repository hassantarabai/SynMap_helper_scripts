import os
from Bio import SeqIO

# Function to calculate contig lengths from the FASTA file
def get_contig_lengths(fasta_file):
    contig_lengths = {}
    with open(fasta_file, "r") as fasta_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            contig_lengths[record.id] = len(record.seq)
    return contig_lengths

# Function to invert GFF3 coordinates and strand
def invert_gff(input_gff, contig_lengths, output_gff):
    with open(input_gff, "r") as infile, open(output_gff, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                # Keep comments and headers intact
                outfile.write(line)
                continue

            parts = line.strip().split("\t")
            contig = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]

            if contig in contig_lengths:
                contig_len = contig_lengths[contig]
                # Invert the coordinates
                new_start = contig_len - end + 1
                new_end = contig_len - start + 1
                parts[3] = str(new_start)
                parts[4] = str(new_end)

                # Flip the strand
                if strand == "+":
                    parts[6] = "-"
                elif strand == "-":
                    parts[6] = "+"

                # Write the updated line to the output GFF file
                outfile.write("\t".join(parts) + "\n")
            else:
                print(f"Contig {contig} not found in FASTA file, skipping line.")

# Ask for the input folder containing the GFF files
input_folder = input("Enter the path to the folder containing the .gff files: ")

# Process all GFF files in the input folder
gff_files = [f for f in os.listdir(input_folder) if f.endswith(".gff") or f.endswith(".gff3")]

for gff_file in gff_files:
    # Ask for the corresponding FASTA file for each GFF
    fasta_file = input(f"Enter the corresponding FASTA file for {gff_file}: ")
    
    # Calculate contig lengths from the provided FASTA file
    contig_lengths = get_contig_lengths(fasta_file)
    
    # Define input and output paths for GFF
    input_file = os.path.join(input_folder, gff_file)
    output_file = os.path.join(input_folder, f"inverted_{gff_file}")
    
    # Perform inversion of coordinates
    invert_gff(input_file, contig_lengths, output_file)
    print(f"Processed: {gff_file} -> {output_file}")

print("All GFF files have been processed and inverted.")
