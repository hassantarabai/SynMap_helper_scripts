import os
from Bio import SeqIO

# Function to reverse complement FASTA files
def reverse_complement_fasta(input_file, output_file):
    with open(output_file, "w") as output_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            record.seq = record.seq.reverse_complement()
            SeqIO.write(record, output_handle, "fasta")

# Ask for the input folder containing the FASTA files
input_folder = input("Enter the path to the folder containing the .fasta files: ")

# Get a list of all .fasta files in the input folder
fasta_files = [f for f in os.listdir(input_folder) if f.endswith(".fasta")]

# Process each .fasta file
for fasta_file in fasta_files:
    input_file = os.path.join(input_folder, fasta_file)
    output_file = os.path.join(input_folder, f"inverted_{fasta_file}")
    
    # Perform reverse complement and save to a new file
    reverse_complement_fasta(input_file, output_file)
    print(f"Processed: {fasta_file} -> {output_file}")

print("All .fasta files have been processed and inverted.")
