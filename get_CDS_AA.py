from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation

# Prompt user for GFF and genome FASTA file paths
gff_file = input("Please enter the path to the GFF file: ")
genome_fasta = input("Please enter the path to the genome FASTA file: ")

# Load genome sequence
genome_dict = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

# Print scaffold names in the genome FASTA file
print("Scaffold names in the genome FASTA:")
for scaffold in genome_dict.keys():
    print(scaffold)

# Lists to store CDS sequences and proteins
cds_sequences = []
protein_sequences = []

# Parsing the GFF file manually
with open(gff_file) as gff_handle:
    for line in gff_handle:
        if line.startswith("#"):
            continue  # Skip header lines
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue  # Skip malformed lines

        seq_id = parts[0]
        feature_type = parts[2]
        start = int(parts[3])
        end = int(parts[4])
        strand = parts[6]
        attributes = parts[8]

        # Extract the CDS if the feature type is CDS
        if feature_type == "CDS":
            # Parse the ID from the attributes (assuming "ID=" exists)
            cds_id = "unknown"
            for attribute in attributes.split(";"):
                if attribute.startswith("ID="):
                    cds_id = attribute.split("=")[1]

            # Extract the relevant sequence from the genome
            if seq_id in genome_dict:
                feature_seq = genome_dict[seq_id].seq[start-1:end]  # Adjust to 0-based index
                if strand == "-":
                    feature_seq = feature_seq.reverse_complement()

                # Store the CDS sequence
                cds_record = SeqRecord(feature_seq, id=cds_id, description="CDS sequence")
                cds_sequences.append(cds_record)

                # Translate the CDS to a protein sequence
                protein_seq = feature_seq.translate(to_stop=True)
                protein_record = SeqRecord(protein_seq, id=cds_id, description="Protein sequence")
                protein_sequences.append(protein_record)

# Output file names
cds_output = input("Please enter the output file name for CDS sequences (e.g., cds_sequences.fasta): ")
protein_output = input("Please enter the output file name for protein sequences (e.g., protein_sequences.fasta): ")

# Save the CDS sequences to a FASTA file
with open(cds_output, "w") as cds_outfile:
    SeqIO.write(cds_sequences, cds_outfile, "fasta")

# Save the protein sequences to a FASTA file
with open(protein_output, "w") as protein_outfile:
    SeqIO.write(protein_sequences, protein_outfile, "fasta")

print("CDS and protein sequences have been extracted and saved.")
