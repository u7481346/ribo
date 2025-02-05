import pandas as pd
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
module load minimap2 samtools pythonlib

def check_output():
    checkpoint_file = "${outputdir}/${sampleid}.ribocop"
    if os.path.exists(checkpoint_file + ".done"):
        print("Task already completed. Exiting")
        sys.exit(0)
    else:
        open(checkpoint_file + ".running", 'w').close()
    
    
def import_modules():
    """Load required modules."""
    modules = ["minimap2", "samtools"]
    for module in modules:
        subprocess.run(f"module load {module}", shell=True, check=True)


def get_frequent_rdna(sampleid):
    paf_file = f"{sampleid}.primary.paf"

    col_names = [
        "query_name", "query_length", "query_start", "query_end", 
        "strand", "target_name", "target_length", "target_start", "target_end", 
        "num_matches", "alignment_length", "mapping_quality"
    ]

    
    data = pd.read_csv(paf_file, sep='\t', names=col_names, usecols=[11, 7, 9, 6])  # Use relevant columns only
    data = data[data["alignment_length"] > data["target_length"] * 0.9]
    counts = data["target_name"].value_counts()
    eighteen = data.loc[data["target_name"].str.contains("_18S"), "target_name"].value_counts().idxmax()
    twoeight = data.loc[data["target_name"].str.contains("_28S"), "target_name"].value_counts().idxmax()
    
    return eighteen, twoeight


def extract_frequent_rdna(rdnalibfa, inputfasta, PBS_NCPUS, sampleid):
    subprocess.run(f"minimap2 -t {PBS_NCPUS} --secondary=no -o {sampleid}.primary.paf {rdnalibfa} {inputfasta}")
    
    eighteen, twoeight = get_frequent_rdna(sampleid)
    record = SeqIO.read(rdnalibfa, 'fasta')

    with open(f'{sampleid}.ssulsurdna.fa', 'w') as outfile:
        if record.id == eighteen:
            SeqIO.write(record, outfile, 'fasta')
        if record.id == twoeight:
            SeqIO.write(record, outfile, 'fasta')
    if record.id == eighteen:
        with open(f'{sampleid}.ssulsurdna.fa', 'w') as outfile:
            SeqIO.write(record, outfile, 'fasta')

    if record.id == twoeight:
        with open(f'{sampleid}.ssulsurdna.fa', 'w') as outfile:
            SeqIO.write(record, outfile, 'fasta')





def align_rdna(PBS_NCPUS, sampleid, inputfasta):
    subprocess.run(f"minimap2 -t {PBS_NCPUS} --secondary=no -o {sampleid}.refined.paf {sampleid}.ssulsurdna.fa {inputfasta}") 


def read_paf_file(file_path):
    # Redefine the column names for the first 12 columns
    column_names = [
        "Query sequence name", "Query sequence length", "Query start", "Query end", "Relative strand",
        "Target sequence name", "Target sequence length", "Target start", "Target end", "Number of residue matches",
        "Alignment block length", "Mapping quality"
    ]

    # Reading the file with more control over parsing
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            # Split the line by tab
            parts = line.strip().split('\t')
            data.append(parts)
    # Convert the list of lists to a DataFrame
    data_df = pd.DataFrame(data)
    # Assign column names to the first 12 columns
    data_df.columns = column_names + [f'Var_col_{i}' for i in range(1, len(data_df.columns) - 12 + 1)]
    # Convert relevant columns to numeric types
    numeric_columns = ["Query sequence length", "Query start", "Query end", 
                       "Target sequence length", "Target start", "Target end", "Number of residue matches",
                       "Alignment block length", "Mapping quality"]
    for col in numeric_columns:
        data_df[col] = pd.to_numeric(data_df[col])
    return data_df


def process_paf_alignments(data_df):
    
    # Filter out short alignments, retaining alignments that cover at least 95% of the target sequence
    data_df = data_df[data_df['Alignment block length'].astype(float) / data_df['Target sequence length'].astype(float) >= 0.95]
    data_df = data_df[data_df['Number of residue matches'].astype(float) / data_df['Target sequence length'].astype(float) >= 0.90]

    # Ensure 18S is followed by 28S or 28S is followed by 18S for each query sequence name
    # Current logic captures 28S-18S units. we need to change it to 18S-28S units.
    morphs = set()
    for query_name, group in data_df.groupby("Query sequence name"):
        group = group.sort_values(by=["Query start"]).reset_index(drop=True)
        valid_group = set()
        i = 0
        while i < len(group):
            # start with a 18S in the plus strand, and 100bp available before this 18S on the contig
            if '_18S' in group.at[i, 'Target sequence name'] and group.at[i, 'Relative strand'] == "+" and group.at[i, 'Query start'] - 100 > 0:
                #check if the next one is a valid 28S
                if i + 1 < len(group) and '_28S' in group.at[i + 1, 'Target sequence name'] and group.at[i + 1, 'Relative strand'] == "+":
                    # check if the next one is a valid 18S
                    if i + 2 < len(group) and '_18S' in group.at[i + 2, 'Target sequence name'] and group.at[i + 2, 'Relative strand'] == "+":
                        morph = group.at[i, 'Query sequence name'], group.at[i, 'Query start'] - 100, group.at[i + 2, 'Query start'] - 100, '+'
                        morphs.add(morph)
            # start with a 18S in the minus strand
            elif '_18S' in group.at[i, 'Target sequence name'] and group.at[i, 'Relative strand'] == "-":
                #check if the next one is a valid 28S
                if i + 1 < len(group) and '_28S' in group.at[i + 1, 'Target sequence name'] and group.at[i, 'Relative strand'] == "-":
                    # check if the next one is a valid 18S
                    if i + 2 < len(group) and '_18S' in group.at[i + 2, 'Target sequence name'] and group.at[i + 2, 'Relative strand'] == "-" and group.at[i + 2, 'Query end'] + 100 <= group.at[i, 'Query sequence length']:
                        morph = group.at[i, 'Query sequence name'], group.at[i, 'Query end'] + 100, group.at[i + 2, 'Query end'] + 100, '-'
                        morphs.add(morph)                    

            i += 1
    # Convert morphs set to DataFrame
    morphs_df = pd.DataFrame(list(morphs), columns=['Query sequence name', 'Start', 'End', 'Strand'])
    return morphs_df
    # return morphs

def extract_sequences(fasta_file, morphs_df, output_dir, sampleid):
    sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    extracted_sequences = []
    valid_rows = []

    for _, row in morphs_df.iterrows():
        seq_id = row['Query sequence name']
        start = int(row['Start'])
        end = int(row['End'])
        strand = row['Strand']

        if seq_id in sequences:
            sequence = sequences[seq_id].seq[start:end]
            if strand == '-':
                sequence = sequence.reverse_complement()
            sequence_str = str(sequence)
            if 'N' not in sequence_str:
                extracted_sequences.append((seq_id + ':' + str(start) + '-' + str(end) + ':' + strand, sequence_str))
                valid_rows.append(row)

    with open(os.path.join(output_dir, sampleid + ".rDNA.morphs.fasta"), 'w') as output_handle:
        for seq_id, sequence in extracted_sequences:
            output_handle.write(f">{seq_id}\n{sequence}\n")

    valid_morphs_df = pd.DataFrame(valid_rows)
    output_file = os.path.join(output_dir, sampleid + ".rDNA.morphs.tsv")
    valid_morphs_df.to_csv(output_file, sep='\t', index=False, header=False)

def rna_builder():
    parser = argparse.ArgumentParser(description="Process PAF alignments and extract sequences.")
    parser.add_argument('-p', '--paf', required=True, help='PAF alignments file')
    parser.add_argument('-f', '--fasta', required=True, help='FASTA file of query sequences')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory for results')
    parser.add_argument('-s', '--sampleid', required=True, help='Sample ID for output files')

    args = parser.parse_args()

    if not os.path.isfile(args.paf):
        raise FileNotFoundError(f"PAF file {args.paf} does not exist.")
    if not os.path.isfile(args.fasta):
        raise FileNotFoundError(f"FASTA file {args.fasta} does not exist.")
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    data_df = read_paf_file(args.paf)
    morphs_df = process_paf_alignments(data_df)
    extract_sequences(args.fasta, morphs_df, args.output_dir, args.sampleid)


def map_individual():
    pass

def process_paf():
    pass

def complete_script():
    pass

def main():
    check_output()
    import_modules()
    align_first()
    align_second()
    read_paf_file()
    process process_paf_alignments()
    extract_sequences()
    rna_builder()
    map_individual()
    process_paf()
    complete_script()

if __name__ == "__main__":
    main()