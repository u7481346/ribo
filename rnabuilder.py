import pandas as pd
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

#script copied from kango2/ausarg

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

def main():
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

if __name__ == "__main__":
    main()
# there is a corner case as follows where it seems like contigs are butterflied. this need to be handled either in morph cadidate generation or once candidates are generated, remove both if they overlap by >1bp

# h2tg000001l     293025  98998   100751  +       AY217912.1.1767_18S     1767    5       1758    1689    1753    60      tp:A:P  cm:i:309        s1:i:1689       s2:i:0  dv:f:0  rl:i:0
# h2tg000001l     293025  101927  104750  +       HACI01071117.32.3921_28S        3890    22      2758    2009    2826    60      tp:A:P  cm:i:318        s1:i:1990       s2:i:0  dv:f:0.0004     rl:i:0
# h2tg000001l     293025  104751  107574  -       HACI01071117.32.3921_28S        3890    22      2758    2025    2826    60      tp:A:P  cm:i:323        s1:i:2006       s2:i:0  dv:f:0.0004     rl:i:0
# h2tg000001l     293025  108757  110510  -       AY217912.1.1767_18S     1767    5       1758    1689    1753    60      tp:A:P  cm:i:309        s1:i:1689       s2:i:0  dv:f:0  rl:i:0
