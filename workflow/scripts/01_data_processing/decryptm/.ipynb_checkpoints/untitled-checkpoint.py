if 'snakemake' in locals():
   int_file = snakemake.input[0]
   ref_proteome_file = snakemake.input[1]
   output_file = snakemake.output[0]
else:
   int_file = '../../../data/decryptm/10_Kinase_Inhibitors/Phosphoproteome/curves_2KI.txt'
   ref_proteome_file = '../../../data/decryptm/uniprot_proteome_up000005640_03112020.fasta'
   output_file = '../../../results/decryptm/protein_mapping/mapped_protein_2KI.csv'

import pandas as pd
import numpy as np
import re
from Bio import SeqIO

def parse_phosphosite_info(prob_peptide):
    ast_peptide = re.sub('\(.*?\)', '*', prob_peptide)
    seq_peptide = re.sub('\(.*?\)', '', prob_peptide)
    # locate the position before every asterisk
    phospho_positions = [m.start() for m in re.finditer('\*', ast_peptide)]
    phospho_positions = [i - n for n, i in enumerate(phospho_positions)]
    # parse the aminoacid in interesting positions
    phospho_aas = [seq_peptide[i - 1] for i in phospho_positions]
    # extract all phospho probabilities
    phospho_probs = [float(i[1:-1]) for i in re.findall('\(.*?\)', prob_peptide)]
    return phospho_positions, phospho_aas, phospho_probs


# read fasta file as a dictionary
ref_proteome = SeqIO.to_dict(SeqIO.parse(ref_proteome_file, "fasta"))

# replace IDs by the second element after splitting by '|'
ref_proteome = {k.split('|')[1]: v for k, v in ref_proteome.items()}

# read tsv file
df = pd.read_csv(int_file, sep='\t')
protein_id_column = 'Leading proteins'
phospho_count_column = 'Phospho (STY)'
phospho_prob_column = 'All Phospho (STY) Probabilities'
seq_column = 'Sequence'

# remove all rows on which phospho count is 0
df = df[df[phospho_count_column] != 0].copy()
filt_df = df.iloc[:, :].copy()

# parse phosphosite info
psite_location = []
peptide_start = []
# iterate over the rows of the df
for index, row in df.iterrows():
    phospho_count = int(row[phospho_count_column])
    phospho_prob_pep = row[phospho_prob_column].split(';')[0]
    lead_prot = row[protein_id_column].split(';')[0]
    pos, aa, prob = parse_phosphosite_info(phospho_prob_pep)
    # keep only the phospho_count positions and aminoacids with the highest probability
    int_psite_index = np.argsort(prob)[-phospho_count:]
    int_pos = [pos[i] for i in int_psite_index]
    int_aa = [aa[i] for i in int_psite_index]
    # locate start position of the peptide in the protein
    # if the lead protein is not a key in the ref_proteome dictionary, skip the row and append an empty string
    if lead_prot not in ref_proteome:
        psite_location.append('')
        peptide_start.append('')
        continue
    ref_seq = ref_proteome[lead_prot].seq
    start_pos = ref_seq.find(row[seq_column]) + 1
    # add the start position to the phosphosite position
    int_pos = [i + start_pos -1 for i in int_pos]
    # concatenate aa and positions without separator
    out_psite = ';'.join([i + str(j) for i, j in zip(int_aa, int_pos)])
    psite_location.append(out_psite)
    peptide_start.append(start_pos)

# create a new column in filt df with the psite location
filt_df['psite_location'] = psite_location
filt_df['peptide_start'] = peptide_start

# write to file after keeping only the columns: Proteins	Leading proteins	Phospho (STY)	All Phospho (STY) Probabilities	psite_location
out_df = filt_df[['Proteins', 'Leading proteins', 'Phospho (STY)', 'All Phospho (STY) Probabilities', 'psite_location', 'peptide_start', 'Sequence']].copy()
out_df.to_csv(output_file, sep='\t', index=False)