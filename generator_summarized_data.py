import os.path
from os import path

file_prefix = ".\observed\OGs by CDS"

overlaps_by_strands = ['diff_stranded', 'same_stranded', 'both_stranded']
overlaps_by_ORF = ['diff_ORF', 'same_ORF', 'both_ORF']
annotations = ['NCBI', 'Ensembl']

species_list = [
    'Human',
    'Mouse',
    'Rat',
    'Drosophila',
    'ZebraFish',
]
print(path.exists("generated_data/xxx.txt"))
file = open("generated_data/xxx.txt", 'r')
lines = file.readlines()
print(lines)
print(path.exists("generated_data/observed/OGs by CDS [diff_stranded, diff_ORF, Mouse, Ensembl].txt"))

def read_necessary_observed_data(file_path):
    # y='./observ'
    # x=".\observed\OGs by CDS [diff_stranded, diff_ORF, Mouse, Ensembl].txt"
    # x=".\observed\OGs by CDS [diff_stranded, diff_ORF, Mouse, Ensembl].txt'"
    if not exists(file_path): return
    print(file_path)
    file = open(file_path, 'r')
    lines = file.readlines()
    for i in range(len(lines)):
        if lines[i].__contains__('nucleotides:'):
            line = lines[i].split('\t')
            print(line)


for annotation in annotations:
    for species in species_list:
        required_files = []
        for strands_type in overlaps_by_strands:
            for ORF_type in overlaps_by_ORF:
                file_path = f'{file_prefix} [{strands_type}, {ORF_type}, {species}, {annotation}].txt'
                read_necessary_observed_data(file_path)
