import os

# os.system("python overlapping_fragments.py Drosophila_melanogaster NCBI diff_stranded both_ORF print_gene_pairs_with_descriptions")
# os.system("python overlapping_fragments.py Danio_rerio NCBI diff_stranded both_ORF print_gene_pairs_with_descriptions")
# os.system("python overlapping_fragments.py Mus_musculus NCBI diff_stranded both_ORF print_gene_pairs_with_descriptions")
# os.system("python overlapping_fragments.py Rattus_norvegicus NCBI diff_stranded both_ORF print_gene_pairs_with_descriptions")
# os.system("python overlapping_fragments.py Homo_sapiens NCBI diff_stranded both_ORF print_gene_pairs_with_descriptions")

# os.system("analyze_og.py Drosophila_melanogaster NCBI")
# os.system("analyze_og.py Caenorhabditis_elegans NCBI")
# os.system("analyze_og.py Danio_rerio NCBI")
# os.system("python analyze_og.py Mus_musculus NCBI")

#
# os.system("python base_occurrences.py Drosophila_melanogaster Ensembl 50 1000")
# os.system("python base_occurrences.py Danio_rerio Ensembl 50 1000")
# os.system("python base_occurrences.py Mus_musculus Ensembl 50 1000")
# os.system("python base_occurrences.py Rattus_norvegicus Ensembl 50 1000")
# os.system("python base_occurrences.py Homo_sapiens Ensembl 50 1000")


os.system("python base_occurrences.py Drosophila_melanogaster NCBI 50 1000")
os.system("python base_occurrences.py Danio_rerio NCBI 50 1000")
os.system("python base_occurrences.py Mus_musculus NCBI 50 1000")
os.system("python base_occurrences.py Rattus_norvegicus NCBI 50 1000")
os.system("python base_occurrences.py Homo_sapiens NCBI 50 1000")

os.system("python base_occurrences_comparative.py NCBI 50 1000")
os.system("python base_occurrences_comparative.py Ensembl 50 1000")
#
