import os

os.system("python ./generator_observed_data.py Homo_sapiens Ensembl diff_stranded diff_ORF")
os.system("python ./generator_observed_data.py Homo_sapiens Ensembl diff_stranded same_ORF")
os.system("python ./generator_observed_data.py Homo_sapiens Ensembl diff_stranded both_ORF")
os.system("python ./generator_observed_data.py Homo_sapiens Ensembl same_stranded diff_ORF")
os.system("python ./generator_observed_data.py Homo_sapiens Ensembl same_stranded same_ORF")
os.system("python ./generator_observed_data.py Homo_sapiens Ensembl same_stranded both_ORF")
os.system("python ./generator_observed_data.py Homo_sapiens Ensembl both_stranded diff_ORF")
os.system("python ./generator_observed_data.py Homo_sapiens Ensembl both_stranded same_ORF")
os.system("python ./generator_observed_data.py Homo_sapiens Ensembl both_stranded both_ORF")

#
# os.system("python ./generator_control_data.py Homo_sapiens Ensembl")
# os.system("python ./generator_observed_data.py Homo_sapiens Ensembl diff_stranded both_ORF")
# os.system("python ./generator_control_data.py Homo_sapiens NCBI")
# os.system("python ./generator_observed_data.py Homo_sapiens NCBI diff_stranded both_ORF")
#
# os.system("python ./generator_control_data.py Mus_musculus Ensembl")
# os.system("python ./generator_observed_data.py Mus_musculus Ensembl diff_stranded both_ORF")
# os.system("python ./generator_control_data.py Mus_musculus NCBI")
# os.system("python ./generator_observed_data.py Mus_musculus NCBI diff_stranded both_ORF")
#
# os.system("python ./generator_control_data.py Rattus_norvegicus Ensembl")
# os.system("python ./generator_observed_data.py Rattus_norvegicus Ensembl diff_stranded both_ORF")
# os.system("python ./generator_control_data.py Rattus_norvegicus NCBI")
# os.system("python ./generator_observed_data.py Rattus_norvegicus NCBI diff_stranded both_ORF")
#
# os.system("python ./generator_control_data.py Danio_rerio Ensembl")
# os.system("python ./generator_observed_data.py Danio_rerio Ensembl diff_stranded both_ORF")
# os.system("python ./generator_control_data.py Danio_rerio NCBI")
# os.system("python ./generator_observed_data.py Danio_rerio NCBI diff_stranded both_ORF")
#
# os.system("python ./generator_control_data.py Drosophila_melanogaster Ensembl")
# os.system("python ./generator_observed_data.py Drosophila_melanogaster Ensembl diff_stranded both_ORF")
# os.system("python ./generator_control_data.py Drosophila_melanogaster NCBI")
# os.system("python ./generator_observed_data.py Drosophila_melanogaster NCBI diff_stranded both_ORF")
#
#
# os.system("python generator_gene_outline.py Drosophila_melanogaster Ensembl 50 1000")
# os.system("python generator_gene_outline.py Danio_rerio Ensembl 50 1000")
# os.system("python generator_gene_outline.py Rattus_norvegicus Ensembl 50 1000")
# os.system("python generator_gene_outline.py Mus_musculus Ensembl 50 1000")
# os.system("python generator_gene_outline.py Homo_sapiens Ensembl 50 1000")
#
# os.system("python generator_gene_outline.py Drosophila_melanogaster NCBI 50 1000")
# os.system("python generator_gene_outline.py Danio_rerio NCBI 50 1000")
# os.system("python generator_gene_outline.py Rattus_norvegicus NCBI 50 1000")
# os.system("python generator_gene_outline.py Mus_musculus NCBI 50 1000")
# os.system("python generator_gene_outline.py Homo_sapiens NCBI 50 1000")
#
# os.system("python generator_gene_outline_comparative.py Ensembl 50 1000")
# os.system("python generator_gene_outline_comparative.py NCBI 50 1000")
#
