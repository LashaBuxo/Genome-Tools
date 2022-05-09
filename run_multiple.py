import os

#os.system("python ./generator_summarized_data.py")

# os.system("python ./generator_control_data.py Ensembl")
# os.system("python ./generator_control_data.py NCBI")

os.system("python ./generator_observed_data.py Homo_sapiens NCBI diff_stranded same_ORF")
os.system("python ./generator_observed_data.py Homo_sapiens NCBI diff_stranded diff_ORF")
os.system("python ./generator_observed_data.py Homo_sapiens NCBI diff_stranded both_ORF")

os.system("python ./generator_observed_data.py Homo_sapiens NCBI same_stranded diff_ORF")
os.system("python ./generator_observed_data.py Homo_sapiens NCBI same_stranded same_ORF")
os.system("python ./generator_observed_data.py Homo_sapiens NCBI same_stranded both_ORF")

os.system("python ./generator_observed_data.py Mus_musculus NCBI diff_stranded diff_ORF")
os.system("python ./generator_observed_data.py Mus_musculus NCBI diff_stranded same_ORF")
os.system("python ./generator_observed_data.py Mus_musculus NCBI diff_stranded both_ORF")

os.system("python ./generator_observed_data.py Mus_musculus NCBI same_stranded diff_ORF")
os.system("python ./generator_observed_data.py Mus_musculus NCBI same_stranded same_ORF")
os.system("python ./generator_observed_data.py Mus_musculus NCBI same_stranded both_ORF")


#os.system("python ./generator_observed_data.py Homo_sapiens Ensembl diff_stranded same_ORF")
# os.system("python ./generator_observed_data.py Homo_sapiens Ensembl diff_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Homo_sapiens Ensembl diff_stranded both_ORF")
#
# os.system("python ./generator_observed_data.py Homo_sapiens Ensembl same_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Homo_sapiens Ensembl same_stranded same_ORF")
# os.system("python ./generator_observed_data.py Homo_sapiens Ensembl same_stranded both_ORF")
#
# os.system("python ./generator_observed_data.py Mus_musculus Ensembl diff_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Mus_musculus Ensembl diff_stranded same_ORF")
# os.system("python ./generator_observed_data.py Mus_musculus Ensembl diff_stranded both_ORF")
#
# os.system("python ./generator_observed_data.py Mus_musculus Ensembl same_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Mus_musculus Ensembl same_stranded same_ORF")
# os.system("python ./generator_observed_data.py Mus_musculus Ensembl same_stranded both_ORF")
#
# os.system("python ./generator_observed_data.py Homo_sapiens NCBI diff_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Homo_sapiens NCBI diff_stranded same_ORF")
# os.system("python ./generator_observed_data.py Homo_sapiens NCBI same_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Homo_sapiens NCBI same_stranded same_ORF")
# os.system("python ./generator_observed_data.py Mus_musculus NCBI diff_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Mus_musculus NCBI diff_stranded same_ORF")
# os.system("python ./generator_observed_data.py Mus_musculus NCBI same_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Mus_musculus NCBI same_stranded same_ORF")

#
# os.system("python ./generator_observed_data.py Drosophila_melanogaster Ensembl diff_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Drosophila_melanogaster Ensembl diff_stranded same_ORF")
# os.system("python ./generator_observed_data.py Drosophila_melanogaster Ensembl diff_stranded both_ORF")
# os.system("python ./generator_observed_data.py Drosophila_melanogaster Ensembl same_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Drosophila_melanogaster Ensembl same_stranded same_ORF")
# os.system("python ./generator_observed_data.py Drosophila_melanogaster Ensembl same_stranded both_ORF")
# os.system("python ./generator_observed_data.py Drosophila_melanogaster Ensembl both_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Drosophila_melanogaster Ensembl both_stranded same_ORF")
# os.system("python ./generator_observed_data.py Drosophila_melanogaster Ensembl both_stranded both_ORF")
#
# os.system("python ./generator_observed_data.py Drosophila_melanogaster NCBI diff_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Drosophila_melanogaster NCBI diff_stranded same_ORF")
# os.system("python ./generator_observed_data.py Drosophila_melanogaster NCBI diff_stranded both_ORF")
# os.system("python ./generator_observed_data.py Drosophila_melanogaster NCBI same_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Drosophila_melanogaster NCBI same_stranded same_ORF")
# os.system("python ./generator_observed_data.py Drosophila_melanogaster NCBI same_stranded both_ORF")
# os.system("python ./generator_observed_data.py Drosophila_melanogaster NCBI both_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Drosophila_melanogaster NCBI both_stranded same_ORF")
# os.system("python ./generator_observed_data.py Drosophila_melanogaster NCBI both_stranded both_ORF")
#
# os.system("python ./generator_observed_data.py Danio_rerio Ensembl diff_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Danio_rerio Ensembl diff_stranded same_ORF")
# os.system("python ./generator_observed_data.py Danio_rerio Ensembl diff_stranded both_ORF")
# os.system("python ./generator_observed_data.py Danio_rerio Ensembl same_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Danio_rerio Ensembl same_stranded same_ORF")
# os.system("python ./generator_observed_data.py Danio_rerio Ensembl same_stranded both_ORF")
# os.system("python ./generator_observed_data.py Danio_rerio Ensembl both_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Danio_rerio Ensembl both_stranded same_ORF")
# os.system("python ./generator_observed_data.py Danio_rerio Ensembl both_stranded both_ORF")
#
# os.system("python ./generator_observed_data.py Danio_rerio NCBI diff_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Danio_rerio NCBI diff_stranded same_ORF")
# os.system("python ./generator_observed_data.py Danio_rerio NCBI diff_stranded both_ORF")
# os.system("python ./generator_observed_data.py Danio_rerio NCBI same_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Danio_rerio NCBI same_stranded same_ORF")
# os.system("python ./generator_observed_data.py Danio_rerio NCBI same_stranded both_ORF")
# os.system("python ./generator_observed_data.py Danio_rerio NCBI both_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Danio_rerio NCBI both_stranded same_ORF")
# os.system("python ./generator_observed_data.py Danio_rerio NCBI both_stranded both_ORF")
#
#
#
# os.system("python ./generator_observed_data.py Homo_sapiens Ensembl diff_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Homo_sapiens Ensembl diff_stranded same_ORF")
# os.system("python ./generator_observed_data.py Homo_sapiens Ensembl diff_stranded both_ORF")
# os.system("python ./generator_observed_data.py Homo_sapiens Ensembl same_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Homo_sapiens Ensembl same_stranded same_ORF")
# os.system("python ./generator_observed_data.py Homo_sapiens Ensembl same_stranded both_ORF")
# os.system("python ./generator_observed_data.py Homo_sapiens Ensembl both_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Homo_sapiens Ensembl both_stranded same_ORF")
# os.system("python ./generator_observed_data.py Homo_sapiens Ensembl both_stranded both_ORF")
#
#
# os.system("python ./generator_observed_data.py Homo_sapiens NCBI diff_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Homo_sapiens NCBI diff_stranded same_ORF")
# os.system("python ./generator_observed_data.py Homo_sapiens NCBI diff_stranded both_ORF")
# os.system("python ./generator_observed_data.py Homo_sapiens NCBI same_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Homo_sapiens NCBI same_stranded same_ORF")
# os.system("python ./generator_observed_data.py Homo_sapiens NCBI same_stranded both_ORF")
# os.system("python ./generator_observed_data.py Homo_sapiens NCBI both_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Homo_sapiens NCBI both_stranded same_ORF")
# os.system("python ./generator_observed_data.py Homo_sapiens NCBI both_stranded both_ORF")
#
#
#
#
# os.system("python ./generator_observed_data.py Mus_musculus Ensembl diff_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Mus_musculus Ensembl diff_stranded same_ORF")
# os.system("python ./generator_observed_data.py Mus_musculus Ensembl diff_stranded both_ORF")
# os.system("python ./generator_observed_data.py Mus_musculus Ensembl same_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Mus_musculus Ensembl same_stranded same_ORF")
# os.system("python ./generator_observed_data.py Mus_musculus Ensembl same_stranded both_ORF")
# os.system("python ./generator_observed_data.py Mus_musculus Ensembl both_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Mus_musculus Ensembl both_stranded same_ORF")
# os.system("python ./generator_observed_data.py Mus_musculus Ensembl both_stranded both_ORF")
#
# os.system("python ./generator_observed_data.py Mus_musculus NCBI diff_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Mus_musculus NCBI diff_stranded same_ORF")
# os.system("python ./generator_observed_data.py Mus_musculus NCBI diff_stranded both_ORF")
# os.system("python ./generator_observed_data.py Mus_musculus NCBI same_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Mus_musculus NCBI same_stranded same_ORF")
# os.system("python ./generator_observed_data.py Mus_musculus NCBI same_stranded both_ORF")
# os.system("python ./generator_observed_data.py Mus_musculus NCBI both_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Mus_musculus NCBI both_stranded same_ORF")
# os.system("python ./generator_observed_data.py Mus_musculus NCBI both_stranded both_ORF")
#
# os.system("python ./generator_observed_data.py Rattus_norvegicus Ensembl diff_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Rattus_norvegicus Ensembl diff_stranded same_ORF")
# os.system("python ./generator_observed_data.py Rattus_norvegicus Ensembl diff_stranded both_ORF")
# os.system("python ./generator_observed_data.py Rattus_norvegicus Ensembl same_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Rattus_norvegicus Ensembl same_stranded same_ORF")
# os.system("python ./generator_observed_data.py Rattus_norvegicus Ensembl same_stranded both_ORF")
# os.system("python ./generator_observed_data.py Rattus_norvegicus Ensembl both_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Rattus_norvegicus Ensembl both_stranded same_ORF")
# os.system("python ./generator_observed_data.py Rattus_norvegicus Ensembl both_stranded both_ORF")
#
# os.system("python ./generator_observed_data.py Rattus_norvegicus NCBI diff_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Rattus_norvegicus NCBI diff_stranded same_ORF")
# os.system("python ./generator_observed_data.py Rattus_norvegicus NCBI diff_stranded both_ORF")
# os.system("python ./generator_observed_data.py Rattus_norvegicus NCBI same_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Rattus_norvegicus NCBI same_stranded same_ORF")
# os.system("python ./generator_observed_data.py Rattus_norvegicus NCBI same_stranded both_ORF")
# os.system("python ./generator_observed_data.py Rattus_norvegicus NCBI both_stranded diff_ORF")
# os.system("python ./generator_observed_data.py Rattus_norvegicus NCBI both_stranded same_ORF")
# os.system("python ./generator_observed_data.py Rattus_norvegicus NCBI both_stranded both_ORF")

