import os

os.system("python ./analyzer_CDS_control_stats.py Homo_sapiens Ensembl")
os.system("python ./analyzer_OGs_by_CDS.py Homo_sapiens Ensembl diff_stranded both_ORF")
os.system("python ./analyzer_CDS_control_stats.py Homo_sapiens NCBI")
os.system("python ./analyzer_OGs_by_CDS.py Homo_sapiens NCBI diff_stranded both_ORF")

os.system("python ./analyzer_CDS_control_stats.py Mus_musculus Ensembl")
os.system("python ./analyzer_OGs_by_CDS.py Mus_musculus Ensembl diff_stranded both_ORF")
os.system("python ./analyzer_CDS_control_stats.py Mus_musculus NCBI")
os.system("python ./analyzer_OGs_by_CDS.py Mus_musculus NCBI diff_stranded both_ORF")

os.system("python ./analyzer_CDS_control_stats.py Rattus_norvegicus Ensembl")
os.system("python ./analyzer_OGs_by_CDS.py Rattus_norvegicus Ensembl diff_stranded both_ORF")
os.system("python ./analyzer_CDS_control_stats.py Rattus_norvegicus NCBI")
os.system("python ./analyzer_OGs_by_CDS.py Rattus_norvegicus NCBI diff_stranded both_ORF")

os.system("python ./analyzer_CDS_control_stats.py Danio_rerio Ensembl")
os.system("python ./analyzer_OGs_by_CDS.py Danio_rerio Ensembl diff_stranded both_ORF")
os.system("python ./analyzer_CDS_control_stats.py Danio_rerio NCBI")
os.system("python ./analyzer_OGs_by_CDS.py Danio_rerio NCBI diff_stranded both_ORF")

os.system("python ./analyzer_CDS_control_stats.py Drosophila_melanogaster Ensembl")
os.system("python ./analyzer_OGs_by_CDS.py Drosophila_melanogaster Ensembl diff_stranded both_ORF")
os.system("python ./analyzer_CDS_control_stats.py Drosophila_melanogaster NCBI")
os.system("python ./analyzer_OGs_by_CDS.py Drosophila_melanogaster NCBI diff_stranded both_ORF")


os.system("python generator_gene_outline.py Drosophila_melanogaster Ensembl 50 1000")
os.system("python generator_gene_outline.py Danio_rerio Ensembl 50 1000")
os.system("python generator_gene_outline.py Rattus_norvegicus Ensembl 50 1000")
os.system("python generator_gene_outline.py Mus_musculus Ensembl 50 1000")
os.system("python generator_gene_outline.py Homo_sapiens Ensembl 50 1000")

os.system("python generator_gene_outline.py Drosophila_melanogaster NCBI 50 1000")
os.system("python generator_gene_outline.py Danio_rerio NCBI 50 1000")
os.system("python generator_gene_outline.py Rattus_norvegicus NCBI 50 1000")
os.system("python generator_gene_outline.py Mus_musculus NCBI 50 1000")
os.system("python generator_gene_outline.py Homo_sapiens NCBI 50 1000")

os.system("python generator_gene_outline_comparative.py Ensembl 50 1000")
os.system("python generator_gene_outline_comparative.py NCBI 50 1000")

