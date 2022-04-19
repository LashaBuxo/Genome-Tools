import genome_tools as genome
import plotly.express as px
import pandas as pd
import time

start_time = time.time()

total_freq_by_regions = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]

chromosomes_dict = [None] * 6
nucleotide_dict = [None] * 6
counts_dict = [None] * 6

# all chromosomes except mitochondria
for chr_id in range(1, genome.number_of_chromosomes):
    genome.preprocess_annotation_for_chr(chr_id)
    genes_cnt = genome.genes_count_on_chr(chr_id)

    print("--- %s seconds ---" % (time.time() - start_time))
    total_freq_by_regions_on_chr = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]

    total_length = 0
    for i in range(0, genes_cnt):
        gene = genome.gene_by_ind(chr_id, i)
        total_length += gene.end - gene.start + 1
        freq_by_regions_for_gene = genome.analyze_gene_occurrences(chr_id, gene)
        for sub_reg in range(0, 6):
            for nucleotide in range(0, 4):
                total_freq_by_regions_on_chr[sub_reg][nucleotide] += freq_by_regions_for_gene[sub_reg][nucleotide]

    print(
        'UTR5_Procession=' + str(total_freq_by_regions_on_chr[0]) + '; UTR5=' + str(
            total_freq_by_regions_on_chr[1]) + '; CDS=' + str(
            total_freq_by_regions_on_chr[2]) + '; introns=' +
        str(total_freq_by_regions_on_chr[3]) + '; UTR3=' + str(
            total_freq_by_regions_on_chr[4]) + '; UTR3_Procession=' + str(
            total_freq_by_regions_on_chr[5]) + ';')

    print("--- %s seconds ---" % (time.time() - start_time))

    for sub_reg in range(0, 6):
        for nucleotide in range(0, 4):
            total_freq_by_regions[sub_reg][nucleotide] += total_freq_by_regions_on_chr[sub_reg][nucleotide]

    for sub_reg in range(0, 6):
        if chromosomes_dict[sub_reg] is None: chromosomes_dict[sub_reg] = []
        if nucleotide_dict[sub_reg] is None: nucleotide_dict[sub_reg] = []
        if counts_dict[sub_reg] is None: counts_dict[sub_reg] = []

        for i in range(0, 4):
            chromosomes_dict[sub_reg].append("chr" + str(chr_id))
            if i == 0: nucleotide_dict[sub_reg].append("C")
            if i == 1: nucleotide_dict[sub_reg].append("G")
            if i == 2: nucleotide_dict[sub_reg].append("A")
            if i == 3: nucleotide_dict[sub_reg].append("T")
            dist = total_freq_by_regions_on_chr[sub_reg][i] / (
                    total_freq_by_regions_on_chr[sub_reg][0] + total_freq_by_regions_on_chr[sub_reg][1] +
                    total_freq_by_regions_on_chr[sub_reg][2] + total_freq_by_regions_on_chr[sub_reg][3])
            counts_dict[sub_reg].append(dist)

# print('UTR5_Procession=' + str(total_freq_by_regions[0]) + '; UTR5=' + str(total_freq_by_regions[1]) + '; CDS=' + str(
#     total_freq_by_regions[2]) + '; introns=' +
#       str(total_freq_by_regions[3]) + '; UTR3=' + str(total_freq_by_regions[4]) + '; UTR3_Procession=' + str(
#     total_freq_by_regions[5]) + ';')

for sub_reg in range(0, 6):
    if chromosomes_dict[sub_reg] is None: chromosomes_dict[sub_reg] = []
    if nucleotide_dict[sub_reg] is None: nucleotide_dict[sub_reg] = []
    if counts_dict[sub_reg] is None: counts_dict[sub_reg] = []

    for i in range(0, 4):
        chromosomes_dict[sub_reg].append("Total")
        if i == 0: nucleotide_dict[sub_reg].append("C")
        if i == 1: nucleotide_dict[sub_reg].append("G")
        if i == 2: nucleotide_dict[sub_reg].append("A")
        if i == 3: nucleotide_dict[sub_reg].append("T")
        dist = total_freq_by_regions[sub_reg][i] / (
                total_freq_by_regions[sub_reg][0] + total_freq_by_regions[sub_reg][1] +
                total_freq_by_regions[sub_reg][2] + total_freq_by_regions[sub_reg][3])
        counts_dict[sub_reg].append(dist)

# x = dict(chromosome=["chr1", "chr1", "chr1", "chr1"], nucleotide=["C", "G", "A", "T"], count=[2, 3, 4, 5])


x = dict(chr=chromosomes_dict[0], nucleotide=nucleotide_dict[0], distribution=counts_dict[0])
df1 = pd.DataFrame(x)
fig = px.bar(df1, x="chr", y="distribution", color="nucleotide", title="UTR5_Procession", text_auto=',.0%'
             , color_discrete_map={'C': 'brown',
                                   'G': 'orange',
                                   'T': 'lightgrey',
                                   'A': 'grey'})

fig.write_image("images/UTR5_Procession.png", scale=3.0)

x = dict(chr=chromosomes_dict[1], nucleotide=nucleotide_dict[1], distribution=counts_dict[1])
df1 = pd.DataFrame(x)
fig = px.bar(df1, x="chr", y="distribution", color="nucleotide", title="UTR5", text_auto=',.0%',
             color_discrete_map={'C': 'brown',
                                 'G': 'orange',
                                 'T': 'lightgrey',
                                 'A': 'grey'})
fig.write_image("images/UTR5.png", scale=3.0)

x = dict(chr=chromosomes_dict[2], nucleotide=nucleotide_dict[2], distribution=counts_dict[2])
df1 = pd.DataFrame(x)
fig = px.bar(df1, x="chr", y="distribution", color="nucleotide", title="CDS", text_auto=',.0%',
             color_discrete_map={'C': 'brown',
                                 'G': 'orange',
                                 'T': 'lightgrey',
                                 'A': 'grey'})
fig.write_image("images/CDS.png", scale=3.0)

x = dict(chr=chromosomes_dict[3], nucleotide=nucleotide_dict[3], distribution=counts_dict[3])
df1 = pd.DataFrame(x)
fig = px.bar(df1, x="chr", y="distribution", color="nucleotide", title="introns", text_auto=',.0%',
             color_discrete_map={'C': 'brown',
                                 'G': 'orange',
                                 'T': 'lightgrey',
                                 'A': 'grey'})
fig.write_image("images/introns.png", scale=3.0)

x = dict(chr=chromosomes_dict[4], nucleotide=nucleotide_dict[4], distribution=counts_dict[4])
df1 = pd.DataFrame(x)
fig = px.bar(df1, x="chr", y="distribution", color="nucleotide", title="UTR3", text_auto=',.0%',
             color_discrete_map={'C': 'brown',
                                 'G': 'orange',
                                 'T': 'lightgrey',
                                 'A': 'grey'})
fig.write_image("images/UTR3.png", scale=3.0)

x = dict(chr=chromosomes_dict[5], nucleotide=nucleotide_dict[5], distribution=counts_dict[5])
df1 = pd.DataFrame(x)
fig = px.bar(df1, x="chr", y="distribution", color="nucleotide", title="UTR3_Procession", text_auto=',.0%',
             color_discrete_map={'C': 'brown',
                                 'G': 'orange',
                                 'T': 'lightgrey',
                                 'A': 'grey'})
fig.write_image("images/UTR3_Procession.png", scale=3.0)
