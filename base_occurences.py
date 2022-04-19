import genome_tools as genome
from gffutils import Feature
import plotly.express as px
import pandas as pd


def add_matrix(a, b, parts):
    for i in range(0, 6):
        for j in range(0, parts):
            for k in range(0, 4):
                a[i][j][k] += b[i][j][k]
    return a


# subregions of interest: UTR'5_Procession, UTR'5, inner CDSs, inner Introns, UTR'3, UTR'3_Procession
# each region divided into k-part
# for each part (a,c,g,t) is calculated

# number of parts
k = 50

# distributions[region][part][base] = [6][k][4]
total_occurrences = []
for i in range(0, 6):
    occurrences_on_region = []
    for j in range(0, k):
        base_distributions_on_part = [0, 0, 0, 0]
        occurrences_on_region.append(base_distributions_on_part)
    total_occurrences.append(occurrences_on_region)

# excluding mitochondria
for chr_id in range(1, genome.number_of_chromosomes):
    print(chr_id)
    genome.preprocess_annotation_for_chr(chr_id)
    genes_cnt = genome.genes_count_on_chr(chr_id)

    for i in range(0, genes_cnt):
        gene = genome.gene_by_ind(chr_id, i)
        gene_occurrences = genome.analyze_gene_occurrences_by_parts(chr_id, gene, k)
        add_matrix(total_occurrences, gene_occurrences, k)

region_by_part = []
distribution = []
nucleotide = []
lineind = []

region_id = 0
for region in range(0, 6):
    if region != 0: region_id += k
    for part in range(0, k):

        region_by_part.append(region_id + part)
        region_by_part.append(region_id + part)
        region_by_part.append(region_id + part)
        region_by_part.append(region_id + part)

        nucleotide.append('C')
        nucleotide.append('G')
        nucleotide.append('A')
        nucleotide.append('T')

        dist1 = total_occurrences[region][part][0] / (
                total_occurrences[region][part][0] + total_occurrences[region][part][1] +
                total_occurrences[region][part][2] + total_occurrences[region][part][3])

        dist2 = total_occurrences[region][part][1] / (
                total_occurrences[region][part][0] + total_occurrences[region][part][1] +
                total_occurrences[region][part][2] + total_occurrences[region][part][3])

        dist3 = -total_occurrences[region][part][2] / (
                total_occurrences[region][part][0] + total_occurrences[region][part][1] +
                total_occurrences[region][part][2] + total_occurrences[region][part][3])

        dist4 = -total_occurrences[region][part][3] / (
                total_occurrences[region][part][0] + total_occurrences[region][part][1] +
                total_occurrences[region][part][2] + total_occurrences[region][part][3])

        distribution.append(dist1)
        distribution.append(dist2)
        distribution.append(dist3)
        distribution.append(dist4)

        ind = 0
        if region == 2: ind = 4
        if region == 3: ind = 8
        if region > 3: ind = 12

        lineind.append(ind)
        lineind.append(ind + 1)
        lineind.append(ind + 2)
        lineind.append(ind + 3)

df = pd.DataFrame(dict(
    regional_parts=region_by_part,
    distribution=distribution,
    nucleotide=nucleotide,
    lineind=lineind
))

fig = px.line(df, x="regional_parts", y="distribution", line_group=lineind,
              color="nucleotide")
#"Each region is divided into " + str(k) + " part"

fig.write_image("results/average_gene_outline.png", scale=4.0)
