# Author: "Lasha Bukhnikashvili"
#
# Description:
#   Calculates 'Average' gene outline in Genome Worker.
#   Genes is divided into 6 region: UTR'5_Procession, UTR'5,
#   inner CDSs, inner Introns, UTR'3, UTR'3_Procession and each
#   region is divided into k subregions, where ATCG
#   occurrences are calculated. Based each nucleotide occurrences,
#   graph is built from 6*k point where X axis corresponds to points,
#   and Y axis corresponds  to nucleotide occurences.
#
# Usage:
#   generator_gene_outline.py <species> <annotation> <subregions> <tail_length>
#
# Params (possible) to run:
#   species:  Homo sapiens / Rattus norvegicus / Mus musculus / Danio rerio / Drosophila melanogaster,
#               Caenorhabditis elegans
#   annotation: NCBI / Ensembl
#   subregions: k - number of subregions within each region
#   procession_length: l - length of processions at the end of gene (before UTR5' or after UTR3')
#
# Example:
#   python generator_gene_outline.py Rattus_norvegicus Ensembl 50 1000
#
# Output:
#   created image in ./results/images folder

from genome_worker import *
import plotly.express as px
import pandas as pd
import time
import sys


def add_matrix(a, b, parts):
    for i in range(0, 6):
        for j in range(0, parts):
            for k in range(0, 4):
                a[i][j][k] += b[i][j][k]
    return a


start_time = time.time()

assert len(sys.argv) == 5
assert sys.argv[2] == 'NCBI' or sys.argv[2] == 'Ensembl'

species = SPECIES.from_string(sys.argv[1])
annotation = ANNOTATIONS.NCBI if sys.argv[2] == 'NCBI' else ANNOTATIONS.ENSEMBL
k = int(sys.argv[3])
procession_length = int(sys.argv[4])

# distributions[region][part][base] = [6][k][4]
total_occurrences = []
for i in range(0, 6):
    occurrences_on_region = []
    for j in range(0, k):
        base_distributions_on_part = [0, 0, 0, 0]
        occurrences_on_region.append(base_distributions_on_part)
    total_occurrences.append(occurrences_on_region)

# Load necessary annotation data, we need all fragments (cds,exons,utrs) to calculate occurrences in each region
genome = GenomeWorker(species, annotation, ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_FRAGMENTS,
                      SEQUENCE_LOAD.LOAD)

# excluding mitochondria DNA
for chr_id in range(1, genome.chromosomes_count() + 1):
    print("calculating occurrences in chromosome: " + str(chr_id))
    genes_cnt = genome.genes_count_on_chr(chr_id)

    for i in range(0, genes_cnt):
        gene = genome.gene_by_indexes(chr_id, i)
        gene_occurrences = genome.analyze_gene_occurrences_by_parts(chr_id, gene, k, procession_length,
                                                                    TRANSCRIPT_CRITERIA.NONE)
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
        total = total_occurrences[region][part][0] + total_occurrences[region][part][1] + \
                total_occurrences[region][part][2] + total_occurrences[region][part][3]
        if total == 0:
            print(region)
            print(part)
            print(total_occurrences)
        dist1 = total_occurrences[region][part][0] / total

        dist2 = total_occurrences[region][part][1] / total

        dist3 = total_occurrences[region][part][2] / total - 0.5

        dist4 = total_occurrences[region][part][3] / total - 0.5

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

# build graph
fig = px.line(df, x="regional_parts", y="distribution", labels={
    "regional_parts": "",
    "distribution": "Frequency",
    'nucleotide': "bases"
}, title='Representation of genes in ' + str(species) + ' by bases', line_group=lineind, color="nucleotide")

fig.update_layout(
    yaxis=dict(
        tickmode='array',
        tickvals=[-0.3, -0.2, -0.1, 0.1, 0.2, 0.3],
        ticktext=['20%', '30%', '40%', "10%", "20%", '30%']
    )
)

fig.update_layout(
    xaxis=dict(
        tickmode='array',
        tickvals=[k - 0.5, 2 * k - 0.5, 3 * k - 0.5, 4 * k - 0.5, 5 * k - 0.5],
        ticktext=[str(k), str(2 * k), str(3 * k), str(4 * k), str(5 * k)]
    )
)

# bottom text
fig.add_annotation(dict(xref='paper', yref='paper', x=0.5, y=-0.1, xanchor='center', yanchor='top',
                        text='Annotation = ' + sys.argv[2] + ', Subregions = ' + sys.argv[3] +
                             ', Tail length = ' + sys.argv[4],
                        font=dict(family='Arial', size=12, color='rgb(150,150,150)'),
                        showarrow=False))

# region texts
fig.add_annotation(x=k / 2, y=0, text="<b>(len " + str(procession_length) + ") >5'</b>", showarrow=False, yshift=10,
                   font=dict(family="Courier New, monospace", size=9, color='rgb(150,150,150)'))

fig.add_annotation(x=3 * k / 2, y=0, text="<b>UTR 5'</b>", showarrow=False, yshift=10,
                   font=dict(family="Courier New, monospace", size=9, color="#ff7f0e"),
                   # bordercolor="#c7c7c7", borderwidth=2, borderpad=4, bgcolor="#ff7f0e", opacity=0.8
                   )

fig.add_annotation(x=5 * k / 2, y=0, text="<b>CDSs</b>", showarrow=False, yshift=10,
                   font=dict(family="Courier New, monospace", size=9, color="#0CE83F"),
                   # bordercolor="#c7c7c7", borderwidth=2, borderpad=4, bgcolor="#0CE83F", opacity=0.8
                   )

fig.add_annotation(x=7 * k / 2, y=0, text="<b>INTRONs</b>", showarrow=False, yshift=10,
                   font=dict(family="Courier New, monospace", size=9, color="#000000"),
                   # bordercolor="#c7c7c7", borderwidth=2, borderpad=4, opacity=0.8
                   )

fig.add_annotation(x=9 * k / 2, y=0, text="<b>UTR 3'</b>", showarrow=False, yshift=10,
                   font=dict(family="Courier New, monospace", size=9, color="#ff7f0e"),
                   # bordercolor="#c7c7c7", borderwidth=2, borderpad=4, bgcolor="#ff7f0e", opacity=0.8
                   )

fig.add_annotation(x=11 * k / 2, y=0, text="<b>3'< (len " + str(procession_length) + ")</b>", showarrow=False,
                   yshift=10,
                   font=dict(family="Courier New, monospace", size=9, color='rgb(150,150,150)'))

# lines
fig.add_hline(y=0, line_width=1.5, line_color='rgb(150,150,150)')
fig.add_hline(y=-0.4, line_width=0.1, line_color='rgb(150,150,150)')
fig.add_hline(y=0.4, line_width=0.1, line_color='rgb(150,150,150)')

fig.add_vline(x=1 * k - 0.5, line_width=0.75, line_dash="dash", line_color='rgb(150,150,150)')
fig.add_vline(x=2 * k - 0.5, line_width=0.75, line_color='rgb(150,150,150)')
fig.add_vline(x=3 * k - 0.5, line_width=0.75, line_color='rgb(150,150,150)')
fig.add_vline(x=4 * k - 0.5, line_width=0.75, line_color='rgb(150,150,150)')
fig.add_vline(x=5 * k - 0.5, line_width=0.75, line_dash="dash", line_color='rgb(150,150,150)')

fig.update_traces(line=dict(width=1))

out_path = f'results/images/gene outline ({species.short_name()}, {annotation.short_name()}, k={sys.argv[3]}, tail={sys.argv[4]}).png'

fig.write_image(out_path, scale=4.0)

print("--- script was running for %s seconds. ---" % (time.time() - start_time))
