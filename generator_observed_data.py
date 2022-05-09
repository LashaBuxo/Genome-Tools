# Author: "Lasha Bukhnikashvili"
# Usage:
#   generator_observed_data.py <species> <annotation>
#
# Params to run:
#   species: Homo_sapiens / Mus_musculus
#   annotation: NCBI / Ensembl
#
# Description:
#   Finds gene clusters (size >1), where 2 gene will be in same cluster if these conditions meet:
#       1) gene1 boundaries are overlapped to gene2 boundaries by specific OVERLAP TYPE
#           a) different strand + convergent
#           b) different strand + nested
#           c) different strand + divergent
#           d) same strand + nested
#           e) same strand + tandem (neither gene is nested)
#       2) at least one transcript(mRNA) of gene1 is in overlap to at least one
#          transcript(mRNA) of gene2 by CDS region (= exons excluding UTR'5 and UTR'3)
#
#   For each  OVERLAP TYPE this data (below) is exported as single file (in ./generated_data/observed/ directory):
#       For every cluster merged overlapped CDS sequence parameters:
#           a. GC content
#           b. HD content - High-Degeneracy amino acids occurrences [Feature 2]
#           *. other parameters... like, di-nucleotide and amino acid occurrences,
#               which are not used later, for comparison to the control group parameters.
#       For each gene in every cluster:
#           a. expression profile - level from all tissue.
#
#
#   Each observed data for specific OVERLAP TYPE is compared to Control Group:
#       For all overlapped genes (in any cluster) we randomly chose non-overlapped gene in a vicinity
#       with almost similar length and for every gene overall expression level (in tissues) was calculated.
#       These expression levels are compared to the overlapped genes expression levels.
#       the overall expression levels of overlapped genes using

from math import log2
from analyzer_worker import *
from genome_worker import *
from genome_worker_enums import *
import time
import sys
import plotly.graph_objects as go
import plotly.express as px


def add_gene_sequence_record(gene_feature: Feature, overlap_type: OVERLAP_TYPE, seq_record, orf):
    if not gene_sequence_records.__contains__(gene_feature.id):
        gene_sequence_records[gene_feature.id] = []
    gene_sequence_records[gene_feature.id].append((overlap_type, (seq_record, orf)))

def get_median_expression_value(expression_data):
    return 0
start_time = time.time()

# region variables and input management
assert len(sys.argv) == 3
assert sys.argv[2] == 'NCBI' or sys.argv[2] == 'Ensembl'

species = SPECIES.from_string(sys.argv[1])
annotation = ANNOTATIONS.NCBI if sys.argv[2] == 'NCBI' else ANNOTATIONS.ENSEMBL

# Load necessary annotation data
genome = GenomeWorker(species, annotation, ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_CDS, SEQUENCE_LOAD.LOAD,
                      EXPRESSIONS_LOAD.LOAD)
observed_data_graph = AnalyzerGraph()
overlapped_isoforms = {}
gene_sequence_records = {}

# endregion

#
# def find_non_overlapped_random_gene(gene_len):
#     # approx. 10% difference
#     genes_with_similar_length = []
#     for chr_id in range(1, genome.chromosomes_count() + 1):
#         genes_cnt = genome.genes_count_on_chr(chr_id)
#         for i in range(0, genes_cnt):
#             random_gene = genome.gene_by_ind(chr_id, i)
#             # we are searching only in non-overlapped genes
#             if stored_gene_features.__contains__(random_gene.id): continue
#             random_gene_length = random_gene.end - random_gene.start + 1
#             if abs(random_gene_length - gene_len) < 10 * gene_len // 100:
#                 genes_with_similar_length.append(random_gene)
#     assert len(genes_with_similar_length) > 0
#     return random.choice(genes_with_similar_length)


# region observed data ( OGs by CDS) collector


for chr_id in range(1, genome.chromosomes_count() + 1):
    genes_cnt = genome.genes_count_on_chr(chr_id)

    for i in range(0, genes_cnt):
        for j in range(i + 1, genes_cnt):
            gene_A = genome.gene_by_indexes(chr_id, i)
            gene_B = genome.gene_by_indexes(chr_id, j)

            overlap_type = genome.get_overlap_type(gene_A, gene_B)

            # if genes don't overlap, there is no point in searching overlaps by fragments (CDS in this case)
            if overlap_type is OVERLAP_TYPE.NONE: continue

            max_overlapped_segments = []
            max_overlapped_length = 0
            order_of_frames_matched = None

            transcripts_A = genome.get_transcripts_from_gene(gene_A.id)
            transcripts_B = genome.get_transcripts_from_gene(gene_B.id)

            for transcript_a in transcripts_A:
                for transcript_b in transcripts_B:
                    segments = genome.get_fragments_overlap_between_transcripts(transcript_a.id, transcript_b.id)
                    length = 0
                    for ind in range(len(segments)):
                        assert segments[ind][0][1] >= segments[ind][0][0]
                        length += segments[ind][0][1] - segments[ind][0][0] + 1

                    if length > 0:
                        overlapped_isoforms[transcript_a.id] = True
                        overlapped_isoforms[transcript_b.id] = True

                    if length > max_overlapped_length:
                        max_overlapped_length = length
                        max_overlapped_segments = segments
                        x = transcript_a
                        y = transcript_b

            if max_overlapped_length == 0: continue

            # we should add edge in analyzer graph, along with corresponding shared sequences
            edge = AnalyzerGraph.GraphEdge(gene_A.id, gene_B.id, overlap_type)
            observed_data_graph.add_edge(edge)

            for interval, frames in max_overlapped_segments:
                assert interval[1] >= interval[0]
                if gene_A.strand == gene_B.strand:
                    seq = genome.retrieve_interval_sequence(chr_id, interval[0], interval[1], gene_A.strand)
                    # add_gene_sequence_record(gene_A, overlap_type, seq, frames[0])
                    # add_gene_sequence_record(gene_B, overlap_type, seq, frames[1])
                else:
                    seq1 = genome.retrieve_interval_sequence(chr_id, interval[0], interval[1], gene_A.strand)
                    seq2 = genome.retrieve_interval_sequence(chr_id, interval[0], interval[1], gene_B.strand)

                    # add_gene_sequence_record(gene_A, overlap_type, seq1, frames[0])
                    # add_gene_sequence_record(gene_B, overlap_type, seq2, frames[1])
# endregion

# region generating control DATA

random.seed(1)

control_genes = []
control_genes_expression_values = []
number_of_genes_in_any_subtype = stored_gene_features.keys()

control_genes_expression_values2 = []
control_genes_expression_values3 = []

for (_, gene_feature) in stored_gene_features.items():
    gene_length = gene_feature.end - gene_feature.start + 1

    control_gene = find_non_overlapped_random_gene(gene_length)
    control_genes.append(control_gene)

    arr = []
    avg_expression = 0
    exp_profile = genome.get_gene_expression_data_by_id(control_gene.id).items()
    for (_, expression) in exp_profile:
        avg_expression += log2(1 + expression)
    avg_expression = 0 if len(exp_profile) == 0 else avg_expression / len(exp_profile)
    control_genes_expression_values.append(avg_expression)

for (_, gene_feature) in stored_gene_features.items():
    gene_length = gene_feature.end - gene_feature.start + 1

    control_gene = find_non_overlapped_random_gene(gene_length)

    arr = []
    avg_expression = 0
    exp_profile = genome.get_gene_expression_data_by_id(control_gene.id).items()
    for (_, expression) in exp_profile:
        avg_expression += log2(1 + expression)
    avg_expression = 0 if len(exp_profile) == 0 else avg_expression / len(exp_profile)
    control_genes_expression_values2.append(avg_expression)

for (_, gene_feature) in stored_gene_features.items():
    gene_length = gene_feature.end - gene_feature.start + 1

    control_gene = find_non_overlapped_random_gene(gene_length)

    arr = []
    avg_expression = 0
    exp_profile = genome.get_gene_expression_data_by_id(control_gene.id).items()
    for (_, expression) in exp_profile:
        avg_expression += log2(1 + expression)
    avg_expression = 0 if len(exp_profile) == 0 else avg_expression / len(exp_profile)
    control_genes_expression_values3.append(avg_expression)

# endregion

expression_comparative_fig = go.Figure()
expression_comparative_fig.update_layout(yaxis_title='log2(nTPM+1)')
expression_comparative_fig.add_trace(go.Box(y=control_genes_expression_values, name='Control 1'))
expression_comparative_fig.add_trace(go.Box(y=control_genes_expression_values2, name='Control 2'))
expression_comparative_fig.add_trace(go.Box(y=control_genes_expression_values2, name='Control 3'))

# region divide genes into clusters for each OVERLAP TYPE
for overlap_type in OVERLAP_TYPE.get_overlap_types():
    gene_clusters = observed_data_graph.get_connected_clusters(overlap_type)

    overall_stats = AnalyzerData()
    cluster_stats = []
    for cluster in gene_clusters:
        stats = AnalyzerData()

        for gene_id in cluster:
            gene = genome.feature_by_id(gene_id)
            transcripts = genome.get_transcripts_from_gene(gene.id)

            for transcript in transcripts:
                expression_values = genome.get_transcript_expression_values(transcript.id)
                overall_stats.analyze_expression_profile(expression_values)

                if not overlapped_isoforms.__contains__(transcript.id):
                    expres

            for record_type, (sequence, frame) in gene_sequence_records[gene_id]:
                if record_type != overlap_type: continue
                stats.analyze_sequence_stats(sequence, frame)
                overall_stats.analyze_sequence_stats(sequence, frame)

        cluster_stats.append(stats)

    # example path './results/overlaps - OGs by CDS [Human, Ensembl, ] .txt'
    output_data_path = f'generated_data/observed/OGs by CDS [{species.short_name()}, {annotation.short_name()}, ' \
                       f'{overlap_type.short_name()}].txt'

    # additional clustering
    overlapped_genes = 0
    overlapped_genes_by_cluster_size = {}
    for cluster in gene_clusters:
        cluster_size = len(cluster)
        overlapped_genes += cluster_size
        if not overlapped_genes_by_cluster_size.__contains__(cluster_size):
            overlapped_genes_by_cluster_size[cluster_size] = 0
        overlapped_genes_by_cluster_size[len(cluster)] += 1

    with open(output_data_path, 'w') as f:
        keys = list(overlapped_genes_by_cluster_size.keys())
        keys.sort()

        f.write(f'Total overlapped (by CDS) genes: {overlapped_genes} ={" 0" if len(keys) == 0 else ""}')

        for i in range(0, len(keys)):
            f.write(f'{" +" if i > 0 else ""} {overlapped_genes_by_cluster_size[keys[i]]}x[{keys[i]}]')
        f.write('\n')

        f.write(f'Total overlapped (by CDS) sequence length: {overall_stats.analyzed_sequence_length}\n\n')
        f.write(f'--------------------------overall mean parameters---------------------------------\n\n')
        f.write(overall_stats.get_formatted_sequence_stats(indent_tabs=1, use_percents=True))
        f.write(overall_stats.get_formatted_expression_profile(indent_tabs=1, use_percents=False))

        f.write(f'--------------------------clusters and their parameters-----------------------------------\n\n')

        for i in range(len(gene_clusters)):
            cluster = gene_clusters[i]
            stats = cluster_stats[i]

            f.write(
                f'cluster #{i + 1}  |  size: {len(cluster)}  |  total overlapped sequence: {stats.analyzed_sequence_length} \n')
            f.write(f'\tGC content: {"{:.2f}".format(stats.get_mean_gc_content())}%\n')
            f.write(f'\tHD content: {"{:.2f}".format(stats.get_mean_hd_content())}%\n\n')

            for gene_index in range(len(cluster)):
                gene_id = cluster[gene_index]
                gene = stored_gene_features[gene_id]

                gene_name = GenomeWorker.get_gene_symbol(gene)
                gene_desc = GenomeWorker.get_gene_description(gene)

                f.write(f'\tGene {gene_index + 1}: {gene_name} | {gene_desc} | {gene_id}\n')

                single_gene_stats = AnalyzerData()
                gene_expression = genome.get_gene_expression_data_by_id(gene_id)
                single_gene_stats.analyze_expression_profile(gene_expression)
                f.write(single_gene_stats.get_formatted_expression_profile(indent_tabs=2, use_percents=False))
            f.write('\n')

    observed_overall_expressions = []
    for cluster in gene_clusters:
        for gene_id in cluster:
            arr = []
            avg_expression = 0
            exp_profile = genome.get_gene_expression_data_by_id(gene_id).items()
            for (_, expression) in exp_profile:
                avg_expression += log2(1 + expression)
            avg_expression = 0 if len(exp_profile) == 0 else avg_expression / len(exp_profile)
            observed_overall_expressions.append(avg_expression)

    Tissue_names = [
        'cerebral cortex', 'olfactory bulb', 'hippocampal formation', 'amygdala', 'basal ganglia', 'thalamus',
        'hypothalamus', 'midbrain', 'pons', 'cerebellum', 'medulla oblongata', 'spinal cord', 'white matter',

        'adipose tissue', 'adrenal gland', 'appendix', 'bone marrow', 'breast', 'cervix', 'choroid plexus', 'colon',
        'duodenum', 'endometrium 1', 'epididymis', 'esophagus', 'fallopian tube', 'gallbladder', 'heart muscle',
        'kidney', 'liver', 'lung', 'lymph node', 'ovary', 'pancreas', 'parathyroid gland', 'pituitary gland',
        'placenta', 'prostate', 'rectum', 'retina', 'salivary gland', 'seminal vesicle', 'skeletal muscle', 'skin 1',
        'small intestine', 'smooth muscle', 'spleen', 'stomach 1', 'testis', 'thymus', 'thyroid gland', 'tongue',
        'tonsil', 'urinary bladder', 'vagina', ]

    expression_profile_data = []
    cluster_names = []
    cluster_ind = 0

    median_lasha = []
    for cluster in gene_clusters:
        cluster_ind += 1
        gene_ind = 0
        for gene_id in cluster:
            gene_ind += 1
            gene_expression_by_tissue = []
            expression_data = genome.get_gene_expression_data_by_id(gene_id)

            lasha_sum = 0

            for tissue_ind in range(len(Tissue_names)):
                tissue_name = Tissue_names[tissue_ind]
                gene_expression_by_tissue.append(0)
                if expression_data.__contains__(f'Tissue RNA - {tissue_name} [nTPM]'):
                    gene_expression_by_tissue[tissue_ind] = log2(
                        expression_data[f'Tissue RNA - {tissue_name} [nTPM]'] + 1)
                    lasha_sum += log2(
                        expression_data[f'Tissue RNA - {tissue_name} [nTPM]'] + 1)
            median_lasha.append(lasha_sum)

            expression_profile_data.append(gene_expression_by_tissue)
            cluster_names.append(f'cluster {cluster_ind} [{gene_ind}]')

    # median_lasha.sort()
    # median_lasha_value = median_lasha[len(median_lasha) // 2]
    #
    # arr = []
    # for i in range(len(cluster_names)):
    #     arr.append((expression_profile_data[i], cluster_names[i], median_lasha_value))
    # arr.sort(key=choose_better_expression)
    #
    # expression_profile_data = []
    # cluster_names = []
    # for i in range(len(arr)):
    #     expression_profile_data.append(arr[i][0])
    #     cluster_names.append(arr[i][1])

    expression_profile_fig = px.imshow(expression_profile_data,
                                       labels=dict(x="Tissues", y="Cluster Index", color="Intensity"),
                                       x=Tissue_names,
                                       y=cluster_names,
                                       # aspect="auto"
                                       )
    expression_profile_fig.update_layout(title="Correlation heatmap",
                                         yaxis_nticks=len(expression_profile_data),
                                         width=1200,
                                         height=min(1200, 1200 // 50 * len(cluster_names)),
                                         )
    expression_profile_fig.update_xaxes(side="top")

    profile_path = f'generated_data/summarized/Expression Profile [{species.short_name()}, {annotation.short_name()},' \
                   f'{overlap_type.short_name()}].png'
    expression_profile_fig.write_image(profile_path, scale=2.0)

    expression_comparative_fig.add_trace(go.Box(y=observed_overall_expressions, name=overlap_type.short_name()))
# endregion

fig_path = f'generated_data/summarized/Expression Comparison [{species.short_name()}, {annotation.short_name()}].png'

expression_comparative_fig.write_image(fig_path, scale=4.0)
print("--- script was running for %s seconds. ---" % (time.time() - start_time))
