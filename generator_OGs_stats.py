# Author: "Lasha Bukhnikashvili"
#
# Description:
#   Calculates overlapping gene pairs count and other general stats in Genome Worker
#   like:
#       number of genes which is overlapped at least  once,
#       number of overlapping gene clusters
#       etc...
#
# Usage:
#   generator_OGs_stats.py <species> <annotation>
#
# Params (possible) to run:
#   species:  Homo sapiens / Rattus norvegicus / Mus musculus / Danio rerio / Drosophila melanogaster,
#               Caenorhabditis elegans
#   annotation: NCBI / Ensembl
#
# Example:
#   python generator_OGs_stats.py 'Homo sapiens' NCBI
#
# Output:
#   prints stats in console
import sys
from overlaps_by_CDS_finder import *
from generator_writer import *

start_time = time.time()

assert len(sys.argv) == 3
assert sys.argv[2] == 'NCBI' or sys.argv[2] == 'Ensembl'

species = SPECIES.from_string(sys.argv[1])
annotation_source = ANNOTATIONS.NCBI if sys.argv[2] == 'NCBI' else ANNOTATIONS.ENSEMBL

# Load only genes. we don't need gene specific fragments as we are calculating general stats
genome = GenomeWorker(species, annotation_source, ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS, SEQUENCE_LOAD.LOAD)
OGs_by_any_ISO_graph, _ = overlaps_finder(genome, SEARCH_METHOD.BY_GENE_BOUNDARIES)
og_count = 0
og_clusters_count = 0
clusters_count_by_length = {}
og_count_by_type = {OVERLAP_TYPE.DIVERGENT: {}, OVERLAP_TYPE.CONVERGENT: {}, OVERLAP_TYPE.DIFF_NESTED: {},
                    OVERLAP_TYPE.TANDEM: {}, OVERLAP_TYPE.SAME_NESTED: {}}
gene_clusters = OGs_by_any_ISO_graph.get_connected_clusters()

og_nuclear_mito_proteome = 0
og_mtDNA_mito_proteome = 0

for cluster in gene_clusters:
    clusters_count_by_length[len(cluster)] = clusters_count_by_length.get(len(cluster), 0) + 1
    if len(cluster) > 1:
        og_clusters_count += 1
        for g_id in cluster:
            is_mito = genome.get_feature_chromosomal_position(g_id)[0] == genome.chromosomes_count()
            if genome.is_gene_MITO(g_id) and is_mito: og_mtDNA_mito_proteome += 1
            if genome.is_gene_MITO(g_id) and not is_mito: og_nuclear_mito_proteome += 1
        for g_id1 in cluster:
            for g_id2 in cluster:
                if g_id1 == g_id2: continue
                gene1 = genome.feature_by_id(g_id1)
                gene2 = genome.feature_by_id(g_id2)
                if genome.get_overlap_type(gene1, gene2) == OVERLAP_TYPE.NONE: continue
                og_count_by_type[genome.get_overlap_type(gene1, gene2)][g_id1] = True
                og_count_by_type[genome.get_overlap_type(gene1, gene2)][g_id2] = True

    og_count += len(cluster)

total_genes = 0
positive_genes = 0
negative_genes = 0
isoforms_count = []
genes_length = []

for chr_id in range(1, genome.chromosomes_count() + 1):
    genes_cnt = genome.genes_count_on_chr(chr_id)

    # just count total genes by adding genes on this chromosome
    total_genes += genes_cnt

    # count genes by strand, by isoform numbers and by gene length for another statistical purposes
    for i in range(0, genes_cnt):
        gene = genome.gene_by_indexes(chr_id, i)
        positive_genes += 1 if gene.strand == '+' else 0
        negative_genes += 1 if gene.strand == '-' else 0
        transcripts = genome.get_transcripts_from_gene(gene.id)
        isoforms_count.append(len(transcripts))
        genes_length.append(gene.end - gene.start + 1)

output_data_path = f'generated_data/texts/OGs and stats [{species.short_name()}, {annotation_source.short_name()}].txt'

with open(output_data_path, 'w') as f:
    f.write(f"Organism: {species}\n")
    f.write(f"Annotation: {annotation_source.short_name()}\n")

    f.write(f'--------------------------Genome General data-----------------------------------\n\n')
    f.write(f"protein coding genes: {total_genes} ({genome.ignored_protein_coding_genes} filtered out)\n")
    f.write(f"\tFiltered out Genes: {str(genome.ignored_genes_by_types)}\n\n")

    f.write(f"Genes on Positive(+) Strand: {positive_genes}\n")
    f.write(f"Genes on Negative(-) Strand: {negative_genes}\n\n")

    f.write(
        f"isoforms per gene (mean, std, median, interquartile): {get_value_estimation(isoforms_count)}\n")
    f.write(
        f"gene length estimation (mean, std, median, interquartile): {get_value_estimation(genes_length)}\n\n")

    f.write(f'--------------------------Overlapping Genes by boundaries-----------------------------------\n\n')

    f.write(f"Overlapping Genes: {og_count}\n")
    f.write(
        f"\tmito proteome coding genes: {og_nuclear_mito_proteome + og_mtDNA_mito_proteome} = {og_nuclear_mito_proteome} (nuclear DNA) + {og_mtDNA_mito_proteome} (mtDNA)\n\n")
    f.write(f"Overlapping Gene clusters (>1 gene): {og_clusters_count}\n")

    cluster_items = list(clusters_count_by_length.items())
    cluster_items.sort(key=lambda x: x[0])

    for c_size, c_count in cluster_items:
        f.write(f"\t{c_size}-length clusters: {c_count}\n")

    f.write(f"\nOverlapping Genes by type of overlap: {og_clusters_count}\n")
    for ov_type, set_of_genes in og_count_by_type.items():
        f.write(f"\t{ov_type.short_name()}\t{len(set_of_genes)}\n")

    write_excel_formatted_general_data(f, genome, gene_clusters)

print("--- %s seconds ---" % (time.time() - start_time))
