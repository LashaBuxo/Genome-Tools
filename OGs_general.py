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
#   OGs_general.py <species> <annotation>
#
# Params (possible) to run:
#   species:  Homo sapiens / Rattus norvegicus / Mus musculus / Danio rerio / Drosophila melanogaster,
#               Caenorhabditis elegans
#   annotation: NCBI / Ensembl
#
# Example:
#   python OGs_general.py 'Homo sapiens' NCBI
#
# Output:
#   prints stats in console
import sys
import time

from genome_worker_enums import *

start_time = time.time()

assert len(sys.argv) == 3
assert sys.argv[2] == 'NCBI' or sys.argv[2] == 'Ensembl'

species = SPECIES.from_string(sys.argv[1])
annotation_source = ANNOTATIONS.NCBI if sys.argv[2] == 'NCBI' else ANNOTATIONS.ENSEMBL

total_genes = 0
positive_genes = 0
negative_genes = 0
og_count = 0
og_clusters_count = 0

genes_by_clusters_length = [0] * 10000

# Load only genes. we don't need gene specific fragments as we are calculating general stats
genome = GenomeWorker(species, annotation_source, ANNOTATION_LOAD.GENES,
                      SEQUENCE_LOAD.NOT_LOAD)

x = 0
# excluding mitochondria
for chr_id in range(1, genome.chromosomes_count() + 1):
    genes_cnt = genome.genes_count_on_chr(chr_id)
    total_overlapped_sequence = ""

    # just count total genes by adding genes on this chromosome
    total_genes += genes_cnt

    # count genes by strand for another statistical purposes
    for i in range(0, genes_cnt):
        gene = genome.gene_by_ind(chr_id, i)
        if gene.strand == '+':
            positive_genes = positive_genes + 1
        if gene.strand == '-':
            negative_genes = negative_genes + 1

    # initialise clusters, for later cluster counting
    cluster_indexes = [0] * genes_cnt
    for i in range(0, genes_cnt):
        cluster_indexes[i] = i

    # for every different pair of genes
    for i in range(0, genes_cnt):
        for j in range(i + 1, genes_cnt):
            gene_a = genome.gene_by_ind(chr_id, i)
            gene_b = genome.gene_by_ind(chr_id, j)

            if GenomeWorker.are_features_overlapped(gene_a, gene_b):

                # for clustering computation
                old_cluster_index = cluster_indexes[i]
                new_cluster_index = cluster_indexes[j]
                for k in range(0, genes_cnt):
                    if cluster_indexes[k] == old_cluster_index:
                        cluster_indexes[k] = new_cluster_index

    # calculating genes count for each cluster
    genes_in_cluster = [0] * genes_cnt
    for k in range(0, genes_cnt):
        genes_in_cluster[cluster_indexes[k]] += 1

    # calculating genes count which at least once overlapped to different gene
    for i in range(0, genes_cnt):
        if genes_in_cluster[cluster_indexes[i]] > 1:
            og_count += 1
        if genes_in_cluster[i] > 1:
            og_clusters_count += 1

    for i in range(0, genes_cnt):
        genes_by_clusters_length[genes_in_cluster[cluster_indexes[i]]] += 1

# print stats
assert total_genes == genome.imported_protein_coding_genes

print("Organism: " + str(species))
print("Annotation: " + annotation_source.name)
print("protein coding genes: " + str(total_genes) + " (" + str(
    genome.ignored_protein_coding_genes) + " filtered out)")
print("")
print("Genes on Positive(+) Strand: " + str(positive_genes))
print("Genes on Negative(-) Strand: " + str(negative_genes))
print("")

print("Overlapping Genes: " + str(og_count))
print("Overlapping Gene clusters (>1 gene): " + str(og_clusters_count))

for i in range(2, 1000):
    if genes_by_clusters_length[i] > 0:
        print("     " + str(i) + "-length clusters: " + str((genes_by_clusters_length[i] // i)))
print("")

print("--- %s seconds ---" % (time.time() - start_time))
