# Author: "Lasha Bukhnikashvili"
#
# Description:
#   Calculates overlapping gene pairs count and other general stats in genome
#   like:
#       number of genes which is overlapped at least  once,
#       number of overlapping gene clusters
#       etc...
#
# Usage:
#   overlapping_genes_stats.py <species> <annotation>
#
# Params (possible) to run:
#   species:  Homo sapiens / Rattus norvegicus / Mus musculus / Danio rerio / Drosophila melanogaster,
#               Caenorhabditis elegans / Saccharomyces cerevisiae
#   annotation: NCBI / Ensembl
#
# Example:
#   python overlapping_genes_stats.py 'Homo sapiens' NCBI
#
# Output:
#   prints stats in console
import sys
import time

from genome_worker import *

start_time = time.time()

assert len(sys.argv) == 3

species = SPECIES.from_string(sys.argv[1])
annotation_source = ANNOTATIONS.NCBI if sys.argv[2] == 'NCBI' else ANNOTATIONS.ENSEMBL

total_genes = 0
positive_genes = 0
negative_genes = 0
og_count = 0
og_pairs_count = 0
og_clusters_count = 0

genes_by_clusters_length = [0] * 1000

# Load only genes. we don't need gene specific fragments as we are calculating general stats
human_genome = GenomeWorker(species, annotation_source, ANNOTATION_LOAD.GENES, SEQUENCE_LOAD.NOT_LOAD)

# excluding mitochondria
for chr_id in range(1, human_genome.chromosomes_count()):
    genes_cnt = human_genome.genes_count_on_chr(chr_id)
    total_overlapped_sequence = ""

    # make sure, that there is no same id over genes on chromosome
    for i in range(0, genes_cnt):
        for j in range(i + 1, genes_cnt):
            if human_genome.gene_by_ind(chr_id, i).id == human_genome.gene_by_ind(chr_id, j).id:
                print("lasha some issue!")

    # just count total genes by adding genes on this chromosome
    total_genes += genes_cnt

    # count genes by strand for another statistical purposes
    for i in range(0, genes_cnt):
        gene = human_genome.gene_by_ind(chr_id, i)
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
            gene_a = human_genome.gene_by_ind(chr_id, i)
            gene_b = human_genome.gene_by_ind(chr_id, j)
            if GenomeWorker.are_genes_overlapped(gene_a, gene_b):

                # for clustering computation
                old_cluster_index = cluster_indexes[i]
                new_cluster_index = cluster_indexes[j]
                for k in range(0, genes_cnt):
                    if cluster_indexes[k] == old_cluster_index:
                        cluster_indexes[k] = new_cluster_index

                # counting just overlapping genes pair
                og_pairs_count += 1

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
print("Number of genes: " + str(total_genes))
print("Number of genes on Positive(+) Strand: " + str(positive_genes))
print("Number of genes on Negative(-) Strand: " + str(negative_genes))
print("")

print("Number of OG: " + str(og_count))
print("Number of OG pairs: " + str(og_pairs_count))
print("Number of OG clusters with >=2 gene: " + str(og_clusters_count))
print("")

print("Genes by clusters length:")
for i in range(0, 1000):
    if genes_by_clusters_length[i] > 0:
        print(str(i) + "-length clusters: " + str((genes_by_clusters_length[i] // i)))
print("")

print("--- %s seconds ---" % (time.time() - start_time))
