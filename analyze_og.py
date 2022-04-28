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
#   analyze_og.py <species> <annotation>
#
# Params (possible) to run:
#   species:  Homo sapiens / Rattus norvegicus / Mus musculus / Danio rerio / Drosophila melanogaster,
#               Caenorhabditis elegans
#   annotation: NCBI / Ensembl
#
# Example:
#   python analyze_og.py 'Homo sapiens' NCBI
#
# Output:
#   prints stats in console
import sys
import time

from genome_worker import *

start_time = time.time()

assert len(sys.argv) == 3
assert sys.argv[2] == 'NCBI' or sys.argv[2] == 'Ensembl'

species = SPECIES.from_string(sys.argv[1])
annotation_source = ANNOTATIONS.NCBI if sys.argv[2] == 'NCBI' else ANNOTATIONS.ENSEMBL

total_genes = 0
positive_genes = 0
negative_genes = 0
og_count = 0
og_pairs_count = 0
og_clusters_count = 0

genes_by_clusters_length = [0] * 10000

# Load only genes. we don't need gene specific fragments as we are calculating general stats
genome = GenomeWorker(species, annotation_source, ANNOTATION_LOAD.GENES,
                      SEQUENCE_LOAD.NOT_LOAD)

x = 0
# excluding mitochondria
for chr_id in range(1, genome.chromosomes_count()):
    genes_cnt = genome.genes_count_on_chr(chr_id)
    total_overlapped_sequence = ""

    removed_gene = [False] * genes_cnt

    # make sure, that there is no readthrough genes
    for i in range(0, genes_cnt):
        gene = genome.gene_by_ind(chr_id, i)
        if not gene.attributes.__contains__('description'):
            removed_gene[i] = True
            continue

        if not gene.attributes.__contains__('Name'):
            removed_gene[i] = True

        if gene.attributes['description'][0].__contains__('readthrough'):
            removed_gene[i] = True

    # make sure, that there is no same name over genes on chromosome
    # for i in range(0, genes_cnt):
    #     for j in range(i + 1, genes_cnt):
    #         if removed_gene[i] or removed_gene[j]: continue
    #         gene_A = genome.gene_by_ind(chr_id, i)
    #         gene_B = genome.gene_by_ind(chr_id, j)
    #
    #         if gene_A.attributes['Name'] == gene_B.attributes['Name']:
    #             removed_gene[i] = True
    #             # removed_gene[j] = True
    #             continue
    #
    #         desc_parts_1 = gene_A.attributes['description'][0].split("Acc:")
    #         ncbi_id_1 = desc_parts_1[1][0:(len(desc_parts_1[1]) - 1)]
    #
    #         desc_parts_2 = gene_B.attributes['description'][0].split("Acc:")
    #         ncbi_id_2 = desc_parts_2[1][0:(len(desc_parts_2[1]) - 1)]
    #
    #         if ncbi_id_1 == ncbi_id_2:
    #             removed_gene[j] = True

    removed_genes = 0
    for i in range(0, genes_cnt):
        if removed_gene[i]: removed_genes += 1

    # just count total genes by adding genes on this chromosome
    total_genes += genes_cnt - removed_genes

    # count genes by strand for another statistical purposes
    for i in range(0, genes_cnt):
        gene = genome.gene_by_ind(chr_id, i)
        if removed_gene[i]: continue
        if gene.strand == '+':
            positive_genes = positive_genes + 1
        if gene.strand == '-':
            negative_genes = negative_genes + 1

    # initialise clusters, for later cluster counting
    cluster_indexes = [0] * genes_cnt
    for i in range(0, genes_cnt):

        if removed_gene[i]: continue
        cluster_indexes[i] = i

    # for every different pair of genes
    for i in range(0, genes_cnt):
        if removed_gene[i]: continue
        for j in range(i + 1, genes_cnt):
            if removed_gene[j]: continue
            gene_a = genome.gene_by_ind(chr_id, i)
            gene_b = genome.gene_by_ind(chr_id, j)

            if GenomeWorker.are_features_overlapped(gene_a, gene_b):
                # feature_b = genome.find_gene_transcript_by_criteria(gene_b, TRANSCRIPT_CRITERIA.LONGEST)
                # feature_a = genome.find_gene_transcript_by_criteria(gene_a, TRANSCRIPT_CRITERIA.LONGEST)

                if True:
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
        if removed_gene[k]: continue
        genes_in_cluster[cluster_indexes[k]] += 1

    # calculating genes count which at least once overlapped to different gene
    for i in range(0, genes_cnt):
        if removed_gene[i]: continue
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
for i in range(1, 1000):
    if genes_by_clusters_length[i] > 0:
        print(str(i) + "-length clusters: " + str((genes_by_clusters_length[i] // i)))
print("")

print("--- %s seconds ---" % (time.time() - start_time))
