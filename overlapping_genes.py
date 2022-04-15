import genome_tools as genome
import time

start_time = time.time()

total_genes = 0
positive_genes = 0
negative_genes = 0
og_count = 0
og_pairs_count = 0
og_clusters_count = 0

for chr_id in range(1, 26):
    genome.preprocess_annotation_for_chr(chr_id)
    genes_cnt = genome.genes_count_on_chr(chr_id)

    # no same id over genes
    for i in range(0, genes_cnt):
        for j in range(i + 1, genes_cnt):
            if genome.gene_by_ind(chr_id, i).id == genome.gene_by_ind(chr_id, j).id:
                print("lasha")

    total_genes += genes_cnt
    for i in range(0, genes_cnt):
        gene = genome.gene_by_ind(chr_id, i)
        if gene.strand == '+':
            positive_genes = positive_genes + 1
        if gene.strand == '-':
            negative_genes = negative_genes + 1

    cluster_indexes = [0] * genes_cnt
    for i in range(0, genes_cnt):
        cluster_indexes[i] = i

    for i in range(0, genes_cnt):
        for j in range(i + 1, genes_cnt):
            if genome.are_genes_overlapped(genome.gene_by_ind(chr_id, i), genome.gene_by_ind(chr_id, j)):
                old_cluster_index = cluster_indexes[i]
                new_cluster_index = cluster_indexes[j]
                for k in range(0, genes_cnt):
                    if cluster_indexes[k] == old_cluster_index:
                        cluster_indexes[k] = new_cluster_index
                og_pairs_count += 1

    genes_in_cluster = [0] * genes_cnt
    for k in range(0, genes_cnt):
        genes_in_cluster[cluster_indexes[k]] += 1

    for i in range(0, genes_cnt):
        if genes_in_cluster[cluster_indexes[i]] > 1:
            og_count += 1
        if genes_in_cluster[i] > 1:
            og_clusters_count += 1

print("Number of genes: " + str(total_genes))
print("Number of genes on Positive(+) Strand: " + str(positive_genes))
print("Number of genes on Negative(-) Strand: " + str(negative_genes))
print("Number of OG: " + str(og_count))
print("Number of OG pairs: " + str(og_pairs_count))
print("Number of OG clusters: " + str(og_clusters_count))

print("--- %s seconds ---" % (time.time() - start_time))
