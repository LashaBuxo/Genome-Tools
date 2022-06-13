from worker_genome import *
from worker_analyzer import *

genome = GenomeWorker(SPECIES.Homo_sapiens, ANNOTATIONS.ENSEMBL,
                      ANNOTATION_LOAD.GENES, SEQUENCE_LOAD.NOT_LOAD)

total_genes = 0
genes_greater_1000000 = 0

for chr_id in range(1, genome.chromosomes_count() + 1):
    genes_cnt = genome.genes_count_on_chr(chr_id)
    total_genes += genes_cnt
    for i in range(0, genes_cnt):
        gene = genome.gene_by_indexes(chr_id, i)
        genes_greater_1000000 += 1 if gene.end - gene.start + 1 > 10000 else 0
print(f'{genes_greater_1000000} genes from {total_genes}')