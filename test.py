import statistics

import numpy
from worker_genome import *
import scipy.stats

specie_list = [
    SPECIES.Homo_sapiens,
    SPECIES.Mus_musculus,
    SPECIES.Drosophila_melanogaster,
    SPECIES.Arabidopsis_thaliana
]

file = open("nested_direction.txt", "w")
for species in specie_list:
    genome = GenomeWorker(species, ANNOTATION_LOAD.GENES, SEQUENCE_LOAD.NOT_LOAD)
    utr3 = []
    utr5 = []
    utr_3_genes = []
    utr_5_genes = []
    for chr_index in range(1, genome.chromosomes_count() + 1):
        genes_cnt = genome.genes_count_on_chr(chr_index)
        for i in range(0, genes_cnt):
            gene1 = genome.gene_by_indexes(chr_index, i)

            for j in range(i + 1, genes_cnt):
                gene2 = genome.gene_by_indexes(chr_index, j)
                ov_type = genome.get_features_overlap_type(gene1, gene2)
                if ov_type == OVERLAP_TYPE.DIFF_NESTED:
                    if gene1.end - gene1.start > gene2.end - gene2.start:
                        l1, r1, s1, g1, l2, r2, s2, g2 = gene1.start, gene1.end, gene1.strand, gene1, gene2.start, gene2.end, gene2.strand, gene2
                    else:
                        l2, r2, s2, g2, l1, r1, s1, g1 = gene1.start, gene1.end, gene1.strand, gene1, gene2.start, gene2.end, gene2.strand, gene2
                    l = l2 - l1
                    assert l2 >= l1
                    r = r1 - r2
                    assert r1 >= r2
                    if s1 == '+':
                        if l > r:
                            utr3.append(r)
                            utr_3_genes.append(g1)
                        else:
                            utr5.append(l)
                            utr_5_genes.append(g1)
                    else:
                        if l > r:
                            utr5.append(r)
                            utr_5_genes.append(g1)
                        else:
                            utr3.append(l)
                            utr_3_genes.append(g1)

    print(f"{len(utr5)}\t{len(utr3)}\t{statistics.median(utr5)}\t{statistics.median(utr3)}")

    file.write(f"{len(utr5)}\t{len(utr3)}\t{statistics.median(utr5)}\t{statistics.median(utr3)}\n")
    if species == SPECIES.Homo_sapiens:
        f = open("near_promotor.txt", "w")
        for gene in utr_5_genes:
            f.write(f"{genome.get_gene_symbol(gene)}\n")
        f.close()

        f = open("near_utr3.txt", "w")
        for gene in utr_3_genes:
            f.write(f"{genome.get_gene_symbol(gene)}\n")
        f.close()
