from GO_graph import *

from worker_genome import *
graph=GOGraph()
genome = GenomeWorker(SPECIES.Homo_sapiens, ANNOTATION_LOAD.GENES, SEQUENCE_LOAD.NOT_LOAD)

genes=graph.get_term_gene_list(genome,"GO:0006281")
nested_longer={}
nested_longer_dna={}
total_genes=0
for chr_index in range(1, genome.chromosomes_count() + 1):
    genes_cnt = genome.genes_count_on_chr(chr_index)
    total_genes+=genes_cnt
    for i in range(0, genes_cnt):
        gene1 = genome.gene_by_indexes(chr_index, i)

        for j in range(i + 1, genes_cnt):
            gene2 = genome.gene_by_indexes(chr_index, j)
            ov_type = genome.get_features_overlap_type(gene1, gene2)
            if ov_type == OVERLAP_TYPE.DIFF_NESTED:
                if gene1.end-gene1.start<gene2.end-gene2.start:
                    nested_longer[genome.get_gene_symbol(gene2)]=True
                    if genes.__contains__(genome.get_gene_symbol(gene2)):
                        nested_longer_dna[genome.get_gene_symbol(gene2)]=True
                else:
                    nested_longer[genome.get_gene_symbol(gene1)] = True
                    if genes.__contains__(genome.get_gene_symbol(gene1)):
                        nested_longer_dna[genome.get_gene_symbol(gene1)] = True

print(f"From {len(nested_longer)}/{total_genes} Nested (diff, longer) {len(nested_longer_dna)}/{len(genes)} genes are DNA repairr genes")

for sym,_ in nested_longer_dna.items():
    print(sym)

