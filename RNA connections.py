import statistics
from GO_graph import *

reference_species = SPECIES.Homo_sapiens
genome = GenomeWorker(reference_species, ANNOTATION_LOAD.GENES, SEQUENCE_LOAD.NOT_LOAD, load_MitoCARTA=True)

rna_genes = genome.get_RNA_genes()

OGs = {}
for chr_index in range(1, genome.chromosomes_count() + 1):
    genes_cnt = genome.genes_count_on_chr(chr_index)
    for i in range(0, genes_cnt):
        gene1 = genome.gene_by_indexes(chr_index, i)
        for j in range(i + 1, genes_cnt):
            gene2 = genome.gene_by_indexes(chr_index, j)
            if genome.get_overlap_type(gene1, gene2) != OVERLAP_TYPE.NONE:
                OGs[gene1.id] = True
                OGs[gene2.id] = True

Genes_overlapped_to_RNA = 0
Genes_cnt = 0

OGs_overlapped_to_RNA = 0
OGs_cnt = 0

Mito_genes_overlapped_to_RNA = 0
Mito_genes_cnt = 0

Mito_OGs_overlapped_to_RNA = 0
Mito_OGs_cnt = 0

GO_graph = GOGraph()
mito_genes = GO_graph.get_term_gene_list(genome, "GO:0005739")
dna_genes = GO_graph.get_term_gene_list(genome, "GO:0006281")


# print(specific_genes)

def is_gene_overlapped_to_RNA(gene: Feature):
    for rna_gene in rna_genes:
        if genome.get_overlap_type(rna_gene, gene) != OVERLAP_TYPE.NONE:
            print(rna_gene.id)
            return True
    return False


print(is_gene_overlapped_to_RNA(genome.feature_by_id("gene:ENSG00000239620")))

for chr_index in range(1, genome.chromosomes_count() + 1):
    genes_cnt = genome.genes_count_on_chr(chr_index)
    for i in range(0, genes_cnt):
        gene1 = genome.gene_by_indexes(chr_index, i)
        sym = genome.get_gene_symbol(gene1)

        if is_gene_overlapped_to_RNA(gene1):
            Genes_overlapped_to_RNA += 1
        Genes_cnt += 1

        if OGs.__contains__(gene1.id):
            OGs_cnt += 1
            if is_gene_overlapped_to_RNA(gene1):
                OGs_overlapped_to_RNA += 1

        if genome.is_gene_MITO(gene1.id):
            Mito_genes_cnt += 1
            if is_gene_overlapped_to_RNA(gene1):
                Mito_genes_overlapped_to_RNA += 1

            if OGs.__contains__(gene1.id):
                Mito_OGs_cnt += 1
                if is_gene_overlapped_to_RNA(gene1):
                    Mito_OGs_overlapped_to_RNA += 1

print(f"Genes overlapped (total = {Genes_cnt}) to at least 1 RNA gene: {Genes_overlapped_to_RNA}")
print(f"OGs overlapped (total = {OGs_cnt}) to at least 1 RNA gene: {OGs_overlapped_to_RNA}")
print(f"Mito genes overlapped (total = {Mito_genes_cnt}) to at least 1 RNA gene: {Mito_genes_overlapped_to_RNA}")
print(f"Mito OGs overlapped (total = {Mito_OGs_cnt}) to at least 1 RNA gene: {Mito_OGs_overlapped_to_RNA}")
