from worker_genome import *

file = open("used_data/viability.csv", "r")
lines = file.readlines()
lethal_genes = {}
for index in range(1, len(lines)):
    arr = lines[index].replace('\n', '').replace('"', '').split(',')
    gene_sym, status = arr[0], arr[len(arr) - 4]
    if status == "viable":
        lethal_genes[gene_sym.upper()] = True
    #print(f"{gene_sym}\t{status}")

genome1 = GenomeWorker(SPECIES.Mus_musculus, ANNOTATION_LOAD.GENES, SEQUENCE_LOAD.NOT_LOAD)

ov_counts = {OVERLAP_TYPE.CONVERGENT: {}, OVERLAP_TYPE.DIVERGENT: {}, OVERLAP_TYPE.DIFF_NESTED_LONGER: {},
             OVERLAP_TYPE.DIFF_NESTED_SHORTER: {}}

for chr_index in range(1, genome1.chromosomes_count() + 1):
    genes_cnt = genome1.genes_count_on_chr(chr_index)
    for i in range(0, genes_cnt):
        gene1 = genome1.gene_by_indexes(chr_index, i)
        for j in range(i + 1, genes_cnt):
            gene2 = genome1.gene_by_indexes(chr_index, j)
            ov_type = genome1.get_features_overlap_type(gene1, gene2, including_PO=True)
            if gene1.strand == gene2.strand or ov_type==OVERLAP_TYPE.NONE: continue

            if ov_type == OVERLAP_TYPE.DIFF_NESTED:
                if gene1.end - gene1.start < gene2.end - gene2.start:
                    ov_counts[OVERLAP_TYPE.DIFF_NESTED_SHORTER][gene1.id] = True
                    ov_counts[OVERLAP_TYPE.DIFF_NESTED_LONGER][gene2.id] = True
                else:
                    ov_counts[OVERLAP_TYPE.DIFF_NESTED_SHORTER][gene2.id] = True
                    ov_counts[OVERLAP_TYPE.DIFF_NESTED_LONGER][gene1.id] = True
            else:
                ov_counts[ov_type][gene1.id] = True
                ov_counts[ov_type][gene2.id] = True
print(len(lethal_genes))
for ov_type in list(ov_counts.keys()):
    lethals = 0
    for gene_id, _ in ov_counts[ov_type].items():
        gene = genome1.feature_by_id(gene_id)
        sym = genome1.get_gene_symbol(gene)
        if lethal_genes.__contains__(sym):
            lethals += 1
    print(f"{ov_type.short_name()}\t{len(ov_counts[ov_type])}\t{lethals}\t")
