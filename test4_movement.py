import numpy
from worker_genome import *
import scipy.stats

specie_list = [SPECIES.Mus_musculus, SPECIES.Drosophila_melanogaster, ]

genome1 = GenomeWorker(specie_list[0], ANNOTATION_LOAD.GENES, SEQUENCE_LOAD.NOT_LOAD)
genome2 = GenomeWorker(specie_list[1], ANNOTATION_LOAD.GENES, SEQUENCE_LOAD.NOT_LOAD)

orth_data = {}

file = open("used_data/orthologous_mouse_drosophila.txt.txt", "r")
lines = file.readlines()
for index in range(1, len(lines)):
    arr = lines[index].replace('\n', '').split(',')
    species1_gene_id = arr[0]
    species2_gene_id = arr[2]

    if genome1.feature_by_id(f"gene:{species1_gene_id}") is None or genome2.feature_by_id(
            f"gene:{species2_gene_id}") is None:
        continue

    if not orth_data.__contains__(species1_gene_id):
        orth_data[species1_gene_id] = []
    if not orth_data.__contains__(species2_gene_id):
        orth_data[species2_gene_id] = []
    orth_data[species1_gene_id].append(species2_gene_id)
    orth_data[species2_gene_id].append(species1_gene_id)

ov_types = [OVERLAP_TYPE.DIVERGENT, OVERLAP_TYPE.DIFF_NESTED_SHORTER, OVERLAP_TYPE.DIFF_NESTED_LONGER,
            OVERLAP_TYPE.CONVERGENT]

species1_OGs = {OVERLAP_TYPE.DIVERGENT: {}, OVERLAP_TYPE.DIFF_NESTED_SHORTER: {}, OVERLAP_TYPE.DIFF_NESTED_LONGER: {},
                OVERLAP_TYPE.CONVERGENT: {}}

species2_OGs = {OVERLAP_TYPE.DIVERGENT: {}, OVERLAP_TYPE.DIFF_NESTED_SHORTER: {}, OVERLAP_TYPE.DIFF_NESTED_LONGER: {},
                OVERLAP_TYPE.CONVERGENT: {}}

for chr_index in range(1, genome1.chromosomes_count() + 1):
    genes_cnt = genome1.genes_count_on_chr(chr_index)
    for i in range(0, genes_cnt):
        gene1 = genome1.gene_by_indexes(chr_index, i)
        for j in range(i + 1, genes_cnt):
            gene2 = genome1.gene_by_indexes(chr_index, j)
            ov_type = genome1.get_features_overlap_type(gene1, gene2)
            if ov_type == OVERLAP_TYPE.DIFF_NESTED:
                if gene1.end - gene1.start > gene2.end - gene2.start:
                    species1_OGs[OVERLAP_TYPE.DIFF_NESTED_LONGER][gene1.id.replace("gene:", "")] = True
                    species1_OGs[OVERLAP_TYPE.DIFF_NESTED_SHORTER][gene2.id.replace("gene:", "")] = True
                else:
                    species1_OGs[OVERLAP_TYPE.DIFF_NESTED_LONGER][gene2.id.replace("gene:", "")] = True
                    species1_OGs[OVERLAP_TYPE.DIFF_NESTED_SHORTER][gene1.id.replace("gene:", "")] = True
            elif ov_type != OVERLAP_TYPE.NONE and gene1.strand != gene2.strand:
                species1_OGs[ov_type][gene1.id.replace("gene:", "")] = True
                species1_OGs[ov_type][gene2.id.replace("gene:", "")] = True

for chr_index in range(1, genome2.chromosomes_count() + 1):
    genes_cnt = genome2.genes_count_on_chr(chr_index)
    for i in range(0, genes_cnt):
        gene1 = genome2.gene_by_indexes(chr_index, i)
        for j in range(i + 1, genes_cnt):
            gene2 = genome2.gene_by_indexes(chr_index, j)
            ov_type = genome2.get_features_overlap_type(gene1, gene2)
            if ov_type == OVERLAP_TYPE.DIFF_NESTED:
                if gene1.end - gene1.start > gene2.end - gene2.start:
                    species2_OGs[OVERLAP_TYPE.DIFF_NESTED_LONGER][gene1.id.replace("gene:", "")] = True
                    species2_OGs[OVERLAP_TYPE.DIFF_NESTED_SHORTER][gene2.id.replace("gene:", "")] = True
                else:
                    species2_OGs[OVERLAP_TYPE.DIFF_NESTED_LONGER][gene2.id.replace("gene:", "")] = True
                    species2_OGs[OVERLAP_TYPE.DIFF_NESTED_SHORTER][gene1.id.replace("gene:", "")] = True
            elif ov_type != OVERLAP_TYPE.NONE and gene1.strand != gene2.strand:
                species2_OGs[ov_type][gene1.id.replace("gene:", "")] = True
                species2_OGs[ov_type][gene2.id.replace("gene:", "")] = True

results = open("movement.txt", "w")
for ov_type in ov_types:
    not_conserved, conserved, conserved_and_1to1 = 0, 0, 0
    counts_by_overlap_type = {OVERLAP_TYPE.DIVERGENT: 0,
                              OVERLAP_TYPE.DIFF_NESTED_SHORTER: 0,
                              OVERLAP_TYPE.DIFF_NESTED_LONGER: 0,
                              OVERLAP_TYPE.CONVERGENT: 0}
    no_overlap = 0

    evidence_genes = {
        OVERLAP_TYPE.DIVERGENT: [],
        OVERLAP_TYPE.DIFF_NESTED_SHORTER: [],
        OVERLAP_TYPE.DIFF_NESTED_LONGER: [],
        OVERLAP_TYPE.CONVERGENT: [], }

    for gene1_id, _ in species1_OGs[ov_type].items():
        if not orth_data.__contains__(gene1_id):
            not_conserved += 1
            continue
        conserved += 1
        if len(orth_data[gene1_id]) != 1:
            continue

        conserved_and_1to1 += 1
        other_gene1_id = orth_data[gene1_id][0]
        other_gene1 = genome2.feature_by_id(f"gene:{other_gene1_id}")

        flag = False
        for new_ov_type in ov_types:
            if species2_OGs[new_ov_type].__contains__(other_gene1_id):
                counts_by_overlap_type[new_ov_type] += 1
                flag = True
        if flag == False:
            no_overlap += 1
    if conserved_and_1to1==0:
        conserved_and_1to1=1
    results.write(f"{ov_type.short_name()}\t{len(species1_OGs[ov_type])}\t"
                  f"{not_conserved}\t{conserved}\t{conserved_and_1to1}\t"
                  f"{'{:.2f}'.format(counts_by_overlap_type[OVERLAP_TYPE.DIVERGENT] / conserved_and_1to1 * 100)}%\t"
                  f"{'{:.2f}'.format(counts_by_overlap_type[OVERLAP_TYPE.DIFF_NESTED_SHORTER] / conserved_and_1to1 * 100)}%\t"
                  f"{'{:.2f}'.format(counts_by_overlap_type[OVERLAP_TYPE.DIFF_NESTED_LONGER] / conserved_and_1to1 * 100)}%\t"
                  f"{'{:.2f}'.format(counts_by_overlap_type[OVERLAP_TYPE.CONVERGENT] / conserved_and_1to1 * 100)}%\t"
                  f"{'{:.2f}'.format(no_overlap / conserved_and_1to1 * 100)}%\t\n")
results.close()
