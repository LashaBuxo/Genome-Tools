import numpy
from worker_genome import *
import scipy.stats

specie_list = [SPECIES.Homo_sapiens, SPECIES.Mus_musculus, ]

genome1 = GenomeWorker(specie_list[0], ANNOTATION_LOAD.GENES, SEQUENCE_LOAD.NOT_LOAD)
genome2 = GenomeWorker(specie_list[1], ANNOTATION_LOAD.GENES, SEQUENCE_LOAD.NOT_LOAD)

orth_data = {}

file = open("used_data/orthologous_human_mouse.txt", "r")
lines = file.readlines()
for index in range(1, len(lines)):
    arr = lines[index].replace('\n', '').split(',')
    human_gene_id = arr[0]
    mouse_gene_id = arr[2]

    if genome1.feature_by_id(f"gene:{human_gene_id}") is None or genome2.feature_by_id(f"gene:{mouse_gene_id}") is None:
        continue

    if not orth_data.__contains__(human_gene_id):
        orth_data[human_gene_id] = []
    if not orth_data.__contains__(mouse_gene_id):
        orth_data[mouse_gene_id] = []
    orth_data[human_gene_id].append(mouse_gene_id)
    orth_data[mouse_gene_id].append(human_gene_id)

main_nested_gene_pairs = []
for chr_index in range(1, genome1.chromosomes_count() + 1):
    genes_cnt = genome1.genes_count_on_chr(chr_index)
    for i in range(0, genes_cnt):
        gene1 = genome1.gene_by_indexes(chr_index, i)
        for j in range(i + 1, genes_cnt):
            gene2 = genome1.gene_by_indexes(chr_index, j)
            ov_type = genome1.get_features_overlap_type(gene1, gene2)
            if ov_type == OVERLAP_TYPE.DIFF_NESTED:
                if gene1.end - gene1.start > gene2.end - gene2.start:
                    main_nested_gene_pairs.append((gene1.id.replace("gene:", ""), gene2.id.replace("gene:", "")))
                else:
                    main_nested_gene_pairs.append((gene2.id.replace("gene:", ""), gene1.id.replace("gene:", "")))

none_conserved, one_conserved, one_and_longer_conserved, both_conserved, both_conserved_only_1to1, both_conserved_only1to1_and_same = 0, 0,0, 0, 0, 0
counts_by_overlap_type = {OVERLAP_TYPE.DIVERGENT: 0,
                          OVERLAP_TYPE.DIFF_NESTED: 0,
                          OVERLAP_TYPE.CONVERGENT: 0}
no_overlap_and_same_chr = 0
no_overlap_and_diff_chr = 0

evidence_genes = {
    OVERLAP_TYPE.DIVERGENT: [],
    OVERLAP_TYPE.DIFF_NESTED: [],
    OVERLAP_TYPE.CONVERGENT: [], }

for gene1_id, gene2_id in main_nested_gene_pairs:
    if not orth_data.__contains__(gene1_id) and not orth_data.__contains__(gene2_id):
        none_conserved += 1
        continue

    if orth_data.__contains__(gene1_id) and orth_data.__contains__(gene2_id):
        both_conserved += 1

        if len(orth_data[gene1_id]) != 1 or len(orth_data[gene2_id]) != 1:
            continue

        both_conserved_only_1to1 += 1
        other_gene1_id = orth_data[gene1_id][0]
        other_gene2_id = orth_data[gene2_id][0]

        if other_gene1_id == other_gene2_id:
            both_conserved_only1to1_and_same += 1
            continue

        other_gene1 = genome2.feature_by_id(f"gene:{other_gene1_id}")
        other_gene2 = genome2.feature_by_id(f"gene:{other_gene2_id}")

        ov_type = genome2.get_features_overlap_type(other_gene1, other_gene2, including_PO=True)
        if ov_type != OVERLAP_TYPE.NONE:
            assert counts_by_overlap_type.__contains__(ov_type)
            counts_by_overlap_type[ov_type] += 1
            evidence_genes[ov_type].append(((gene1_id, gene2_id), (other_gene1_id, other_gene2_id)))
        else:
            if other_gene1.chrom == other_gene2.chrom:
                no_overlap_and_same_chr += 1
            else:
                no_overlap_and_diff_chr += 1
        continue

    if orth_data.__contains__(gene1_id) or orth_data.__contains__(gene2_id):
        one_conserved += 1
        if orth_data.__contains__(gene1_id):
            one_and_longer_conserved += 1

results = open("movement.txt", "w")

results.write(f"{len(main_nested_gene_pairs)}\t"
              f"{none_conserved}\t{one_conserved}\t{one_and_longer_conserved}\t{both_conserved}\t"
              f"{both_conserved_only_1to1}\t{both_conserved_only1to1_and_same}\t"
              f"{counts_by_overlap_type[OVERLAP_TYPE.DIVERGENT]}\t"
              f"{counts_by_overlap_type[OVERLAP_TYPE.DIFF_NESTED]}\t"
              f"{counts_by_overlap_type[OVERLAP_TYPE.CONVERGENT]}\t"
              f"{no_overlap_and_same_chr}\t"
              f"{no_overlap_and_diff_chr}\t")

for ov_type in list(evidence_genes.keys()):
    results.write(f"\n{ov_type.short_name()}\n")
    for (main_gene1, main_gene2), (other_gene1, other_gene2) in evidence_genes[ov_type]:
        results.write(f"{main_gene1}\t{main_gene2}\t{other_gene1}\t{other_gene2}\t\n")

results.close()
