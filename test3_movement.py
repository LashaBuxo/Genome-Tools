import numpy
from worker_genome import *
import scipy.stats

specie_list = [SPECIES.Homo_sapiens, SPECIES.Mus_musculus, SPECIES.Drosophila_melanogaster]

genome1 = GenomeWorker(specie_list[0], ANNOTATION_LOAD.GENES, SEQUENCE_LOAD.NOT_LOAD)
genome2 = GenomeWorker(specie_list[1], ANNOTATION_LOAD.GENES, SEQUENCE_LOAD.NOT_LOAD)
genome3 = GenomeWorker(specie_list[2], ANNOTATION_LOAD.GENES, SEQUENCE_LOAD.NOT_LOAD)

human_mouse_orth_data = {}
file = open("used_data/orthologous_human_mouse.txt", "r")
lines = file.readlines()
for index in range(1, len(lines)):
    arr = lines[index].replace('\n', '').split(',')
    human_gene_id = arr[0]
    mouse_gene_id = arr[2]

    if genome1.feature_by_id(f"gene:{human_gene_id}") is None or genome2.feature_by_id(f"gene:{mouse_gene_id}") is None:
        continue

    if not human_mouse_orth_data.__contains__(human_gene_id):
        human_mouse_orth_data[human_gene_id] = []
    if not human_mouse_orth_data.__contains__(mouse_gene_id):
        human_mouse_orth_data[mouse_gene_id] = []
    human_mouse_orth_data[human_gene_id].append(mouse_gene_id)
    human_mouse_orth_data[mouse_gene_id].append(human_gene_id)

human_droso_orth_data = {}
file = open("used_data/orthologous_human_drosophila.txt", "r")
lines = file.readlines()
for index in range(1, len(lines)):
    arr = lines[index].replace('\n', '').split(',')
    human_gene_id = arr[0]
    mouse_gene_id = arr[2]

    if genome1.feature_by_id(f"gene:{human_gene_id}") is None or genome3.feature_by_id(f"gene:{mouse_gene_id}") is None:
        continue

    if not human_droso_orth_data.__contains__(human_gene_id):
        human_droso_orth_data[human_gene_id] = []
    if not human_droso_orth_data.__contains__(mouse_gene_id):
        human_droso_orth_data[mouse_gene_id] = []
    human_droso_orth_data[human_gene_id].append(mouse_gene_id)
    human_droso_orth_data[mouse_gene_id].append(human_gene_id)

mouse_droso_orth_data = {}
file = open("used_data/orthologous_mouse_drosophila.txt.txt", "r")
lines = file.readlines()
for index in range(1, len(lines)):
    arr = lines[index].replace('\n', '').split(',')
    human_gene_id = arr[0]
    mouse_gene_id = arr[2]

    if genome2.feature_by_id(f"gene:{human_gene_id}") is None or genome3.feature_by_id(f"gene:{mouse_gene_id}") is None:
        continue

    if not mouse_droso_orth_data.__contains__(human_gene_id):
        mouse_droso_orth_data[human_gene_id] = []
    if not mouse_droso_orth_data.__contains__(mouse_gene_id):
        mouse_droso_orth_data[mouse_gene_id] = []
    mouse_droso_orth_data[human_gene_id].append(mouse_gene_id)
    mouse_droso_orth_data[mouse_gene_id].append(human_gene_id)

print(len(human_mouse_orth_data))
print(len(human_droso_orth_data))
print(len(mouse_droso_orth_data))
cnt_mouse_droso = 0
cnt_not_found = 0
for chr_index in range(1, genome1.chromosomes_count() + 1):
    genes_cnt = genome1.genes_count_on_chr(chr_index)
    for i in range(0, genes_cnt):
        gene1 = genome1.gene_by_indexes(chr_index, i)
        gene1_id = gene1.id.replace('gene:', '')

        if not human_mouse_orth_data.__contains__(gene1_id): continue
        if not human_droso_orth_data.__contains__(gene1_id): continue
        if len(human_mouse_orth_data[gene1_id]) != 1 or len(human_droso_orth_data[gene1_id]) != 1:
            continue
        mouse_orth_id = human_mouse_orth_data[gene1_id][0]
        droso_orth_id = human_droso_orth_data[gene1_id][0]

        if mouse_droso_orth_data.__contains__(mouse_orth_id) and \
                mouse_droso_orth_data[mouse_orth_id].__contains__(droso_orth_id):
            cnt_mouse_droso += 1
        else:
            cnt_not_found += 1
            print(gene1_id)
            print(mouse_orth_id)
            print(droso_orth_id)
print(f"found: {cnt_mouse_droso} not found: {cnt_not_found}")
