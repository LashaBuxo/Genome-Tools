import numpy

from worker_genome import *
from GO_graph import *
import scipy.stats
from AgeData import *

specie_list = [SPECIES.Homo_sapiens,
               SPECIES.Mus_musculus,
               SPECIES.Drosophila_melanogaster,
               SPECIES.Arabidopsis_thaliana
               ]

ov_types = [OVERLAP_TYPE.CONVERGENT,
            OVERLAP_TYPE.DIFF_NESTED,
            OVERLAP_TYPE.FAR_DIVERGENT,
            OVERLAP_TYPE.NEAR_DIVERGENT,
            OVERLAP_TYPE.PROMOTOR]

for species in specie_list:
    # Load necessary data
    genome = GenomeWorker(species, ANNOTATION_LOAD.GENES, SEQUENCE_LOAD.NOT_LOAD)
    age_data = AgeData(species, genome)

    OGs_by_type = {OVERLAP_TYPE.CONVERGENT: {}, OVERLAP_TYPE.DIFF_NESTED: {},
                   OVERLAP_TYPE.FAR_DIVERGENT: {}, OVERLAP_TYPE.NEAR_DIVERGENT: {},
                   OVERLAP_TYPE.PROMOTOR: {}}

    for chr_index in range(1, genome.chromosomes_count() + 1):
        genes_cnt = genome.genes_count_on_chr(chr_index)
        for i in range(0, genes_cnt):
            gene1 = genome.gene_by_indexes(chr_index, i)
            for j in range(i + 1, genes_cnt):
                gene2 = genome.gene_by_indexes(chr_index, j)
                ov_type = genome.get_features_overlap_type(gene1, gene2, including_PO=True)
                if ov_type==OVERLAP_TYPE.CONVERGENT or ov_type==OVERLAP_TYPE.FAR_DIVERGENT\
                        or ov_type==OVERLAP_TYPE.NEAR_DIVERGENT or ov_type==OVERLAP_TYPE.DIFF_NESTED \
                        or ov_type==OVERLAP_TYPE.PROMOTOR:
                    OGs_by_type[ov_type][gene1.id] = True
                    OGs_by_type[ov_type][gene2.id] = True

    outfile = open(f"./generated_data/age_vs_OGs_data_{species.short_name()}.txt", "w")

    age_group_type = AgeData.AgeGroupsType.annotated_groups
    groups = age_data.calculate_groups(age_group_type)
    for group_index in range(len(groups)):
        group = groups[group_index]
        group_label = age_data.get_group_label(group, group_index, age_group_type)
        gene_ids = age_data.get_gene_ids_in_group(group)

        extra_data = ""
        for ov_type in ov_types:
            overlapped_ids = OGs_by_type[ov_type]
            ov_count = 0
            for gene_id in gene_ids:
                if overlapped_ids.__contains__(gene_id):
                    ov_count += 1
            extra_data += f"\t{ov_count}"
        outfile.write(f"{group_label}\t{len(gene_ids)}{extra_data}\n")
    outfile.close()
