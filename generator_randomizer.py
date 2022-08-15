import numpy

from worker_genome import *
from GO_graph import *
import scipy.stats

specie_list = [SPECIES.Homo_sapiens,
               SPECIES.Mus_musculus,
               SPECIES.Drosophila_melanogaster,
               SPECIES.Arabidopsis_thaliana
               ]


def boundaries_standard(gene: Feature):
    return (gene.start, gene.end, gene.strand)


def boundaries_randomized(gene: Feature):
    locus = randomized_locus[gene.id]
    return locus[1], locus[2], locus[0]


total = {}

observed_OGs = {}
observed_OGs_conv = {}
observed_OGs_nested = {}
observed_OGs_div = {}
observed_OGs_div_large = {}
observed_OGs_prom = {}

expected_OGs = {}
expected_OGs_conv = {}
expected_OGs_nested = {}
expected_OGs_div = {}
expected_OGs_div_large = {}
expected_OGs_prom = {}

for species in specie_list:
    OGs = {}
    OGs_conv = {}
    OGs_nested = {}
    OGs_div = {}
    OGs_div_large = {}
    OGs_prom = {}
    genome = GenomeWorker(species, ANNOTATION_LOAD.GENES, SEQUENCE_LOAD.LOAD)

    total[species] = 0

    for chr_index in range(1, genome.chromosomes_count() + 1):
        genes_cnt = genome.genes_count_on_chr(chr_index)
        total[species] += genes_cnt

        for i in range(0, genes_cnt):
            gene1 = genome.gene_by_indexes(chr_index, i)

            for j in range(i + 1, genes_cnt):
                gene2 = genome.gene_by_indexes(chr_index, j)
                segment1, segment2 = boundaries_standard(gene1), boundaries_standard(gene2)
                ov_type = genome.get_segments_overlap_type(segment1, segment2, including_PO=True)
                if ov_type == OVERLAP_TYPE.CONVERGENT or ov_type == OVERLAP_TYPE.NEAR_DIVERGENT or ov_type == OVERLAP_TYPE.FAR_DIVERGENT\
                        or ov_type == OVERLAP_TYPE.PROMOTOR or ov_type == OVERLAP_TYPE.DIFF_NESTED:
                    OGs[gene1.id] = True
                    OGs[gene2.id] = True

                if ov_type == OVERLAP_TYPE.CONVERGENT:
                    OGs_conv[gene1.id] = True
                    OGs_conv[gene2.id] = True

                if ov_type == OVERLAP_TYPE.DIFF_NESTED:
                    OGs_nested[gene1.id] = True
                    OGs_nested[gene2.id] = True

                if ov_type == OVERLAP_TYPE.NEAR_DIVERGENT:
                        OGs_div[gene1.id] = True
                        OGs_div[gene2.id] = True
                if ov_type == OVERLAP_TYPE.FAR_DIVERGENT:
                        OGs_div_large[gene1.id] = True
                        OGs_div_large[gene2.id] = True
                if ov_type == OVERLAP_TYPE.PROMOTOR:
                    OGs_prom[gene1.id] = True
                    OGs_prom[gene2.id] = True

    observed_OGs[species] = len(OGs)
    observed_OGs_conv[species] = len(OGs_conv)
    observed_OGs_nested[species] = len(OGs_nested)
    observed_OGs_div[species] = len(OGs_div)
    observed_OGs_div_large[species] = len(OGs_div_large)
    observed_OGs_prom[species] = len(OGs_prom)

    expected_OGs[species] = []
    expected_OGs_conv[species] = []
    expected_OGs_nested[species] = []
    expected_OGs_div[species] = []
    expected_OGs_div_large[species] = []
    expected_OGs_prom[species] = []

    turns = 3

    while turns > 0:
        turns -= 1

        randomized_locus = {}
        for chr_index in range(1, genome.chromosomes_count() + 1):
            genes_cnt = genome.genes_count_on_chr(chr_index)
            chr_base_size = len(genome.retrieve_sequence_record(chr_index))
            for i in range(0, genes_cnt):
                gene1 = genome.gene_by_indexes(chr_index, i)
                size = gene1.end - gene1.start + 1
                strand = '+' if random.randint(0, 100) % 2 == 0 else '-'
                random_start = random.randint(0, chr_base_size - size)
                randomized_locus[gene1.id] = (strand, random_start, random_start + size - 1)

        OGs = {}
        OGs_conv = {}
        OGs_nested = {}
        OGs_div = {}
        OGs_div_large = {}
        OGs_prom = {}

        for chr_index in range(1, genome.chromosomes_count() + 1):
            genes_cnt = genome.genes_count_on_chr(chr_index)
            for i in range(0, genes_cnt):
                gene1 = genome.gene_by_indexes(chr_index, i)
                for j in range(i + 1, genes_cnt):
                    gene2 = genome.gene_by_indexes(chr_index, j)
                    segment1, segment2 = boundaries_randomized(gene1), boundaries_randomized(gene2)
                    ov_type = genome.get_segments_overlap_type(segment1, segment2, including_PO=True)
                    if ov_type == OVERLAP_TYPE.CONVERGENT or ov_type == OVERLAP_TYPE.FAR_DIVERGENT or ov_type == OVERLAP_TYPE.NEAR_DIVERGENT\
                            or ov_type == OVERLAP_TYPE.PROMOTOR or ov_type == OVERLAP_TYPE.DIFF_NESTED:
                        OGs[gene1.id] = True
                        OGs[gene2.id] = True

                    if ov_type == OVERLAP_TYPE.CONVERGENT:
                        OGs_conv[gene1.id] = True
                        OGs_conv[gene2.id] = True

                    if ov_type == OVERLAP_TYPE.DIFF_NESTED:
                        OGs_nested[gene1.id] = True
                        OGs_nested[gene2.id] = True

                    if ov_type == OVERLAP_TYPE.NEAR_DIVERGENT:
                        OGs_div[gene1.id] = True
                        OGs_div[gene2.id] = True

                    if ov_type == OVERLAP_TYPE.FAR_DIVERGENT:
                        OGs_div_large[gene1.id] = True
                        OGs_div_large[gene2.id] = True

                    if ov_type == OVERLAP_TYPE.PROMOTOR:
                        OGs_prom[gene1.id] = True
                        OGs_prom[gene2.id] = True

        expected_OGs[species].append(len(OGs))
        expected_OGs_conv[species].append(len(OGs_conv))
        expected_OGs_nested[species].append(len(OGs_nested))
        expected_OGs_div[species].append(len(OGs_div))
        expected_OGs_div_large[species].append(len(OGs_div_large))
        expected_OGs_prom[species].append(len(OGs_prom))

file = open("generated_data/Expected_OGs_counts.txt", "w")

file.write("Expected OGs:\n")
for species in specie_list:
    file.write(f"{species.short_name()}\t{total[species]}\t{observed_OGs[species]}")
    for i in range(len(expected_OGs[species])):
        file.write(f"\t{expected_OGs[species][i]}")
    file.write("\n")

file.write("\nExpected Convergent Count:\n")
for species in specie_list:
    file.write(f"{species.short_name()}\t{total[species]}\t{observed_OGs_conv[species]}")
    for i in range(len(expected_OGs_conv[species])):
        file.write(f"\t{expected_OGs_conv[species][i]}")
    file.write("\n")

file.write("\nExpected Diff Nested Count:\n")
for species in specie_list:
    file.write(f"{species.short_name()}\t{total[species]}\t{observed_OGs_nested[species]}")
    for i in range(len(expected_OGs_nested[species])):
        file.write(f"\t{expected_OGs_nested[species][i]}")
    file.write("\n")

file.write("\nExpected Divergent Count:\n")
for species in specie_list:
    file.write(f"{species.short_name()}\t{total[species]}\t{observed_OGs_div[species]}")
    for i in range(len(expected_OGs_div[species])):
        file.write(f"\t{expected_OGs_div[species][i]}")
    file.write("\n")

file.write("\nExpected Divergent (Large) Count:\n")
for species in specie_list:
    file.write(f"{species.short_name()}\t{total[species]}\t{observed_OGs_div_large[species]}")
    for i in range(len(expected_OGs_div_large[species])):
        file.write(f"\t{expected_OGs_div_large[species][i]}")
    file.write("\n")

file.write("\nExpected Promotor sharing Count:\n")
for species in specie_list:
    file.write(f"{species.short_name()}\t{total[species]}\t{observed_OGs_prom[species]}")
    for i in range(len(expected_OGs_prom[species])):
        file.write(f"\t{expected_OGs_prom[species][i]}")
    file.write("\n")

file.close()
