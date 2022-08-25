import numpy
from worker_genome import *
import scipy.stats

specie_list = [SPECIES.Homo_sapiens,
               SPECIES.Mus_musculus,
               SPECIES.Drosophila_melanogaster,
               SPECIES.Arabidopsis_thaliana
               ]


class Collected_Data:
    class Chromosome_Data:
        def __init__(self, label, genes_cnt):
            self.label = label
            self.genes_cnt = genes_cnt
            self.OGs_by_Reviewed_TYPE = {OVERLAP_TYPE.CONVERGENT: {},
                                         OVERLAP_TYPE.DIFF_NESTED_SHORTER: {},
                                         OVERLAP_TYPE.DIFF_NESTED_LONGER: {},
                                         OVERLAP_TYPE.DIVERGENT: {}}

    def __init__(self, genome: GenomeWorker):
        self.loaded_genes = 0
        self.filtered_genes = {"nested": 0, "novel": 0, "repeat": 0}

        self.median_iso = 0
        self.median_length = 0
        self.median_incomplete_utr5 = 0
        self.median_incomplete_utr3 = 0

        self.chrs_data = []
        for chr_index in range(1, genome.chromosomes_count() + 1):
            genes_cnt = genome.genes_count_on_chr(chr_index)
            chr_label = genome.chr_index2_seq_id(chr_index)
            self.chrs_data.append(Collected_Data.Chromosome_Data(chr_label, genes_cnt))

        self.pearson_tests = {OVERLAP_TYPE.CONVERGENT: (0, 0),
                              OVERLAP_TYPE.DIFF_NESTED_SHORTER: (0, 0),
                              OVERLAP_TYPE.DIFF_NESTED_LONGER: (0, 0),
                              OVERLAP_TYPE.DIVERGENT: (0, 0), }

        self.background_genes = []

        self.overlap_distances = {OVERLAP_TYPE.CONVERGENT: [], OVERLAP_TYPE.DIVERGENT: []}

    def get_OGs_by_type(self, ov_type: OVERLAP_TYPE, detailed_chr_data=False):
        label_arr = []
        size_arr = []
        cnt_arr = []

        gene_ids = []
        for chr_data in self.chrs_data:
            for gene_id, _ in chr_data.OGs_by_Reviewed_TYPE[ov_type].items():
                gene_ids.append(gene_id)

            label_arr.append(chr_data.label)
            size_arr.append(chr_data.genes_cnt)
            cnt_arr.append(len(chr_data.OGs_by_Reviewed_TYPE[ov_type]))
        return gene_ids if detailed_chr_data == False else (label_arr, size_arr, cnt_arr)

    def get_OR_union_sets(self):
        gene_set = {}
        for ov_type, _ in data.pearson_tests.items():
            for chr_data in self.chrs_data:
                for gene_id, _ in chr_data.OGs_by_Reviewed_TYPE[ov_type].items():
                    gene_set[gene_id] = True
        return list(gene_set.keys())

    def get_AND_union_sets(self, ov_type1, ov_type2):
        gene_set = {}
        final_gene_set = {}
        for chr_data in self.chrs_data:
            for gene_id, _ in chr_data.OGs_by_Reviewed_TYPE[ov_type1].items():
                gene_set[gene_id] = True
        for chr_data in self.chrs_data:
            for gene_id, _ in chr_data.OGs_by_Reviewed_TYPE[ov_type2].items():
                if gene_set.__contains__(gene_id):
                    final_gene_set[gene_id] = True
        return list(final_gene_set.keys())

    def calculate_filter_types(self, gen: GenomeWorker):
        ignores = gen.ignored_genes_by_types
        self.filtered_genes["nested"] = ignores["nested"] if ignores.__contains__("nested") else 0
        self.filtered_genes["novel"] = (ignores["is_readthrough"] if ignores.__contains__("is_readthrough") else 0) + \
                                       (ignores["is_pseudogene"] if ignores.__contains__("is_pseudogene") else 0) \
                                       + (ignores["is_novel"] if ignores.__contains__("is_novel") else 0) \
                                       + (ignores["is_predicted"] if ignores.__contains__("is_predicted") else 0)
        self.filtered_genes["repeat"] = (ignores["name_duplicated"] if ignores.__contains__("name_duplicated") else 0) \
                                        + (ignores["acc_duplicated"] if ignores.__contains__("acc_duplicated") else 0)


species_data = {}
max_sets_in_species = 0
for species in specie_list:
    genome = GenomeWorker(species, ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_FRAGMENTS, SEQUENCE_LOAD.NOT_LOAD)

    data = Collected_Data(genome)

    data.loaded_genes = genome.imported_protein_coding_genes
    data.calculate_filter_types(genome)

    isoforms_count = []
    incomplete_scores_utr5 = []
    incomplete_scores_utr3 = []
    genes_length = []

    for chr_index in range(1, genome.chromosomes_count() + 1):
        genes_cnt = genome.genes_count_on_chr(chr_index)
        for i in range(0, genes_cnt):
            gene1 = genome.gene_by_indexes(chr_index, i)

            data.background_genes.append(gene1.id)

            isoforms_count.append(len(genome.get_transcripts_from_gene(gene1.id)))
            genes_length.append(gene1.end - gene1.start + 1)
            incomplete_scores_utr5.append(genome.get_gene_incomplete_level(gene1.id, "utr5"))
            incomplete_scores_utr3.append(genome.get_gene_incomplete_level(gene1.id, "utr3"))

            for j in range(i + 1, genes_cnt):
                gene2 = genome.gene_by_indexes(chr_index, j)
                ov_type = genome.get_features_overlap_type(gene1, gene2, including_PO=True)
                if ov_type != OVERLAP_TYPE.NONE and gene1.strand != gene2.strand:
                    if ov_type == OVERLAP_TYPE.DIFF_NESTED:
                        if gene1.end - gene1.start < gene2.end - gene2.start:
                            data.chrs_data[chr_index - 1].OGs_by_Reviewed_TYPE[OVERLAP_TYPE.DIFF_NESTED_SHORTER][
                                gene1.id] = True
                            data.chrs_data[chr_index - 1].OGs_by_Reviewed_TYPE[OVERLAP_TYPE.DIFF_NESTED_LONGER][
                                gene2.id] = True
                        else:
                            data.chrs_data[chr_index - 1].OGs_by_Reviewed_TYPE[OVERLAP_TYPE.DIFF_NESTED_SHORTER][
                                gene2.id] = True
                            data.chrs_data[chr_index - 1].OGs_by_Reviewed_TYPE[OVERLAP_TYPE.DIFF_NESTED_LONGER][
                                gene1.id] = True
                    else:
                        data.chrs_data[chr_index - 1].OGs_by_Reviewed_TYPE[ov_type][
                            gene2.id] = True
                        data.chrs_data[chr_index - 1].OGs_by_Reviewed_TYPE[ov_type][
                            gene1.id] = True

                    if ov_type == OVERLAP_TYPE.CONVERGENT or ov_type == OVERLAP_TYPE.DIVERGENT:
                        data.overlap_distances[ov_type].append(
                            genome.get_features_overlap_length(gene1, gene2, including_PO=True))

    # Perform Pearson Correlation Tests
    for ov_type, _ in data.pearson_tests.items():
        labels, sizes, cnts = data.get_OGs_by_type(ov_type, detailed_chr_data=True)
        data.pearson_tests[ov_type] = scipy.stats.pearsonr(sizes, cnts)

    data.median_iso = numpy.median(isoforms_count)
    data.median_length = numpy.median(genes_length)
    data.median_incomplete_utr5 = numpy.average(incomplete_scores_utr5)
    data.median_incomplete_utr3 = numpy.average(incomplete_scores_utr3)

    species_data[species] = data

####################################################################################################

file = open("generated_data/general_stats.txt", 'w')

file.write("Filtered Genes in species (Nested,readthrough/predicted/novel/pseudo,repeat,)\n")
for species in specie_list:
    all_filtered = 0
    for _, cnt in species_data[species].filtered_genes.items():
        all_filtered += cnt
    file.write(f'{species.short_name()}\t{species_data[species].loaded_genes + all_filtered}\t{all_filtered}\t'
               f'{species_data[species].filtered_genes["nested"]}\t'
               f'{species_data[species].filtered_genes["novel"]}\t'
               f'{species_data[species].filtered_genes["repeat"]}\n')
file.write("\n")

file.write("General Data\n")
for species in specie_list:
    data = species_data[species]

    file.write(f'{species.short_name()}\t{data.loaded_genes}\t{data.median_iso}\t{data.median_length}\t'
               f'{data.median_incomplete_utr5}\t{data.median_incomplete_utr3}\t')

    for ov_type, (pearson_r, p_value) in data.pearson_tests.items():
        gene_ids = data.get_OGs_by_type(ov_type)
        ogs_percent = "{:.2f}".format(len(gene_ids) / data.loaded_genes * 100)
        beauty_result = "{:.2f}".format(pearson_r * pearson_r)
        beauty_p_value = "{:.2e}".format(p_value)
        file.write(f"{ogs_percent}%\t({beauty_result},{beauty_p_value})\t")
    file.write(f"\n")
file.write("\n")

file.write(f"Genes in Set (% Total Genes)\n")
all_gene_sets = {}
for species in specie_list:
    data = species_data[species]

    gene_sets = [data.get_OR_union_sets(),
                 data.get_OGs_by_type(OVERLAP_TYPE.CONVERGENT),
                 data.get_OGs_by_type(OVERLAP_TYPE.DIFF_NESTED_SHORTER),
                 data.get_OGs_by_type(OVERLAP_TYPE.DIFF_NESTED_LONGER),
                 data.get_OGs_by_type(OVERLAP_TYPE.DIVERGENT),
                 data.get_AND_union_sets(OVERLAP_TYPE.CONVERGENT, OVERLAP_TYPE.DIFF_NESTED_SHORTER),
                 data.get_AND_union_sets(OVERLAP_TYPE.CONVERGENT, OVERLAP_TYPE.DIFF_NESTED_LONGER),
                 data.get_AND_union_sets(OVERLAP_TYPE.CONVERGENT, OVERLAP_TYPE.DIVERGENT),
                 data.get_AND_union_sets(OVERLAP_TYPE.DIFF_NESTED_SHORTER, OVERLAP_TYPE.DIFF_NESTED_LONGER),
                 data.get_AND_union_sets(OVERLAP_TYPE.DIFF_NESTED_SHORTER, OVERLAP_TYPE.DIVERGENT),
                 data.get_AND_union_sets(OVERLAP_TYPE.DIFF_NESTED_LONGER, OVERLAP_TYPE.DIVERGENT),
                 ]

    all_gene_sets[species] = gene_sets

    file.write(f'{species.short_name()}\t')
    for gene_set in gene_sets:
        num = "{:.2f}".format(len(gene_set) / data.loaded_genes * 100)
        file.write(f'{num}%\t')
    file.write(f'\n')
file.write("\n")

file.close()

####################################################################################################
file = open("generated_data/OGs_by_types.txt", "w")

gene_row_index = 0
while True:
    flag = False
    for set_index in range(0, 5):
        for species in specie_list:
            if len(all_gene_sets[species][set_index]) <= gene_row_index:
                gene_id = ""
            else:
                flag = True
                gene_id = all_gene_sets[species][set_index][gene_row_index].replace("gene:", "")
            file.write(f"{gene_id}\t")
    file.write(f"\n")
    if not flag: break
    gene_row_index += 1

file.close()

####################################################################################################
file = open("generated_data/Background_genes.txt", "w")

gene_row_index = 0
while True:
    flag = False
    for species in specie_list:
        gene_ids = species_data[species].background_genes
        if len(gene_ids) <= gene_row_index:
            gene_id = ""
        else:
            flag = True
            gene_id = gene_ids[gene_row_index].replace("gene:", "")
        file.write(f"{gene_id}\t")
    file.write(f"\n")
    if not flag: break
    gene_row_index += 1

file.close()

####################################################################################################
file = open("generated_data/distances.txt", "w")

distance_row_index = 0
while True:
    flag = False
    for species in specie_list:
        data = species_data[species]
        ov_types = list(data.overlap_distances.keys())
        for ov_type in ov_types:
            distances = data.overlap_distances[ov_type]
            if len(distances) <= distance_row_index:
                distance = ""
            else:
                flag = True
                distance = distances[distance_row_index]
            file.write(f"{distance}\t")
    file.write(f"\n")
    if not flag: break
    distance_row_index += 1

file.close()
