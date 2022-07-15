import numpy

from worker_genome import *
from GO_graph import *
import scipy.stats
from analyzer_graph import *


class AgeData:
    class INFO:
        def __init__(self, gene_id, gene_age, group_label):
            self.gene_id = gene_id
            self.gene_age = gene_age
            self.group_label = group_label

        def get_gene_age_numeric(self):
            if self.gene_age.__contains__(">"):
                return int(self.gene_age.replace(">", ""))
            return int(self.gene_age)

    def __init__(self, spec: SPECIES, gen: GenomeWorker):
        self.spec = spec
        self.gen = gen

        self.data_by_gene = {}

        self.load_age_data()

    def load_age_data(self):
        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5174733/?fbclid=IwAR2JHOlsC36nJGXCyiY_8MNWQ5xszKf4BGBVBSbeZoW1sz8bJPT3TNyGc7I
        # ages_file = open(f"./used_data/age_data/ages.txt", "r")
        # lines = ages_file.readlines()
        #
        # for age_line in lines:
        #     if age_line.__contains__("ensembl_id"): continue
        #     arr = age_line.replace("\n", '').split("\t")
        #     gene_id, group_index = arr[0], arr[1]
        #
        #     # we don't need gene age data, if gene is not loaded in genome
        #     if self.gen.feature_by_id(f"gene:{gene_id}") is None: continue
        #
        #     data = self.INFO(gene_id, gene_age=group_index, group_label=f"{group_index}")
        #     assert not self.data_by_gene.__contains__(gene_id)
        #     self.data_by_gene[gene_id] = data

        # http://genorigin.chenzxlab.cn/
        ages_file = open(f"./used_data/age_data/{self.spec.name}.csv", "r")
        lines = ages_file.readlines()

        for age_line in lines:
            if age_line.__contains__("ensembl_id"): continue
            arr = age_line.replace("\n", '').split(",")
            gene_id, gene_age, age_interval = arr[0], arr[1], arr[2]

            # we don't need gene age data, if gene is not loaded in genome
            if self.gen.feature_by_id(f"gene:{gene_id}") is None: continue

            data = self.INFO(gene_id, gene_age=gene_age, group_label=age_interval)
            assert not self.data_by_gene.__contains__(gene_id)
            self.data_by_gene[gene_id] = data

        print(f"Loaded AgeData (http://genorigin.chenzxlab.cn/) for {self.spec.short_name()}:\n\t"
              f"{len(self.data_by_gene.keys())} genes/records")

    class AgeGroupsType(enum.Enum):
        same_size_age_steps = 0,
        same_size_groups = 1,
        annotated_groups = 2

    gene_group_size = 800
    age_step_size = 50

    def calculate_groups(self, groups_type: AgeGroupsType):
        data = list(self.data_by_gene.values())
        groups = []

        if groups_type == self.AgeGroupsType.annotated_groups:
            data_by_annotated_label = {}
            for age_info in data:
                annotated_age_group_label = age_info.group_label
                if not data_by_annotated_label.__contains__(annotated_age_group_label):
                    data_by_annotated_label[annotated_age_group_label] = []
                data_by_annotated_label[annotated_age_group_label].append(age_info)
            arr = list(data_by_annotated_label.keys())
            arr.sort(key=lambda x: int(x.replace(">", "")) if x.__contains__(">") else int(
                x.split("-")[0]) if x.__contains__("-") else int(x))
            for interval in arr:
                groups.append(data_by_annotated_label[interval])

        if groups_type == self.AgeGroupsType.same_size_groups:
            data.sort(key=lambda x: x.get_gene_age_numeric())
            groups_count = (len(data) - 1) // self.gene_group_size + 1
            for index in range(groups_count):
                arr = []
                for i in range(index * self.gene_group_size, (index + 1) * self.gene_group_size):
                    if len(data) > i: arr.append(data[i])
                groups.append(arr)

        if groups_type == self.AgeGroupsType.same_size_age_steps:
            data.sort(key=lambda x: x.get_gene_age_numeric())
            for age_info in data:
                age = age_info.get_gene_age_numeric()
                index = (age - 1) // self.age_step_size
                while len(groups) <= index:
                    groups.append([])
                groups[index].append(age_info)
        return groups

    def get_group_label(self, group, group_index, age_group_type: AgeGroupsType):
        if age_group_type == self.AgeGroupsType.annotated_groups:
            assert len(group) > 0
            return group[0].group_label
        elif age_group_type == self.AgeGroupsType.same_size_groups:
            min_age, max_age = -1, -1
            for age_info in group:
                age = age_info.get_gene_age_numeric()
                min_age = age if min_age == -1 else min(age, min_age)
                max_age = age if max_age == -1 else max(age, max_age)
            return f"{min_age}-{max_age}"
        else:
            return f"{group_index * self.age_step_size + 1}-{(group_index + 1) * self.age_step_size}"

    def get_gene_ids_in_group(self, group):
        arr = []
        for age_info in group:
            gene_id = age_info.gene_id
            gene = self.gen.feature_by_id(f"gene:{gene_id}")
            assert gene is not None
            arr.append(f"gene:{gene_id}")
        return arr

    def get_gene_ids_sorted_by_age(self):
        arr = list(self.data_by_gene.items())
        arr.sort(key=lambda x: x[1].get_gene_age_numeric())
        data = []
        for g_id, info in arr:
            data.append((g_id, info.get_gene_age_numeric()))
        return data

    def get_gene_age(self, gene_id):
        gene_id = gene_id.replace("gene:", "")
        if self.data_by_gene.__contains__(gene_id):
            return self.data_by_gene[gene_id].get_gene_age_numeric()
        return None


def calculate_OGs(gen: GenomeWorker):
    graph = AnalyzerGraph()

    for chr_index in range(1, gen.chromosomes_count() + 1):
        genes_cnt = gen.genes_count_on_chr(chr_index)
        for i in range(0, genes_cnt):
            gene1 = gen.gene_by_indexes(chr_index, i)
            for j in range(i + 1, genes_cnt):
                gene2 = gen.gene_by_indexes(chr_index, j)
                ov_type = gen.get_overlap_type(gene1, gene2)
                if ov_type != OVERLAP_TYPE.NONE:
                    edge = AnalyzerGraph.GraphEdge(gene1.id, gene2.id)
                    graph.add_edge(edge)

    return graph


import statistics


# region Gene Parameters
def get_gene_GC(g_id, gen: GenomeWorker):
    gene = gen.feature_by_id(g_id)
    chr_id = gen.get_feature_chromosomal_position(g_id)[0]
    transcript = gen.get_transcript_from_gene_by_criteria(g_id, criteria=TRANSCRIPT_CRITERIA.LONGEST_CDS,
                                                          tie_breaker_criteria=TRANSCRIPT_CRITERIA.LONGEST_CDS_AND_UTRs)

    if transcript is None:
        print(gene.id)
        assert transcript is not None
    frags = gen.get_fragments_from_transcript(transcript.id)
    composition = [0, 0, 0, 0]
    for frag in frags:
        if frag.featuretype == "CDS":
            seq = gen.retrieve_feature_sequence(chr_id, frag).upper()
            stats = gen.sequence_composition(seq)
            composition[0] += stats[0]
            composition[1] += stats[1]
            composition[2] += stats[2]
            composition[3] += stats[3]
    return (composition[0] + composition[1]) / (composition[0] + composition[1] + composition[2] + composition[3])


def get_gene_CDS_length(g_id, gen: GenomeWorker):
    gene = gen.feature_by_id(g_id)
    chr_id = gen.get_feature_chromosomal_position(g_id)[0]
    transcript = gen.get_transcript_from_gene_by_criteria(g_id, criteria=TRANSCRIPT_CRITERIA.LONGEST_CDS,
                                                          tie_breaker_criteria=TRANSCRIPT_CRITERIA.LONGEST_CDS_AND_UTRs)

    if transcript is None:
        print(gene.id)
        assert transcript is not None

    frags = gen.get_fragments_from_transcript(transcript.id)
    length = 0
    for frag in frags:
        if frag.featuretype == "CDS":
            length += frag.end - frag.start + 1
    return length


def get_gene_length(g_id, gen: GenomeWorker):
    gene = gen.feature_by_id(g_id)
    return gene.end - gene.start + 1


# endregion

def get_genes_GC(g_ids, gen: GenomeWorker):
    values = []
    for g_id in g_ids:
        gene = gen.feature_by_id(g_id)
        if gene is None: continue
        val = get_gene_GC(g_id, gen)
        values.append(val)
    return statistics.median(values)


def get_genes_length(g_ids, gen: GenomeWorker):
    values = []
    for g_id in g_ids:
        gene = gen.feature_by_id(g_id)
        if gene is None: continue
        val = get_gene_length(g_id, gen)
        values.append(val)
    return statistics.median(values)


def get_genes_CDS_length(g_ids, gen: GenomeWorker):
    values = []
    for g_id in g_ids:
        gene = gen.feature_by_id(g_id)
        if gene is None: continue
        val = get_gene_CDS_length(g_id, gen)
        values.append(val)
    return statistics.median(values)


def get_genes_CUB(g_ids, gen: GenomeWorker):
    return 1


def get_genes_expression(g_ids, gen: GenomeWorker):
    return 1


def get_genes_methylation(g_ids, gen: GenomeWorker):
    return 1


def get_genes_pressure(g_ids, gen: GenomeWorker):
    # Ka/Ks ratio
    return 1


def get_genes_PPIN(g_ids, gen: GenomeWorker):
    return 1


def get_extra_ratio_data(g_ids, ogs_ids, gen):
    val1 = get_genes_GC(g_ids, gen)  # get_genes_GC(ogs_ids, gen) / get_genes_GC(g_ids, gen)
    val2 = get_genes_length(g_ids, gen)  # get_genes_length(ogs_ids, gen) / get_genes_length(g_ids, gen)
    val3 = get_genes_CDS_length(g_ids, gen)  # get_genes_length(ogs_ids, gen) / get_genes_length(g_ids, gen)
    val4 = get_genes_CUB(ogs_ids, gen) / get_genes_CUB(g_ids, gen)
    val5 = get_genes_expression(ogs_ids, gen) / get_genes_expression(g_ids, gen)
    val6 = get_genes_methylation(ogs_ids, gen) / get_genes_methylation(g_ids, gen)
    val7 = get_genes_pressure(ogs_ids, gen) / get_genes_pressure(g_ids, gen)
    val8 = get_genes_PPIN(ogs_ids, gen) / get_genes_PPIN(g_ids, gen)
    return f"{val1}\t{val2}\t{val3}\t{val4}\t{val5}\t{val6}\t{val7}\t{val8}"


#############################################################################################


specie_list = [SPECIES.Homo_sapiens,
               # SPECIES.Rattus_norvegicus,
               # SPECIES.Mus_musculus,
               # SPECIES.Sus_scrofa,
               # SPECIES.Danio_rerio,
               # SPECIES.Drosophila_melanogaster,
               # SPECIES.Arabidopsis_thaliana,
               ]

GO_graph = GOGraph()

for species in specie_list:
    # Load necessary data
    genome = GenomeWorker(species, ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_FRAGMENTS, SEQUENCE_LOAD.LOAD)
    age_data = AgeData(species, genome)

    graph = calculate_OGs(genome)
    gene_clusters = graph.get_connected_clusters()
    OGs = {}
    for cluster in gene_clusters:
        for gene_id in cluster:
            OGs[gene_id] = True

    outfile = open(f"./generated_data/agedata_{species.short_name()}.txt", "w")

    age_group_type = AgeData.AgeGroupsType.annotated_groups
    groups = age_data.calculate_groups(age_group_type)
    for group_index in range(len(groups)):
        group = groups[group_index]
        group_label = age_data.get_group_label(group, group_index, age_group_type)
        gene_ids = age_data.get_gene_ids_in_group(group)
        overlapped_ids = []
        for gene_id in gene_ids:
            if OGs.__contains__(gene_id):
                overlapped_ids.append(gene_id)
        extra_data = get_extra_ratio_data(gene_ids, overlapped_ids, genome)
        outfile.write(f"{group_index + 1}\t{group_label}\t{len(gene_ids)}\t{len(overlapped_ids)}\t{extra_data}\n")
    outfile.close()

    file = open("generated_data/all_genes_by_comps.txt", "w")
    data = age_data.get_gene_ids_sorted_by_age()
    arr2 = []
    arr_gc = []
    arr_len1 = []
    arr_len2 = []
    index = 0
    unknown = 0
    for g_id, age in data:
        g_id = "gene:" + g_id
        if genome.feature_by_id(g_id) == None:
            unknown += 1
            continue
        gc = get_gene_GC(g_id, genome)
        len1 = get_gene_length(g_id, genome)
        len2 = get_gene_CDS_length(g_id, genome)
        index += 1
        arr2.append(index)
        arr_gc.append(gc)
        arr_len1.append(len1)
        arr_len2.append(len2)
        file.write(f"{g_id.replace('gene:', '')}\t{age}\t{gc}\t{len1}\t{len2}\n")
    file.close()
    print(f"unknown: {unknown}")
    print(f"GC spearman: {scipy.stats.spearmanr(arr2, arr_gc)}")
    print(f"GC spearman: {scipy.stats.spearmanr(arr2, arr_len1)}")
    print(f"GC spearman: {scipy.stats.spearmanr(arr2, arr_len2)}")
#
# dna_genes = GO_graph.get_term_gene_list(genome, "GO:0006281")
# mito_genes = GO_graph.get_term_gene_list(genome, "GO:0005739")
# arr = []
# all_ogs = []
# for age_class, gene_ids in age2gene.items():
#     ignored = 0
#     used = 0
#     ogs = []
#     dna = 0
#     mito = 0
#     for gene_id in gene_ids:
#         gene = genome.feature_by_id("gene:" + gene_id)
#
#         if gene is None:
#             ignored += 1
#             continue
#         gene_sym = genome.get_gene_symbol(gene)
#         used += 1
#         if OGs.__contains__(gene_sym):
#             ogs.append(gene_sym)
#             all_ogs.append(gene_sym)
#         dna += 1 if dna_genes.__contains__(gene_sym) else 0
#         mito += 1 if mito_genes.__contains__(gene_sym) else 0
#     if age_class > 24:
#         otfile = open(f"{age_class}.txt", "w")
#         for gene_sym in ogs:
#             otfile.write(f"{gene_sym}\n")
#         otfile.close()
#     arr.append((age_class, ignored, used, len(ogs), dna, mito))
#
# otfile = open(f"all.txt", "w")
# for gene_sym in all_ogs:
#     otfile.write(f"{gene_sym}\n")
# otfile.close()
#
# arr.sort(key=lambda x: x[0])
# for age_class, ignored, used, ogs, dna, mito in arr:
#     out_file.write(f"{age_class}\t{ignored}\t{used}\t{ogs}\t{dna}\t{mito}\t\n")
