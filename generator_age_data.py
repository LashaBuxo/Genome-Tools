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


#############################################################################################


specie_list = [SPECIES.Homo_sapiens,
               SPECIES.Rattus_norvegicus,
               SPECIES.Mus_musculus,
               SPECIES.Sus_scrofa,
               SPECIES.Danio_rerio,
               SPECIES.Drosophila_melanogaster,
               SPECIES.Arabidopsis_thaliana,
               ]

GO_graph = GOGraph()

######################
OGs_by_type = {}
ov_types = [OVERLAP_TYPE.MULTI, OVERLAP_TYPE.CONVERGENT, OVERLAP_TYPE.DIVERGENT, OVERLAP_TYPE.DIFF_NESTED,
            OVERLAP_TYPE.TANDEM,
            OVERLAP_TYPE.SAME_NESTED]
for specie in specie_list:
    OGs_by_type[specie] = {}
    for ov_type in ov_types:
        OGs_by_type[specie][ov_type] = []

file = open("generated_data/OGs_data.txt", "r")
lines = file.readlines()
for line in lines:
    arr = line.split('\t')
    for index in range(42):
        specie = specie_list[index // 6]
        ov_type = ov_types[index % 6]
        OGs_by_type[specie][ov_type].append(arr[index])

######################
for species in specie_list:
    # Load necessary data
    genome = GenomeWorker(species, ANNOTATION_LOAD.GENES, SEQUENCE_LOAD.NOT_LOAD)
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

        # overlapped_ids = []
        # for gene_id in gene_ids:
        #     if OGs.__contains__(gene_id):
        #         overlapped_ids.append(gene_id)
        extra_data = ""
        for ov_type in ov_types:
            if ov_type==OVERLAP_TYPE.SAME_NESTED: continue
            overlapped_ids = OGs_by_type[species][ov_type]
            ov_count = 0
            for gene_id in gene_ids:
                gene_id = gene_id.replace("gene:", "")
                if overlapped_ids.__contains__(gene_id):
                    ov_count += 1
            extra_data += f"\t{ov_count}"
        outfile.write(f"{group_index + 1}\t{group_label}\t{len(gene_ids)}{extra_data}\n")
    outfile.close()

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
