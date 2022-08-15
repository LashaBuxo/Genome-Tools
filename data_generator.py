import numpy

from worker_genome import *
from GO_graph import *
import scipy.stats

specie_list = [SPECIES.Homo_sapiens,
               SPECIES.Rattus_norvegicus,
               SPECIES.Mus_musculus,
               SPECIES.Sus_scrofa,
               SPECIES.Danio_rerio]

GO_graph = GOGraph()


class Collected_Data:
    def __init__(self, species: SPECIES):
        self.species = species
        self.filtered_genes = 0
        self.loaded_genes = 0
        self.median_iso = 0
        self.median_incomplete = 0
        self.homologues_to_reference_species = 0
        self.OGs_count = 0
        self.pearson_tests = (0, 0)
        self.OGs_count_by_type = {OVERLAP_TYPE.CONVERGENT: {},
                                  OVERLAP_TYPE.DIVERGENT: {},
                                  OVERLAP_TYPE.DIFF_NESTED: {},
                                  OVERLAP_TYPE.TANDEM: {},
                                  OVERLAP_TYPE.SAME_NESTED: {}
                                  }

        self.OGs = []
        self.OGs_MITO = []
        self.OGS_DNA = []
        self.OGs_MITO_1000 = []
        self.OGS_DNA_1000 = []

        self.ontology_sets_by_OGs_percentage = []
        self.OGs_by_type = {}

        self.OGs_by_length = {"=0": {}, "<100": {}, "<500": {}, "<1000": {}}

        self.edges = {}
        self.edges_1000 = {}
        self.OGs_by_cluster_size = {}

    def add_edge_1000(self, sym1, sym2):
        if not self.edges_1000.__contains__(sym1):
            self.edges_1000[sym1] = {}
        if not self.edges_1000.__contains__(sym2):
            self.edges_1000[sym2] = {}
        self.edges_1000[sym1][sym2] = True
        self.edges_1000[sym2][sym1] = True

    def add_edge(self, sym1, sym2):
        if not self.edges.__contains__(sym1):
            self.edges[sym1] = {}
        if not self.edges.__contains__(sym2):
            self.edges[sym2] = {}
        self.edges[sym1][sym2] = True
        self.edges[sym2][sym1] = True

    def calculate_clusters(self):
        nodes = self.edges.keys()
        assert len(nodes) == len(self.OGs)
        for node in nodes:
            viz = {}
            self.dfs(node, viz)
            sz = len(viz)
            self.OGs_by_cluster_size[sz] = self.OGs_by_cluster_size.get(sz, 0) + 1

    def get_OGs_cluster_ratio(self, size):
        answer = 0
        for sz, cnt in self.OGs_by_cluster_size.items():
            if sz > 5 and size > 5:
                answer += cnt
            elif sz == size:
                answer += cnt
        val = "{:.1f}".format(answer / self.loaded_genes * 100)
        return f'{val}%'

    def dfs(self, node, viz):
        viz[node] = True
        if not self.edges.__contains__(node): return
        neighbors = list(self.edges[node].keys())
        for neighbor in neighbors:
            if viz.__contains__(neighbor): continue
            self.dfs(neighbor, viz)

    def get_OGs_by_length_ratio(self, category):
        cnt = len(self.OGs_by_length[category])
        val = "{:.1f}".format(cnt / self.loaded_genes * 100)
        return f'{val}%'

    def add_OG_by_length(self, gene, dist):
        if dist == 0:
            self.OGs_by_length["=0"][gene] = True
            self.OGs_by_length["<100"][gene] = True
            self.OGs_by_length["<500"][gene] = True
            self.OGs_by_length["<1000"][gene] = True
        elif dist < 100:
            self.OGs_by_length["<100"][gene] = True
            self.OGs_by_length["<500"][gene] = True
            self.OGs_by_length["<1000"][gene] = True
        elif dist < 500:
            self.OGs_by_length["<500"][gene] = True
            self.OGs_by_length["<1000"][gene] = True
        elif dist < 1000:
            self.OGs_by_length["<1000"][gene] = True

    def get_INCOMPLETE_ratio(self):
        val = "{:.1f}".format(self.median_incomplete * 100)
        return f'{val}%'

    def get_Homolog_ratio(self):
        val = "{:.1f}".format(self.homologues_to_reference_species / self.loaded_genes * 100)
        return f'{val}%'

    def get_OGs_ratio(self):
        val = "{:.1f}".format(self.OGs_count / self.loaded_genes * 100)
        return f'{val}%'

    def get_gene_by__index(self, index):
        if len(self.OGs) <= index: return ""
        return self.OGs[index]

    def get_gene_by_ov_and_index(self, ov_type: OVERLAP_TYPE, index):
        if not self.OGs_by_type.__contains__(ov_type):
            self.OGs_by_type[ov_type] = []
            x = self.OGs_count_by_type[ov_type].items()
            for gene_sym, _ in x:
                self.OGs_by_type[ov_type].append(gene_sym)
        if len(self.OGs_by_type[ov_type]) <= index: return ""
        return self.OGs_by_type[ov_type][index]

    def get_overlap_type_ratio(self, ov_type: OVERLAP_TYPE):
        total = 0
        for _, dic in self.OGs_count_by_type.items():
            total += len(dic)

        val = "{:.1f}".format(len(self.OGs_count_by_type[ov_type]) / total * 100)
        return f'{val}%'


reference_gene_symbols = {}
reference_species = SPECIES.Homo_sapiens
genome = GenomeWorker(reference_species, ANNOTATION_LOAD.GENES, SEQUENCE_LOAD.NOT_LOAD)
for chr_index in range(1, genome.chromosomes_count() + 1):
    genes_cnt = genome.genes_count_on_chr(chr_index)
    for i in range(0, genes_cnt):
        gene1 = genome.gene_by_indexes(chr_index, i)
        reference_gene_symbols[genome.get_gene_symbol(gene1)] = True

species_data = {}
max_sets_in_species = 0
for species in specie_list:
    data = Collected_Data(species)

    OGs = {}
    genome = GenomeWorker(species, ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_FRAGMENTS, SEQUENCE_LOAD.NOT_LOAD)

    genes_on_chr = []
    OGs_on_chr = []
    homologues_cnt = 0
    incomplete_scores = []
    isoforms_count = []

    for chr_index in range(1, genome.chromosomes_count() + 1):
        genes_cnt = genome.genes_count_on_chr(chr_index)
        genes_on_chr.append(genes_cnt)
        og = {}
        for i in range(0, genes_cnt):
            gene1 = genome.gene_by_indexes(chr_index, i)
            sym1 = genome.get_gene_symbol(gene1)
            homologues_cnt += 1 if reference_gene_symbols.__contains__(genome.get_gene_symbol(gene1)) else 0
            incomplete_scores.append(genome.get_gene_incomplete_level(gene1.id))
            isoforms_count.append(len(genome.get_transcripts_from_gene(gene1.id)))
            for j in range(i + 1, genes_cnt):
                gene2 = genome.gene_by_indexes(chr_index, j)
                sym2 = genome.get_gene_symbol(gene2)

                dist = genome.get_features_distance(gene1, gene2)
                data.add_OG_by_length(sym1, dist)
                data.add_OG_by_length(sym2, dist)

                if dist < 1000:
                    data.add_edge_1000(sym1, sym2)

                ov_type = genome.get_overlap_type(gene1, gene2)
                if ov_type != OVERLAP_TYPE.NONE:
                    OGs[sym1] = True
                    OGs[sym2] = True
                    og[sym1] = True
                    og[sym2] = True
                    data.add_edge(sym1, sym2)
                    data.OGs_count_by_type[ov_type][sym1] = True
                    data.OGs_count_by_type[ov_type][sym2] = True

        OGs_on_chr.append(len(og))

    data.filtered_genes = genome.ignored_protein_coding_genes
    data.loaded_genes = genome.imported_protein_coding_genes
    data.homologues_to_reference_species = homologues_cnt
    data.OGs_count = len(OGs)
    data.pearson_tests = scipy.stats.pearsonr(genes_on_chr, OGs_on_chr)[0] * \
                         scipy.stats.pearsonr(genes_on_chr, OGs_on_chr)[0]

    data.median_iso = numpy.median(isoforms_count)
    data.median_incomplete = numpy.median(incomplete_scores)

    for gene_sym, _ in OGs.items():
        data.OGs.append(gene_sym)

    data.ontology_sets_by_OGs_percentage = GO_graph.calculate_OGs_in_GOs(genome, OGs, data.OGs_count_by_type, 100)
    data.ontology_sets_by_OGs_percentage.sort(key=lambda x: x[5], reverse=True)
    max_sets_in_species = max(max_sets_in_species, len(data.ontology_sets_by_OGs_percentage))

    data.OGS_DNA = GO_graph.get_term_genes(genome, "GO:0006281", data.OGs)
    data.OGS_DNA_1000 = GO_graph.get_term_genes(genome, "GO:0006281", data.OGs_by_length["<1000"])
    data.OGs_MITO = GO_graph.get_term_genes(genome, "GO:0005739", data.OGs)
    data.OGs_MITO_1000 = GO_graph.get_term_genes(genome, "GO:0005739", data.OGs_by_length["<1000"])

    data.calculate_clusters()
    species_data[species] = data

OGs_file = open(f"generated_data/OGs_data.txt", 'w')
for index in range(0, len(species_data[SPECIES.Homo_sapiens].OGs)):
    for species in specie_list:
        OGs_file.write(f"{species_data[species].get_gene_by__index(index)}\t")
        OGs_file.write(f"{species_data[species].get_gene_by_ov_and_index(OVERLAP_TYPE.CONVERGENT, index)}\t")
        OGs_file.write(f"{species_data[species].get_gene_by_ov_and_index(OVERLAP_TYPE.DIVERGENT, index)}\t")
        OGs_file.write(f"{species_data[species].get_gene_by_ov_and_index(OVERLAP_TYPE.DIFF_NESTED, index)}\t")
        OGs_file.write(f"{species_data[species].get_gene_by_ov_and_index(OVERLAP_TYPE.TANDEM, index)}\t")
        OGs_file.write(f"{species_data[species].get_gene_by_ov_and_index(OVERLAP_TYPE.SAME_NESTED, index)}\t")
    OGs_file.write(f"\n")
OGs_file.close()

file = open("generated_data/old2/general_stats.txt", 'w')
for species in specie_list:
    file.write(f'{species.short_name()}\t{species_data[species].filtered_genes}\t{species_data[species].loaded_genes}\t'
               f'{species_data[species].median_iso}\t'
               f'{species_data[species].get_INCOMPLETE_ratio()}\t'
               f'{species_data[species].homologues_to_reference_species}\t'
               f'{species_data[species].get_Homolog_ratio()}\t'
               f'{species_data[species].OGs_count}\t'
               f'{species_data[species].get_OGs_ratio()}\t'
               f'{species_data[species].pearson_tests}\t'
               f'{species_data[species].get_overlap_type_ratio(OVERLAP_TYPE.CONVERGENT)}\t'
               f'{species_data[species].get_overlap_type_ratio(OVERLAP_TYPE.DIVERGENT)}\t'
               f'{species_data[species].get_overlap_type_ratio(OVERLAP_TYPE.DIFF_NESTED)}\t'
               f'{species_data[species].get_overlap_type_ratio(OVERLAP_TYPE.TANDEM)}\t'
               f'{species_data[species].get_overlap_type_ratio(OVERLAP_TYPE.SAME_NESTED)}\t\n')
file.write("\n")

for species in specie_list:
    file.write(f'{species.short_name()}\t'
               f'{species_data[species].get_OGs_by_length_ratio("=0")}\t'
               f'{species_data[species].get_OGs_by_length_ratio("<100")}\t'
               f'{species_data[species].get_OGs_by_length_ratio("<500")}\t'
               f'{species_data[species].get_OGs_by_length_ratio("<1000")}\t\n'
               )
file.write("\n")

for species in specie_list:
    file.write(f'{species.short_name()}\t'
               f'{species_data[species].get_OGs_cluster_ratio(2)}\t'
               f'{species_data[species].get_OGs_cluster_ratio(3)}\t'
               f'{species_data[species].get_OGs_cluster_ratio(4)}\t'
               f'{species_data[species].get_OGs_cluster_ratio(5)}\t'
               f'{species_data[species].get_OGs_cluster_ratio(6)}\t\n'
               )
file.write("\n")

for i in range(0, max_sets_in_species):
    file.write('\n')
    for species in specie_list:
        if len(species_data[species].ontology_sets_by_OGs_percentage) <= i:
            file.write(f'\t\t\t\t\t\t\t\t\t\t\t')
        else:
            file.write(f'{species_data[species].ontology_sets_by_OGs_percentage[i][0]}\t'
                       f'{species_data[species].ontology_sets_by_OGs_percentage[i][1]}\t'
                       f'{species_data[species].ontology_sets_by_OGs_percentage[i][2]}\t'
                       f'{species_data[species].ontology_sets_by_OGs_percentage[i][3]}\t'
                       f'{species_data[species].ontology_sets_by_OGs_percentage[i][4]}\t'
                       f'{species_data[species].ontology_sets_by_OGs_percentage[i][5]}%\t'
                       f'{species_data[species].ontology_sets_by_OGs_percentage[i][6]}\t'
                       f'{species_data[species].ontology_sets_by_OGs_percentage[i][7]}\t'
                       f'{species_data[species].ontology_sets_by_OGs_percentage[i][8]}\t'
                       f'{species_data[species].ontology_sets_by_OGs_percentage[i][9]}\t'
                       f'{species_data[species].ontology_sets_by_OGs_percentage[i][10]}\t')

file.close()


def calculate_cons(genes1, genes2, data1: Collected_Data, data2: Collected_Data, ov_type1, ov_type2, paired_same):
    cnt = 0
    if paired_same == False:
        for gene1 in genes1:
            if genes2.__contains__(gene1):
                cnt += 1
    else:
        for gene in genes1:
            if genes2.__contains__(gene):
                neighbors1 = list(data1.edges[gene].keys()) if ov_type1 == "=0" else list(
                    data1.edges_1000[gene].keys())
                neighbors2 = list(data2.edges[gene].keys()) if ov_type2 == "=0" else list(
                    data2.edges_1000[gene].keys())
                flag = False
                for neighbor in neighbors1:
                    if neighbors2.__contains__(neighbor):
                        flag = True
                        break
                if flag:
                    cnt += 1
    return cnt


def OGs_cons(data1: Collected_Data, data2: Collected_Data, ov_type1, ov_type2, paired_same):
    genes1 = list(data1.OGs_by_length[ov_type1].keys())
    genes2 = list(data2.OGs_by_length[ov_type2].keys())
    return calculate_cons(genes1, genes2, data1, data2, ov_type1, ov_type2, paired_same)


def OGs_DNA_cons(data1: Collected_Data, data2: Collected_Data, ov_type1, ov_type2, paired_same):
    genes1 = data1.OGS_DNA if ov_type1 == "=0" else data1.OGS_DNA_1000
    genes2 = data2.OGS_DNA if ov_type2 == "=0" else data2.OGS_DNA_1000
    return calculate_cons(genes1, genes2, data1, data2, ov_type1, ov_type2, paired_same)


def OGs_MITO_cons(data1: Collected_Data, data2: Collected_Data, ov_type1, ov_type2, paired_same):
    genes1 = data1.OGs_MITO if ov_type1 == "=0" else data1.OGs_MITO_1000
    genes2 = data2.OGs_MITO if ov_type2 == "=0" else data2.OGs_MITO_1000
    return calculate_cons(genes1, genes2, data1, data2, ov_type1, ov_type2, paired_same)


file = open("generated_data/old/OGs_cons.txt", 'w')
for species_A in specie_list:
    for species_B in specie_list:
        file.write(f'{OGs_cons(species_data[species_A], species_data[species_B], "=0", "=0", paired_same=False)}\t')
    file.write("\n")
file.write("\n")

for species_A in specie_list:
    for species_B in specie_list:
        file.write(f'{OGs_cons(species_data[species_A], species_data[species_B], "=0", "<1000", paired_same=False)}\t')
    file.write("\n")
file.write("\n")

for species_A in specie_list:
    for species_B in specie_list:
        file.write(f'{OGs_cons(species_data[species_A], species_data[species_B], "=0", "=0", paired_same=True)}\t')
    file.write("\n")
file.write("\n")

for species_A in specie_list:
    for species_B in specie_list:
        file.write(f'{OGs_cons(species_data[species_A], species_data[species_B], "=0", "<1000", paired_same=True)}\t')
    file.write("\n")
file.write("\n")
file.close()

file = open("generated_data/old/OGs_DNA_cons.txt", 'w')
for species_A in specie_list:
    for species_B in specie_list:
        file.write(f'{OGs_DNA_cons(species_data[species_A], species_data[species_B], "=0", "=0", paired_same=False)}\t')
    file.write("\n")
file.write("\n")

for species_A in specie_list:
    for species_B in specie_list:
        file.write(
            f'{OGs_DNA_cons(species_data[species_A], species_data[species_B], "=0", "<1000", paired_same=False)}\t')
    file.write("\n")
file.write("\n")

for species_A in specie_list:
    for species_B in specie_list:
        file.write(f'{OGs_DNA_cons(species_data[species_A], species_data[species_B], "=0", "=0", paired_same=True)}\t')
    file.write("\n")
file.write("\n")

for species_A in specie_list:
    for species_B in specie_list:
        file.write(
            f'{OGs_DNA_cons(species_data[species_A], species_data[species_B], "=0", "<1000", paired_same=True)}\t')
    file.write("\n")
file.write("\n")
file.close()

file = open("generated_data/old/OGs_MITO_cons.txt", 'w')
for species_A in specie_list:
    for species_B in specie_list:
        file.write(
            f'{OGs_MITO_cons(species_data[species_A], species_data[species_B], "=0", "=0", paired_same=False)}\t')
    file.write("\n")
file.write("\n")

for species_A in specie_list:
    for species_B in specie_list:
        file.write(
            f'{OGs_MITO_cons(species_data[species_A], species_data[species_B], "=0", "<1000", paired_same=False)}\t')
    file.write("\n")
file.write("\n")

for species_A in specie_list:
    for species_B in specie_list:
        file.write(f'{OGs_MITO_cons(species_data[species_A], species_data[species_B], "=0", "=0", paired_same=True)}\t')
    file.write("\n")
file.write("\n")

for species_A in specie_list:
    for species_B in specie_list:
        file.write(
            f'{OGs_MITO_cons(species_data[species_A], species_data[species_B], "=0", "<1000", paired_same=True)}\t')
    file.write("\n")


def get_paired_GOs(out_file, ov_genes):
    fixed = {}
    for og in ov_genes:
        if fixed.__contains__(og): continue
        out_file.write(f"{og}\t")
        neighbors = list(species_data[species].edges[og].keys())
        for neighbor in neighbors:
            fixed[neighbor] = True
            out_file.write(f"{neighbor.upper() if ov_genes.__contains__(neighbor) else neighbor.lower()}\t")
            term = GO_graph.get_gene_term(neighbor)
            out_file.write(f"{term}\t")
        out_file.write("\n")


for species in specie_list:
    file = open(f"generated_data/Enriched GOs of Pairs[{species.short_name()}].txt", 'w')
    get_paired_GOs(file, species_data[species].OGS_DNA)
    file.write("\n")
    file.write("\n")
    get_paired_GOs(file, species_data[species].OGs_MITO)
    file.close()
