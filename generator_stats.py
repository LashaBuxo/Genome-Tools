import numpy

from worker_genome import *
from GO_graph import *
import scipy.stats

specie_list = [SPECIES.Homo_sapiens,
               SPECIES.Rattus_norvegicus,
               SPECIES.Mus_musculus,
               SPECIES.Sus_scrofa,
               SPECIES.Danio_rerio,
               SPECIES.Drosophila_melanogaster,
               SPECIES.Arabidopsis_thaliana
               ]


class Collected_Data:
    def __init__(self, species: SPECIES):
        self.species = species

        self.loaded_genes = 0
        self.filtered_genes = 0
        self.filtered_novel = 0
        self.filtered_repeat = 0
        self.filtered_nested = 0

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

        self.OGs_count_on_chr_by_type = []
        self.chr_label_and_genes_count = []

        self.used_genes = []
        self.OGs = []
        self.OGs_by_type = {}

    def filter_OGs_by_multi_type_overlap(self, data_by_type):
        types = [OVERLAP_TYPE.CONVERGENT, OVERLAP_TYPE.DIVERGENT, OVERLAP_TYPE.DIFF_NESTED, OVERLAP_TYPE.TANDEM,
                 OVERLAP_TYPE.SAME_NESTED, ]

        all_ids = {}
        for ov_type in types:
            g_ids = data_by_type[ov_type].keys()
            for g_id in g_ids:
                all_ids[g_id] = True
        all_g_ids = list(all_ids.keys())

        for g_id in all_g_ids:
            types_contains_id = []
            for ov_type in types:
                if data_by_type[ov_type].__contains__(g_id):
                    types_contains_id.append(ov_type)
            if len(types_contains_id) > 1:
                for type in types_contains_id:
                    data_by_type[type].pop(g_id)
        return data_by_type

    def get_gene_by_index(self, index):
        if len(self.used_genes) <= index: return ""
        return self.used_genes[index]

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

    def get_OG_gene_by_ov_and_index(self, ov_type: OVERLAP_TYPE, index):
        if not self.OGs_by_type.__contains__(ov_type):
            self.OGs_by_type[ov_type] = []
            x = self.OGs_count_by_type[ov_type].items()
            for gene_id, _ in x:
                self.OGs_by_type[ov_type].append(gene_id)
        if len(self.OGs_by_type[ov_type]) <= index: return ""
        return self.OGs_by_type[ov_type][index]

    def get_overlap_type_ratio(self, ov_type: OVERLAP_TYPE):
        val = "{:.1f}".format(len(self.OGs_count_by_type[ov_type]) / self.OGs_count * 100)
        return f'{val}%'

    def _filtered_novel(self, ignores):
        return (ignores["is_readthrough"] if ignores.__contains__("is_readthrough") else 0) + \
               (ignores["is_pseudogene"] if ignores.__contains__("is_pseudogene") else 0) \
               + (ignores["is_novel"] if ignores.__contains__("is_novel") else 0) \
               + (ignores["is_predicted"] if ignores.__contains__("is_predicted") else 0)

    def _filtered_repeat(self, ignores):
        return (ignores["name_duplicated"] if ignores.__contains__("name_duplicated") else 0) \
               + (ignores["acc_duplicated"] if ignores.__contains__("acc_duplicated") else 0)

    def filtered_by_types(self, gen: GenomeWorker):
        ignores = gen.ignored_genes_by_types
        self.filtered_novel = self._filtered_novel(ignores)
        self.filtered_repeat = self._filtered_repeat(ignores)
        self.filtered_nested = ignores["nested"] if ignores.__contains__("nested") else 0


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

    used_genes = {}

    OGs = {}
    genome = GenomeWorker(species, ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_FRAGMENTS, SEQUENCE_LOAD.NOT_LOAD)

    genes_on_chr = []
    OGs_on_chr_by_type = {OVERLAP_TYPE.MULTI: [],
                          OVERLAP_TYPE.CONVERGENT: [],
                          OVERLAP_TYPE.DIVERGENT: [],
                          OVERLAP_TYPE.DIFF_NESTED: [],
                          OVERLAP_TYPE.TANDEM: [],
                          OVERLAP_TYPE.SAME_NESTED: []}

    homologues_cnt = 0
    incomplete_scores = []
    isoforms_count = []

    for chr_index in range(1, genome.chromosomes_count() + 1):
        genes_cnt = genome.genes_count_on_chr(chr_index)
        genes_on_chr.append(genes_cnt)
        og = {OVERLAP_TYPE.MULTI: {},
              OVERLAP_TYPE.CONVERGENT: {},
              OVERLAP_TYPE.DIVERGENT: {},
              OVERLAP_TYPE.DIFF_NESTED: {},
              OVERLAP_TYPE.TANDEM: {},
              OVERLAP_TYPE.SAME_NESTED: {}}

        for i in range(0, genes_cnt):
            gene1 = genome.gene_by_indexes(chr_index, i)
            used_genes[gene1.id] = True
            homologues_cnt += 1 if reference_gene_symbols.__contains__(genome.get_gene_symbol(gene1)) else 0
            incomplete_scores.append(genome.get_gene_incomplete_level(gene1.id))
            isoforms_count.append(len(genome.get_transcripts_from_gene(gene1.id)))
            for j in range(i + 1, genes_cnt):
                gene2 = genome.gene_by_indexes(chr_index, j)
                ov_type = genome.get_overlap_type(gene1, gene2)
                if ov_type != OVERLAP_TYPE.NONE:
                    og[ov_type][gene1.id] = True
                    og[ov_type][gene2.id] = True
                    og[OVERLAP_TYPE.MULTI][gene1.id] = True
                    og[OVERLAP_TYPE.MULTI][gene2.id] = True

                    OGs[gene1.id] = True
                    OGs[gene2.id] = True

                    data.OGs_count_by_type[ov_type][gene1.id.replace("gene:", "")] = True
                    data.OGs_count_by_type[ov_type][gene2.id.replace("gene:", "")] = True

        og = data.filter_OGs_by_multi_type_overlap(og)  # filter multi types
        for ov_type in list(og.keys()):
            OGs_on_chr_by_type[ov_type].append(len(og[ov_type]))
        data.chr_label_and_genes_count.append((genes_cnt, genome.chr_index2_seq_id(chr_index)))

    data.OGs_count_on_chr_by_type = OGs_on_chr_by_type

    data.filtered_genes = genome.ignored_protein_coding_genes
    data.loaded_genes = genome.imported_protein_coding_genes

    data.filtered_by_types(genome)

    data.homologues_to_reference_species = homologues_cnt
    data.OGs_count = len(OGs)
    data.pearson_tests = scipy.stats.pearsonr(genes_on_chr, OGs_on_chr_by_type[OVERLAP_TYPE.MULTI])[0] * \
                         scipy.stats.pearsonr(genes_on_chr, OGs_on_chr_by_type[OVERLAP_TYPE.MULTI])[0]

    data.median_iso = numpy.median(isoforms_count)
    data.median_incomplete = numpy.median(incomplete_scores)

    for gene_id, _ in OGs.items():
        data.OGs.append(gene_id.replace("gene:", ""))

    for gene_id, _ in used_genes.items():
        data.used_genes.append(gene_id.replace("gene:", ""))

    data.OGs_count_by_type = data.filter_OGs_by_multi_type_overlap(data.OGs_count_by_type)

    species_data[species] = data

# print chromosome genes cnt, OGs
for species in specie_list:
    chr_file = open(f"./generated_data/chrdata_{species.short_name()}.txt", "w")
    chr_index = 0
    assert len(species_data[species].chr_label_and_genes_count) == len(
        species_data[species].OGs_count_on_chr_by_type[OVERLAP_TYPE.MULTI])
    for chr_genes_cnt, chr_label in species_data[species].chr_label_and_genes_count:
        ogs_data = ""
        for ov_type in list(species_data[species].OGs_count_on_chr_by_type.keys()):
            ogs_data += f"\t{species_data[species].OGs_count_on_chr_by_type[ov_type][chr_index]}"
        chr_file.write(f"{chr_label}\t{chr_genes_cnt}{ogs_data}\n")
        chr_index += 1
    chr_file.close()

# print OGs
OGs_file = open(f"generated_data/OGs_data.txt", 'w')
max_entries = 0
for species in specie_list:
    max_entries = max(max_entries, len(species_data[species].OGs))

for index in range(max_entries):
    for species in specie_list:
        OGs_file.write(f"{species_data[species].get_gene_by__index(index)}\t")
        OGs_file.write(f"{species_data[species].get_OG_gene_by_ov_and_index(OVERLAP_TYPE.CONVERGENT, index)}\t")
        OGs_file.write(f"{species_data[species].get_OG_gene_by_ov_and_index(OVERLAP_TYPE.DIVERGENT, index)}\t")
        OGs_file.write(f"{species_data[species].get_OG_gene_by_ov_and_index(OVERLAP_TYPE.DIFF_NESTED, index)}\t")
        OGs_file.write(f"{species_data[species].get_OG_gene_by_ov_and_index(OVERLAP_TYPE.TANDEM, index)}\t")
        OGs_file.write(f"{species_data[species].get_OG_gene_by_ov_and_index(OVERLAP_TYPE.SAME_NESTED, index)}\t")
    OGs_file.write(f"\n")
OGs_file.close()

# print All Genes
Genes_file = open(f"generated_data/Genes_data.txt", 'w')
max_entries = 0
for species in specie_list:
    max_entries = max(max_entries, len(species_data[species].used_genes))

for index in range(max_entries):
    for species in specie_list:
        Genes_file.write(f"{species_data[species].get_gene_by_index(index)}\t")
    Genes_file.write(f"\n")
Genes_file.close()

# print General Stats
file = open("generated_data/general_stats.txt", 'w')
file.write("General Data\n")
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

file.write("Filtered Genes in species (Incomplete,readthrough,predicted,repeat,nested,=)\n")
for species in specie_list:
    file.write(f'{species_data[species].loaded_genes + species_data[species].filtered_genes}\t'
               f'{species_data[species].filtered_novel}\t'
               f'{species_data[species].filtered_repeat}\t'
               f'{species_data[species].filtered_nested}\t'
               f'{species_data[species].loaded_genes}\n')

file.close()
