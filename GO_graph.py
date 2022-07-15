from worker_genome import *
import matplotlib.pyplot as plt
import networkx as nx
import pydot
from networkx.drawing.nx_pydot import graphviz_layout
import numpy


class GeneAnnotations:
    def __init__(self, gene_symbol, terms, gene2term_relationships, evidence_codes):
        self.symbol = gene_symbol
        self.related_terms = terms
        self.relationships = gene2term_relationships
        self.evidence_codes = evidence_codes


class GONode:
    def __init__(self, go_id, go_name, go_namespace):
        self.id = go_id
        self.name = go_name
        self.namespace = go_namespace

        self.related_gene_annotations = []
        self.all_related_lower_gene_symbols = {}

    def add_related_gene_record(self, gene_record: GeneAnnotations):
        self.related_gene_annotations.append(gene_record)
        self.all_related_lower_gene_symbols[gene_record.symbol] = True


class GOGraph:
    graph_file = "used_data/GO_data/GO_graph/go-basic.obo"
    annotations_file = "used_data/GO_data/human_GO_annotations/GOAnnotations_all.txt"

    def __init__(self):
        # self.genome = genome

        self.GO_nodes = {}
        self.inner_edges = {}
        self.outer_edges = {}
        self.root_nodes = []
        self._load_OG_graph()

        self.GENE_annotations = {}
        self._load_OG_annotations()

        for root_node in self.root_nodes:
            print(f'{root_node.name} - related lower genes {len(root_node.all_related_lower_gene_symbols)} ')

        # slim_generic_ids = "GO:1901135,GO:0140053,GO:0140014,GO:0140013,GO:0098754,GO:0098542,GO:0072659,GO:0071941,GO:0071554,GO:0065003,GO:0061024,GO:0061007,GO:0055086,GO:0055085,GO:0055065,GO:0051604,GO:0050886,GO:0050877,GO:0048870,GO:0048856,GO:0044782,GO:0042254,GO:0042060,GO:0036211,GO:0034330,GO:0032200,GO:0031047,GO:0030198,GO:0030163,GO:0030154,GO:0023052,GO:0022600,GO:0022414,GO:0016192,GO:0016073,GO:0016071,GO:0015979,GO:0012501,GO:0007568,GO:0007163,GO:0007155,GO:0007059,GO:0007040,GO:0007031,GO:0007018,GO:0007010,GO:0007005,GO:0006954,GO:0006914,GO:0006913,GO:0006886,GO:0006790,GO:0006766,GO:0006629,GO:0006575,GO:0006520,GO:0006486,GO:0006457,GO:0006399,GO:0006355,GO:0006351,GO:0006325,GO:0006310,GO:0006260,GO:0006091,GO:0005975,GO:0003016,GO:0003014,GO:0003013,GO:0003012,GO:0002376,GO:0002181,GO:0000910,GO:0000278,GO:0140691,GO:0140657,GO:0140313,GO:0140299,GO:0140223,GO:0140110,GO:0140104,GO:0140098,GO:0140097,GO:0140096,GO:0120274,GO:0098772,GO:0098631,GO:0090729,GO:0060090,GO:0060089,GO:0048018,GO:0045735,GO:0045182,GO:0044183,GO:0042393,GO:0038024,GO:0031386,GO:0016874,GO:0016853,GO:0016829,GO:0016787,GO:0016740,GO:0016491,GO:0016209,GO:0009975,GO:0008289,GO:0008092,GO:0005215,GO:0005198,GO:0003924,GO:0003824,GO:0003774,GO:0003723,GO:0003677,GO:0001618,GO:0043226,GO:0031410,GO:0031012,GO:0030312,GO:0009579,GO:0009536,GO:0005929,GO:0005886,GO:0005856,GO:0005840,GO:0005829,GO:0005815,GO:0005811,GO:0005794,GO:0005783,GO:0005777,GO:0005773,GO:0005768,GO:0005764,GO:0005730,GO:0005694,GO:0005654,GO:0005635,GO:0005634,GO:0005618,GO:0005615,GO:0005576,GO:0000228,GO:0005739,GO:0006281"
        # terms = slim_generic_ids.split(',')
        # s = ""
        # for term in terms:
        #     go_node = self.GO_nodes[term]
        #     if go_node.namespace == "cellular_component":
        #         s = s + f",{term}"
        # print(s)
        # G = nx.DiGraph()
        #
        # for term, node in self.GO_nodes.items():
        #     if len(node.all_related_lower_gene_symbols) > 1:
        #         G.add_node(term, color="green")
        #
        # for upper_term, edges in self.outer_edges.items():
        #     for lower_term, relationship in edges:
        #         if len(self.GO_nodes[lower_term].all_related_lower_gene_symbols) > 1:
        #             G.add_edge(upper_term, lower_term)
        #         # G.add_edge(self.GO_nodes[0][0], self.GO_nodes[1][0])
        #
        # nx.draw(G, with_labels=False, arrows=False)
        # plt.show()

    def get_term_gene_list(self, genome: GenomeWorker, term):
        arr = []
        go_node = self.GO_nodes[term]
        for gene_symbol, _ in go_node.all_related_lower_gene_symbols.items():
            gene = genome.gene_by_symbol(gene_symbol)
            if gene is not None:
                arr.append(gene_symbol)
        return arr

    def get_term_genes(self, genome: GenomeWorker, term, OGs):
        arr = []
        go_node = self.GO_nodes[term]
        for gene_symbol, _ in go_node.all_related_lower_gene_symbols.items():
            gene = genome.gene_by_symbol(gene_symbol)
            if gene is not None and OGs.__contains__(gene_symbol):
                arr.append(gene_symbol)
        return arr

    def get_gene_term(self, gene_sym):
        if not self.GENE_annotations.__contains__(gene_sym): return "Unknown"
        possible_terms = self.GENE_annotations[gene_sym].related_terms
        best = ""
        mini = 100000
        for term in possible_terms:
            go_node = self.GO_nodes[term]
            l = len(go_node.all_related_lower_gene_symbols.items())
            if l < mini and l >= 100:
                best = f"{go_node.name} ({term})"
                mini = l
        return best

    def calculate_OGs_in_GOs(self, genome: GenomeWorker, OGs, OGs_count_by_type, minimum_size):
        results = []

        gene_types = {}

        for type, dic in OGs_count_by_type.items():
            for gene_sym, _ in dic.items():
                if not gene_types.__contains__(gene_sym):
                    gene_types[gene_sym] = []
                gene_types[gene_sym].append(type)

        slim_generic_ids = "GO:1901135,GO:0140053,GO:0140014,GO:0140013,GO:0098754,GO:0098542,GO:0072659,GO:0071941,GO:0071554,GO:0065003,GO:0061024,GO:0061007,GO:0055086,GO:0055085,GO:0055065,GO:0051604,GO:0050886,GO:0050877,GO:0048870,GO:0048856,GO:0044782,GO:0042254,GO:0042060,GO:0036211,GO:0034330,GO:0032200,GO:0031047,GO:0030198,GO:0030163,GO:0030154,GO:0023052,GO:0022600,GO:0022414,GO:0016192,GO:0016073,GO:0016071,GO:0015979,GO:0012501,GO:0007568,GO:0007163,GO:0007155,GO:0007059,GO:0007040,GO:0007031,GO:0007018,GO:0007010,GO:0007005,GO:0006954,GO:0006914,GO:0006913,GO:0006886,GO:0006790,GO:0006766,GO:0006629,GO:0006575,GO:0006520,GO:0006486,GO:0006457,GO:0006399,GO:0006355,GO:0006351,GO:0006325,GO:0006310,GO:0006260,GO:0006091,GO:0005975,GO:0003016,GO:0003014,GO:0003013,GO:0003012,GO:0002376,GO:0002181,GO:0000910,GO:0000278,GO:0140691,GO:0140657,GO:0140313,GO:0140299,GO:0140223,GO:0140110,GO:0140104,GO:0140098,GO:0140097,GO:0140096,GO:0120274,GO:0098772,GO:0098631,GO:0090729,GO:0060090,GO:0060089,GO:0048018,GO:0045735,GO:0045182,GO:0044183,GO:0042393,GO:0038024,GO:0031386,GO:0016874,GO:0016853,GO:0016829,GO:0016787,GO:0016740,GO:0016491,GO:0016209,GO:0009975,GO:0008289,GO:0008092,GO:0005215,GO:0005198,GO:0003924,GO:0003824,GO:0003774,GO:0003723,GO:0003677,GO:0001618,GO:0043226,GO:0031410,GO:0031012,GO:0030312,GO:0009579,GO:0009536,GO:0005929,GO:0005886,GO:0005856,GO:0005840,GO:0005829,GO:0005815,GO:0005811,GO:0005794,GO:0005783,GO:0005777,GO:0005773,GO:0005768,GO:0005764,GO:0005730,GO:0005694,GO:0005654,GO:0005635,GO:0005634,GO:0005618,GO:0005615,GO:0005576,GO:0000228,GO:0005739,GO:0006281"
        terms = slim_generic_ids.split(',')
        for term in terms:
            OGs_local = {OVERLAP_TYPE.CONVERGENT: 0,
                         OVERLAP_TYPE.DIVERGENT: 0,
                         OVERLAP_TYPE.DIFF_NESTED: 0,
                         OVERLAP_TYPE.SAME_NESTED: 0,
                         OVERLAP_TYPE.TANDEM: 0
                         }

            go_node = self.GO_nodes[term]
            ignored_genes = 0
            used_genes = 0
            inc_scores = []
            iso_count = []
            ov_genes = 0
            for gene_symbol, _ in go_node.all_related_lower_gene_symbols.items():
                gene = genome.gene_by_symbol(gene_symbol)
                if gene is None:
                    ignored_genes += 1
                else:
                    used_genes += 1
                    ov_genes += 1 if OGs.__contains__(gene_symbol) else 0
                    if gene_types.__contains__(gene_symbol):
                        for typ in gene_types[gene_symbol]:
                            OGs_local[typ] += 1
                    inc_scores.append(genome.get_gene_incomplete_level(gene.id))
                    iso_count.append(len(genome.get_transcripts_from_gene(gene.id)))

            if used_genes >= minimum_size:
                results.append((go_node.name,
                                ignored_genes,
                                used_genes,
                                numpy.median(iso_count),
                                self.__ratio_to_percent(numpy.median(inc_scores)),
                                float("{:.1f}".format(ov_genes / used_genes * 100)),
                                self.__ratio_ov_type(OVERLAP_TYPE.CONVERGENT, OGs_local),
                                self.__ratio_ov_type(OVERLAP_TYPE.DIVERGENT, OGs_local),
                                self.__ratio_ov_type(OVERLAP_TYPE.DIFF_NESTED, OGs_local),
                                self.__ratio_ov_type(OVERLAP_TYPE.TANDEM, OGs_local),
                                self.__ratio_ov_type(OVERLAP_TYPE.SAME_NESTED, OGs_local),
                                ))
        return results

    def __ratio_ov_type(self, ov_type: OVERLAP_TYPE, OGs_count_by_type):
        total = 0
        for _, cnt in OGs_count_by_type.items():
            total += cnt

        val = "{:.1f}".format(OGs_count_by_type[ov_type] / total * 100)
        return f'{val}%'

    def __ratio_to_percent(self, val):
        val = "{:.1f}".format(val * 100)
        return f'{val}%'

    def _load_OG_graph(self):
        file = open(self.graph_file, 'r')
        lines = file.readlines()
        lower_term, lower_name, lower_namespace = "", "", ""
        for index in range(0, len(lines)):
            lines[index] = lines[index].replace('\n', '')
            if index < 3: continue
            if lines[index].startswith("[Typedef]"):
                break
            if lines[index - 3].startswith("[Term]"):
                lower_term = lines[index - 2].split('id: ')[1]
                lower_name = lines[index - 1].split('name: ')[1]
                lower_namespace = lines[index].split('namespace: ')[1]
                self.__add_node(lower_term, lower_name, lower_namespace)

            if lines[index].startswith("is_a:"):
                relationship = "is_a"
                upper_term = lines[index].split(' ')[1]
                self.__add_edges(lower_term, upper_term, relationship)
            if lines[index].startswith("relationship:"):
                relationship = lines[index].split(' ')[1]
                upper_term = lines[index].split(' ')[2]
                self.__add_edges(lower_term, upper_term, relationship)
        for term, node in self.GO_nodes.items():
            if len(self.outer_edges[term]) == 0:
                self.root_nodes.append(node)

    def __add_node(self, term, name, namespace):
        if name.startswith("obsolete"): return
        node = GONode(term, name, namespace)
        assert not self.GO_nodes.__contains__(term)
        self.GO_nodes[term] = node

    def __create_edge_records_if_not_exist(self, term):
        if not self.outer_edges.__contains__(term):
            self.outer_edges[term] = []
        if not self.inner_edges.__contains__(term):
            self.inner_edges[term] = []

    def __add_edges(self, lower_term, upper_term, relationship):
        self.__create_edge_records_if_not_exist(lower_term)
        self.__create_edge_records_if_not_exist(upper_term)

        self.outer_edges[lower_term].append((upper_term, relationship))
        self.inner_edges[upper_term].append((lower_term, relationship))

    def __array_string_to_items(self, s: str):
        items = s.split(',')
        for i in range(0, len(items)):
            items[i] = items[i].replace('[', '').replace(']', '').replace(' ', '').replace("'", "")
        return items

    def _load_OG_annotations(self):
        file = open(self.annotations_file, 'r')
        lines = file.readlines()
        for line in lines:
            items = line.replace('\n', '').split('\t')

            gene_symbol = items[0]
            go_terms = self.__array_string_to_items(items[3])
            gene2term_relationships = self.__array_string_to_items(items[5])
            evidence_codes = self.__array_string_to_items(items[6])

            gene_annotations = GeneAnnotations(gene_symbol, go_terms, gene2term_relationships, evidence_codes)
            self.GENE_annotations[gene_symbol] = gene_annotations

            for term in go_terms:
                assert self.GO_nodes.__contains__(term)
                self.GO_nodes[term].add_related_gene_record(gene_annotations)

        # calculate_nodes_all_related_inner_genes
        for root_node in self.root_nodes:
            self.__DFS(root_node)

    def __DFS(self, node: GONode):
        # print(f'{node.name} {len(node.all_related_lower_gene_symbols)}')
        for lower_term, relationship in self.inner_edges[node.id]:
            lower_node = self.GO_nodes[lower_term]
            self.__DFS(lower_node)
            for gene_symbol, _ in lower_node.all_related_lower_gene_symbols.items():
                node.all_related_lower_gene_symbols[gene_symbol] = True
