from Bio.Seq import Seq
from gffutils import Feature

from genome_worker_enums import *


class AnalyzerGraph:
    class GraphEdge:
        def __init__(self, gene1_id, gene2_id, edge_type: OVERLAP_TYPE):
            self.node1 = gene1_id
            self.node2 = gene2_id
            self.edge_type = edge_type

    def __init__(self):
        self.neighbors = {}

    def add_edge(self, edge: GraphEdge):
        self.__add_neighbor(edge.node1, (edge.node2, edge.edge_type))
        self.__add_neighbor(edge.node2, (edge.node1, edge.edge_type))

    def __add_neighbor(self, node, neighbor_data):
        if not self.neighbors.__contains__(node):
            self.neighbors[node] = []
        self.neighbors[node].append(neighbor_data)

    _visited = {}

    def dfs(self, node, allowed_edge_type):
        self._visited[node] = True
        if not self.neighbors.__contains__(node): return []
        cluster_members = [node]
        for (neighbor, edge_type) in self.neighbors[node]:
            if edge_type is allowed_edge_type and not self._visited[neighbor]:
                cluster_members += self.dfs(neighbor, allowed_edge_type)
        return cluster_members

    def get_connected_clusters(self, overlap_type: OVERLAP_TYPE):
        nodes = self.neighbors.keys()

        for node in nodes:
            self._visited[node] = False

        overlapped_gene_clusters = []
        for node in nodes:
            if not self._visited[node]:
                cluster_members = self.dfs(node, overlap_type)
                # we only needs clusters with >1 members, coz we are finding overlapped genes
                if len(cluster_members) > 1:
                    overlapped_gene_clusters.append(cluster_members)

        return overlapped_gene_clusters
        # (overall_stats, overlapped_gene_clusters, cluster_stats_array)
        # # analyze overall overlapped data (genes expression and sequences stats)
        # overall_stats = AnalyzerData()
        # for edge_info in self.graph_edges_info:
        #     if edge_info.edge_type != category: continue
        #
        #     seq = edge_info.shared_sequences[0]
        #     frame = edge_info.shared_sequences[1]
        #
        #     overall_stats.analyze_sequence_stats(seq, frame)
        # for cluster in overlapped_gene_clusters:
        #     for (gene_id, _) in cluster:
        #         overall_stats.analyze_expression_profile(self.gene_expressions[gene_id])
        #
        # # analyze individual clusters stats
        # cluster_stats_array = []
        # for cluster in overlapped_gene_clusters:
        #     cluster_stats = AnalyzerData()
        #     for edge_info in self.graph_edges_info:
        #         if edge_info.edge_type != category: continue
        #
        #         # if this edge is not connecting nodes from this cluster
        #         cluster_contains_edge_node = False
        #         for (gene_id, gene_expression) in cluster:
        #             if gene_id == edge_info.node1: cluster_contains_edge_node = True
        #         if not cluster_contains_edge_node: continue
        #
        #         seq = edge_info.shared_sequences[0]
        #         frame = edge_info.shared_sequences[1]
        #
        #         cluster_stats.analyze_sequence_stats(seq, frame)
        #     cluster_stats_array.append(cluster_stats)

        # (overall_stats,clusters,cluster_stats)
        # for each cluster find stats


class AnalyzerData:
    # region Class/Constant Variables and Instance Variables

    # Codons Degeneracy/alternatives values
    CODONS_DEGENERACY = {'L': 6, 'R': 6, 'S': 6,  # high
                         'A': 4, 'G': 4, 'P': 4, 'T': 3, 'V': 4, 'I': 4,  # medium
                         'C': 2, 'D': 2, 'E': 2, 'F': 2, 'H': 2, 'K': 2, 'N': 2, 'Q': 2, 'Y': 2, 'M': 1, 'W': 1}  # low

    # Amino acids by charge or polarity
    AMINO_ACID_TYPE = {'R': '+', 'H': '+', 'K': '+',
                       'D': '-', 'E': '-',
                       'S': 'polar', 'T': 'polar', 'Y': 'polar', 'N': 'polar', 'Q': 'polar',
                       'G': 'non-polar', 'A': 'non-polar', 'V': 'non-polar', 'C': 'non-polar', 'P': 'non-polar',
                       'L': 'non-polar', 'I': 'non-polar', 'M': 'non-polar', 'W': 'non-polar', 'F': 'non-polar',
                       }

    SUB_REGIONS = 50

    def __init__(self):
        self.analyzed_sequence_length = 0
        self.analyzed_expression_profiles = {}
        self.nucleotide_frequency = {'C': 0, 'G': 0, 'A': 0, 'T': 0}
        self.di_nucleotide_frequency = {'AA': 0, 'AT': 0, 'AG': 0, 'AC': 0,
                                        'TA': 0, 'TT': 0, 'TC': 0, 'TG': 0,
                                        'CA': 0, 'CT': 0, 'CC': 0, 'CG': 0,
                                        'GA': 0, 'GT': 0, 'GC': 0, 'GG': 0}
        self.amino_acid_frequency = {'L': 0, 'R': 0, 'S': 0, 'A': 0, 'G': 0, 'P': 0, 'T': 0, 'V': 0, 'I': 0, 'C': 0,
                                     'D': 0, 'E': 0, 'F': 0, 'H': 0, 'K': 0, 'N': 0, 'Q': 0, 'Y': 0, 'M': 0, 'W': 0,
                                     '*': 0}  # stop codon frequency
        self.degeneracy_frequency = {'high': 0, 'medium': 0, 'low': 0}
        self.amino_acid_types_frequency = {'+': 0, '-': 0, 'polar': 0, 'non-polar': 0}
        self.nucleotide_frequency_by_subregions = []
        self.degeneracy_frequency_by_subregions = []
        for i in range(0, self.SUB_REGIONS):
            self.nucleotide_frequency_by_subregions.append({'C': 0, 'G': 0, 'A': 0, 'T': 0})
            self.degeneracy_frequency_by_subregions.append({'high': 0, 'medium': 0, 'low': 0})

    @staticmethod
    def is_valid_nucleotide(char):
        return char == 'A' or char == 'T' or char == 'C' or char == 'G'

    @staticmethod
    def degeneracy2level(degeneracy):
        return 'high' if degeneracy > 4 else 'medium' if degeneracy > 2 else 'low'

    # endregion

    # region Main Methods

    def analyze_expression_profile(self, expressions_dict):
        for (key, value) in expressions_dict.items():
            if not self.analyzed_expression_profiles.__contains__(key):
                self.analyzed_expression_profiles[key] = 0
            self.analyzed_expression_profiles[key] += value

    # additive method, which analyzes sequence and adds stats to the object variables
    def analyze_sequence_stats(self, sequence: str, frame):
        assert frame != '-'
        frame = int(frame)  # make sure it is converted integer, if it's not

        self.analyzed_sequence_length += len(sequence)

        sequence = sequence.upper()
        unframed_seq = sequence[frame: len(sequence)]
        unframed_seq = unframed_seq[0: len(unframed_seq) - len(unframed_seq) % 3]

        self.__analyze_sequence_by_nucleotides(sequence)
        self.__analyze_sequence_by_amino_acids(unframed_seq)

    # updates nucleotide_frequency and di_nucleotide_frequency dictionary accordingly
    # also calculates nucleotides frequency in subdivided k region
    def __analyze_sequence_by_nucleotides(self, seq):
        # analyze overall nucleotide & di-nucleotide composition
        for i in range(0, len(seq)):
            if self.is_valid_nucleotide(seq[i]):
                self.nucleotide_frequency[seq[i]] += 1
            if i > 0:
                if self.is_valid_nucleotide(seq[i]) and self.is_valid_nucleotide(seq[i - 1]):
                    self.di_nucleotide_frequency[seq[i - 1:i + 1]] += 1

        # analyze nucleotide frequency by sub-regions
        sequence = AnalyzerData.make_sequence_k_dividable(seq, self.SUB_REGIONS)
        part_length = len(sequence) // self.SUB_REGIONS  # it is now 100% dividable by k
        for i in range(0, self.SUB_REGIONS):
            l = i * part_length
            r = l + part_length - 1
            for j in range(l, r + 1):
                if self.is_valid_nucleotide(sequence[j]):
                    self.nucleotide_frequency_by_subregions[i][sequence[j]] += 1

    # updates nucleotide_frequency and di_nucleotide_frequency dictionary accordingly
    def __analyze_sequence_by_amino_acids(self, seq):
        rec_sec = Seq(seq)
        amino_acids_seq = str(rec_sec.translate(table="Standard"))

        # analyze amino acid repeat
        # amino_acids_seq.replace('*','')
        max_repeated_seq = ''

        for i in range(0, len(amino_acids_seq)):
            used = {}
            cnt = 0
            for j in range(i, len(amino_acids_seq)):
                if not used.__contains__(amino_acids_seq[j]):
                    cnt += 1
                used[amino_acids_seq[j]] = 1
                if cnt > 3: break
                if len(max_repeated_seq) < j - i + 1:
                    max_repeated_seq = amino_acids_seq[i:j + 1]

        # print(max_repeated_seq)
        # if debug == 'gene:ENSG00000095713':
        #     print(amino_acids_seq)
        # if len(max_repeated_seq) >= 10:
        #     print(f'repeat: {max_repeated_seq}  gene: {debug}')
        #
        # for repeat_length in range(1, len(amino_acids_seq)):
        #     l = 0
        #     while l + repeat_length <= len(amino_acids_seq):
        #         my_seq = amino_acids_seq[l:l + repeat_length]
        #         r = l + repeat_length
        #         while r + repeat_length <= len(amino_acids_seq):
        #             if amino_acids_seq[r:r + repeat_length] == my_seq:
        #                 r = r + repeat_length
        #             else:
        #                 break
        #         if r == l + repeat_length: r = l
        #
        #         if r - l + 1 > len(max_repeated_seq):
        #             max_repeated_seq = my_seq[l:r]
        #
        #         l += 1

        # analyze overall amino acid frequency
        for i in range(0, len(amino_acids_seq)):
            amino_acid = amino_acids_seq[i]
            if not self.amino_acid_frequency.__contains__(amino_acid): continue
            self.amino_acid_frequency[amino_acid] += 1
            if amino_acid != '*':
                level = self.degeneracy2level(self.CODONS_DEGENERACY[amino_acid])
                self.degeneracy_frequency[level] += 1
                self.amino_acid_types_frequency[self.AMINO_ACID_TYPE[amino_acid]] += 1

        # analyze amino acid high-degeneracy frequency by parts
        amino_acids_seq = AnalyzerData.make_sequence_k_dividable(amino_acids_seq, self.SUB_REGIONS)
        part_length = len(amino_acids_seq) // self.SUB_REGIONS  # it is now 100% dividable by k
        for i in range(0, self.SUB_REGIONS):
            l = i * part_length
            r = l + part_length - 1
            for j in range(l, r + 1):
                if not self.CODONS_DEGENERACY.__contains__(amino_acids_seq[i]): continue
                level = self.degeneracy2level(self.CODONS_DEGENERACY[amino_acids_seq[i]])
                self.degeneracy_frequency_by_subregions[i][level] += 1

    # endregion

    # region Static methods
    @staticmethod
    def make_sequence_k_dividable(sequence, k):
        # sliding window with size k.
        if len(sequence) % k != 0:  # place equally distanced '*'s to make it k-dividable for later smooth analyze
            total_gaps_to_place = k - (len(sequence) % k)
            gaps_placed = 0
            new_seq = ""
            length = len(sequence)
            for index in range(length):
                new_seq += sequence[index]
                # index / len * gaps_to_place = gaps_placed   index*(gaps_to_place+1) // len= 15
                if round((index + 1) * total_gaps_to_place / length) > gaps_placed:
                    n = round((index + 1) * total_gaps_to_place / length) - gaps_placed
                    new_seq += '*' * n  # unknown
                    gaps_placed += n
            assert len(new_seq) % k == 0  # just, make sure
            return new_seq
        return sequence

    # endregion

    # region Get Data  (+ Formatters) in percents
    #
    # def get_mean_expression_profile(self):
    #     return self.analyzed_expression_profiles

    def get_debug_gc_content(self):
        return f"{self.nucleotide_frequency['C']} {self.nucleotide_frequency['G']}"

    def get_mean_gc_content(self):
        total = self.nucleotide_frequency['C'] + self.nucleotide_frequency['G'] + self.nucleotide_frequency['A'] + \
                self.nucleotide_frequency['T']
        gc = self.nucleotide_frequency['C'] + self.nucleotide_frequency['G']
        return -1 if total == 0 else gc * 100 / total

    def get_mean_hd_content(self):
        total = self.degeneracy_frequency['low'] + self.degeneracy_frequency['medium'] + self.degeneracy_frequency[
            'high']
        hd = self.degeneracy_frequency['high']
        return -1 if total == 0 else hd * 100 / total

    def get_regional_gc_contents(self):
        output = []
        sum = 0
        my = 0
        for i in range(0, self.SUB_REGIONS):
            gc_content = self.nucleotide_frequency_by_subregions[i]['C'] + self.nucleotide_frequency_by_subregions[i][
                'G']
            total_sum = self.nucleotide_frequency_by_subregions[i]['A'] + self.nucleotide_frequency_by_subregions[i][
                'T'] + self.nucleotide_frequency_by_subregions[i]['C'] + self.nucleotide_frequency_by_subregions[i]['G']
            sum += total_sum
            my += gc_content
            value = 0 if total_sum == 0 else gc_content * 100 / total_sum
            output.append(value)
        return output

    def get_regional_hd_contents(self):
        output = []
        for i in range(0, self.SUB_REGIONS):
            hd_content = self.degeneracy_frequency_by_subregions[i]['high']
            total_sum = self.degeneracy_frequency_by_subregions[i]['low'] + self.degeneracy_frequency_by_subregions[i][
                'high'] + self.degeneracy_frequency_by_subregions[i]['medium']

            value = 0 if total_sum == 0 else hd_content * 100 / total_sum
            output.append(value)
        return output

    def beauty_dict(self, data_storage, use_percents: bool) -> list:
        sorted_dict = sorted(data_storage.items(), key=lambda kv: kv[1], reverse=True)
        if not use_percents:
            # just use only 2 point decimal
            for i in range(len(sorted_dict)):
                sorted_dict[i] = (sorted_dict[i][0], float("{:.2f}".format(sorted_dict[i][1])))
            return sorted_dict

        total_sum = 0
        for i in range(len(sorted_dict)):
            total_sum += sorted_dict[i][1]
        for i in range(len(sorted_dict)):
            value = 0 if total_sum == 0 else float("{:.2f}".format(sorted_dict[i][1] * 100 / total_sum))
            sorted_dict[i] = (sorted_dict[i][0], value)
        return sorted_dict

    def get_formatted_expression_profile(self, indent_tabs=0, use_percents=False):
        expression_profile = self.beauty_dict(self.analyzed_expression_profiles, use_percents)

        indent = '\t' * indent_tabs
        percent_symbol = '%' if use_percents else ''
        results = f"{indent}expression profile [nTPM]:\n"

        index = 0
        for (key, value) in expression_profile:
            key = key.replace('Tissue RNA - ', '').replace(' [nTPM]', '')
            index += 1
            end_symbols = '' if index % 3 != 0 else '\n'
            front_symbols = indent if index % 3 == 1 else ''
            results += f"{front_symbols}\t{key} - {value}{percent_symbol}{end_symbols}"
        results += "\n\n"

        return results

    def get_formatted_sequence_stats(self, indent_tabs=0, use_percents=False):
        nucleotide_frequency = self.beauty_dict(self.nucleotide_frequency, use_percents)
        di_nucleotide_frequency = self.beauty_dict(self.di_nucleotide_frequency, use_percents)
        amino_acid_frequency = self.beauty_dict(self.amino_acid_frequency, use_percents)
        amino_acid_degeneracy_frequency = self.beauty_dict(self.degeneracy_frequency, use_percents)
        amino_acid_type_frequency = self.beauty_dict(self.amino_acid_types_frequency, use_percents)

        percent_symbol = '%' if use_percents else ''
        indent = '\t' * indent_tabs

        results = f"{indent}nucleotides:\n"
        index = 0
        for (key, value) in nucleotide_frequency:
            index += 1
            end_symbols = '' if index % 4 != 0 else '\n'
            front_symbols = indent if index % 4 == 1 else ''
            results += f"{front_symbols}\t{key} - {value}{percent_symbol}{end_symbols}"

        results += f"\n\n{indent}di-nucleotides:\n"
        index = 0
        for (key, value) in di_nucleotide_frequency:
            index += 1
            end_symbols = '' if index % 4 != 0 else '\n'
            front_symbols = indent if index % 4 == 1 else ''
            results += f"{front_symbols}\t{key} - {value}{percent_symbol}{end_symbols}"

        results += f"\n\n{indent}amino-acids:\n"
        index = 0
        for (key, value) in amino_acid_frequency:
            index += 1
            end_symbols = '' if index % 4 != 0 else '\n'
            front_symbols = indent if index % 4 == 1 else ''
            results += f"{front_symbols}\t{key} - {value}{percent_symbol}{end_symbols}"

        results += f"\n\n{indent}amino-acids by codon degeneracy:\n"
        index = 0
        for (key, value) in amino_acid_degeneracy_frequency:
            index += 1
            front_symbols = indent if index % 4 == 1 else ''
            results += f"{front_symbols}\t{key} - {value}{percent_symbol}"

        results += f"\n\n{indent}amino-acids by type:\n"
        index = 0
        for (key, value) in amino_acid_type_frequency:
            index += 1
            front_symbols = indent if index % 4 == 1 else ''
            results += f"{front_symbols}\t{key} - {value}{percent_symbol}"

        results += "\n\n"

        return results
    # endregion
