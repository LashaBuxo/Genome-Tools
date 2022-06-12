from Bio.Seq import Seq
import numpy as np


# (mean, std, median, interquartile)
def get_value_estimation(values, with_sum_max=False, use_percent_mark=False):
    if len(values) == 0:
        if with_sum_max: return 0, 0, 0, 0, 0, 0
        return 0, 0, 0, 0
    mean_value = float('{:.2f}'.format(np.mean(values)))
    std_value = float('{:.2f}'.format(np.std(values)))
    median_q1, median, median_q3 = np.percentile(values, [25, 50, 75])
    median = float('{:.2f}'.format(median))
    interquartile = float('{:.2f}'.format(median_q3 - median_q1))

    if with_sum_max:
        sum_value = float('{:.2f}'.format(np.sum(values)))
        max_value = float('{:.2f}'.format(np.max(values)))
        if use_percent_mark:
            return f'{mean_value}%', f'{std_value}%', f'{median}%', f'{interquartile}%', f'{sum_value}%', f'{max_value}%'
        return mean_value, std_value, median, interquartile, sum_value, max_value
    else:
        if use_percent_mark:
            return f'{mean_value}%', f'{std_value}%', f'{median}%', f'{interquartile}%'
        return mean_value, std_value, median, interquartile


class AnalyzerGraph:
    class GraphEdge:
        def __init__(self, gene1_id, gene2_id, edge_type=None):
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
            if (allowed_edge_type is None or (edge_type is allowed_edge_type)) and not self._visited[neighbor]:
                cluster_members += self.dfs(neighbor, allowed_edge_type)
        return cluster_members

    def get_connected_clusters(self, overlap_interaction=None):
        nodes = self.neighbors.keys()

        for node in nodes:
            self._visited[node] = False

        overlapped_gene_clusters = []
        for node in nodes:
            if not self._visited[node]:
                cluster_members = self.dfs(node, overlap_interaction)
                # we only needs clusters with >1 members, coz we are finding overlapped genes
                if len(cluster_members) > 1:
                    overlapped_gene_clusters.append(cluster_members)

        return overlapped_gene_clusters


class AnalyzerData:
    # region Class/Constant Variables and Instance Variables

    SUB_REGIONS = 50

    def __init__(self):
        self.analyzed_peptide = ''
        self.analyzed_sequences = ''
        self.total_analyzed_sequence_length = 0
        self.max_analyzed_sequence_length = 0
        self.nucleotide_frequency = {'C': 0, 'G': 0, 'A': 0, 'T': 0}

        self.nucleotide_frequency_by_subregions = []

        for i in range(0, self.SUB_REGIONS):
            self.nucleotide_frequency_by_subregions.append({'C': 0, 'G': 0, 'A': 0, 'T': 0})

    @staticmethod
    def is_valid_nucleotide(char):
        return char == 'A' or char == 'T' or char == 'C' or char == 'G'

    # endregion

    # region Main Methods

    # additive method, which analyzes sequence and adds stats to the object variables
    def analyze_sequence_stats(self, seq: str, frame, peptide_indexes=None):
        assert frame != '-'
        frame = int(frame)  # make sure it is converted integer, if it's not

        seq = seq.upper()
        unframed_seq = seq[frame: len(seq)]
        unframed_seq = unframed_seq[0: len(unframed_seq) - len(unframed_seq) % 3]

        self.__analyze_sequence_by_amino_acids(unframed_seq, peptide_indexes)
        self.__analyze_sequence_by_nucleotides(seq, frame)

    # updates nucleotide_frequency and di_nucleotide_frequency dictionary accordingly
    # also calculates nucleotides frequency in subdivided k region
    def __analyze_sequence_by_nucleotides(self, seq, frame):
        self.total_analyzed_sequence_length += len(seq)
        self.max_analyzed_sequence_length = max(self.max_analyzed_sequence_length, len(seq))
        self.analyzed_sequences += f"{'...' if len(self.analyzed_sequences) > 0 else ''}{seq if frame == 0 else seq[0:frame] + '|' + seq[frame: len(seq)]}"

        # analyze overall nucleotide & di-nucleotide composition
        for i in range(0, len(seq)):
            if self.is_valid_nucleotide(seq[i]):
                self.nucleotide_frequency[seq[i]] += 1

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
    def __analyze_sequence_by_amino_acids(self, seq, peptide_indexes=None):
        rec_sec = Seq(seq)
        amino_acids_seq = str(rec_sec.translate(table="Standard"))

        formatted_aa = amino_acids_seq if peptide_indexes is None else f'({peptide_indexes[0]}-{peptide_indexes[1]}){amino_acids_seq}'
        self.analyzed_peptide += f"{'...' if len(self.analyzed_sequences) > 0 else ''}{formatted_aa}"

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

    def get_gc_content(self):
        total = self.nucleotide_frequency['C'] + self.nucleotide_frequency['G'] + self.nucleotide_frequency['A'] + \
                self.nucleotide_frequency['T']
        gc = self.nucleotide_frequency['C'] + self.nucleotide_frequency['G']
        return -1 if total == 0 else gc / total

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
            value = 0 if total_sum == 0 else gc_content / total_sum
            output.append(value)
        return output

    def beauty_dict(self, data_storage) -> list:
        sorted_dict = sorted(data_storage.items(), key=lambda kv: kv[1], reverse=True)

        total_sum = 0
        for i in range(len(sorted_dict)):
            total_sum += sorted_dict[i][1]

        for i in range(len(sorted_dict)):
            value = 0 if total_sum == 0 else float("{:.3f}".format(sorted_dict[i][1] / total_sum))
            sorted_dict[i] = (sorted_dict[i][0], value)
        return sorted_dict

    def get_short_sequence_stats(self):
        # GC - 61.36%
        GC = "{:.2f}".format(self.get_gc_content()) + '%'
        result = f'GC - {GC}'
        return result

    # endregion
