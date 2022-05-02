from Bio.Seq import Seq


class AnalyzerData:
    # region Variables and Tables
    nucleotide_frequency = {'C': 0, 'G': 0, 'A': 0, 'T': 0}

    di_nucleotide_frequency = {'AA': 0, 'AT': 0, 'AG': 0, 'AC': 0,
                               'TA': 0, 'TT': 0, 'TC': 0, 'TG': 0,
                               'CA': 0, 'CT': 0, 'CC': 0, 'CG': 0,
                               'GA': 0, 'GT': 0, 'GC': 0, 'GG': 0}

    amino_acid_frequency = {'L': 0, 'R': 0, 'S': 0, 'A': 0, 'G': 0, 'P': 0, 'T': 0, 'V': 0, 'I': 0, 'C': 0,
                            'D': 0, 'E': 0, 'F': 0, 'H': 0, 'K': 0, 'N': 0, 'Q': 0, 'Y': 0, 'M': 0, 'W': 0,
                            '*': 0}  # stop codon frequency

    degeneracy_frequency = {'high': 0, 'medium': 0, 'low': 0}

    # Codons Degeneracy/alternatives values
    CODONS_DEGENERACY = {'L': 6, 'R': 6, 'S': 6,  # high
                         'A': 4, 'G': 4, 'P': 4, 'T': 3, 'V': 4, 'I': 4,  # medium
                         'C': 2, 'D': 2, 'E': 2, 'F': 2, 'H': 2, 'K': 2, 'N': 2, 'Q': 2, 'Y': 2, 'M': 1, 'W': 1}  # low

    SUB_REGIONS = 50

    nucleotide_frequency_by_subregions = [None] * SUB_REGIONS

    degeneracy_frequency_by_subregions = [None] * SUB_REGIONS

    def release_memory(self):
        for key in self.nucleotide_frequency.keys():
            self.nucleotide_frequency[key] = 0
        for key in self.di_nucleotide_frequency.keys():
            self.di_nucleotide_frequency[key] = 0
        for key in self.amino_acid_frequency.keys():
            self.amino_acid_frequency[key] = 0
        for key in self.degeneracy_frequency.keys():
            self.degeneracy_frequency[key] = 0
        for i in range(0, self.SUB_REGIONS):
            self.nucleotide_frequency_by_subregions[i] = 0
            self.degeneracy_frequency_by_subregions[i] = 0

    def __init__(self):
        for i in range(0, self.SUB_REGIONS):
            self.nucleotide_frequency_by_subregions[i] = {'C': 0, 'G': 0, 'A': 0, 'T': 0}
            self.degeneracy_frequency_by_subregions[i] = {'high': 0, 'medium': 0, 'low': 0}

    @staticmethod
    def is_valid_nucleotide(char):
        return char == 'A' or char == 'T' or char == 'C' or char == 'G'

    @staticmethod
    def degeneracy2level(degeneracy):
        return 'high' if degeneracy > 4 else 'medium' if degeneracy > 2 else 'low'

    # endregion

    # region Main Methods
    # additive method, which analyzes sequence and adds stats to the object variables
    def analyze_sequence_stats(self, sequence: str, frame):
        assert frame != '-'
        frame = int(frame)  # make sure it is converted integer, if it's not

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
        my_seq = Seq(seq)
        amino_acids_seq = my_seq.translate(table="Standard")

        # analyze overall amino acid frequency
        for i in range(0, len(amino_acids_seq)):
            if not self.amino_acid_frequency.__contains__(amino_acids_seq[i]): continue
            self.amino_acid_frequency[amino_acids_seq[i]] += 1
            if amino_acids_seq[i] != '*':
                level = self.degeneracy2level(self.CODONS_DEGENERACY[amino_acids_seq[i]])
                self.degeneracy_frequency[level] += 1

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

    def get_mean_gc_content(self):
        total = self.nucleotide_frequency['C'] + self.nucleotide_frequency['G'] + self.nucleotide_frequency['A'] + \
                self.nucleotide_frequency['T']
        gc = self.nucleotide_frequency['C'] + self.nucleotide_frequency['G']
        return gc * 100 / total

    def get_mean_hd_content(self):
        total = self.degeneracy_frequency['low'] + self.degeneracy_frequency['medium'] + self.degeneracy_frequency[
            'high']
        hd = self.degeneracy_frequency['high']
        return hd * 100 / total

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
            value = gc_content * 100 / total_sum
            output.append(value)
        return output

    def get_regional_hd_contents(self):
        output = []
        for i in range(0, self.SUB_REGIONS):
            hd_content = self.degeneracy_frequency_by_subregions[i]['high']
            total_sum = self.degeneracy_frequency_by_subregions[i]['low'] + self.degeneracy_frequency_by_subregions[i][
                'high'] + self.degeneracy_frequency_by_subregions[i]['medium']

            value = hd_content * 100 / total_sum
            output.append(value)
        return output

    def beauty_dict(self, data_storage, use_percents: bool) -> list:
        sorted_dict = sorted(data_storage.items(), key=lambda kv: kv[1], reverse=True)
        if not use_percents:
            return sorted_dict

        total_sum = 0
        for i in range(len(sorted_dict)):
            total_sum += sorted_dict[i][1]
        for i in range(len(sorted_dict)):
            value = float("{:.2f}".format(sorted_dict[i][1] * 100 / total_sum))
            sorted_dict[i] = (sorted_dict[i][0], value)
        return sorted_dict

    def get_formatted_results(self, indent_tabs=0, use_percents=False):
        nucleotide_frequency = self.beauty_dict(self.nucleotide_frequency, use_percents)
        di_nucleotide_frequency = self.beauty_dict(self.di_nucleotide_frequency, use_percents)
        amino_acid_frequency = self.beauty_dict(self.amino_acid_frequency, use_percents)
        amino_acid_degeneracy_frequency = self.beauty_dict(self.degeneracy_frequency, use_percents)

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

        results += "\n"

        return results
    # endregion
