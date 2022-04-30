from Bio.Seq import Seq


class AnalyzerData:
    # Codons Degeneracy/alternatives values
    CODONS_DEGENERACY = {'L': 6, 'R': 6, 'S': 6,  # high
                         'A': 4, 'G': 4, 'P': 4, 'T': 3, 'V': 4, 'I': 4,  # medium
                         'C': 2, 'D': 2, 'E': 2, 'F': 2, 'H': 2, 'K': 2, 'N': 2, 'Q': 2, 'Y': 2, 'M': 1, 'W': 1}  # low

    @staticmethod
    def degeneracy2level(degeneracy):
        return 'high' if degeneracy > 4 else 'medium' if degeneracy > 2 else 'low'

    # region Variables
    nucleotide_frequency = {'C': 0, 'G': 0, 'A': 0, 'T': 0, '*': 0}

    di_nucleotide_frequency = {'AA': 0, 'AT': 0, 'AG': 0, 'AC': 0,
                               'TA': 0, 'TT': 0, 'TC': 0, 'TG': 0,
                               'CA': 0, 'CT': 0, 'CC': 0, 'CG': 0,
                               'GA': 0, 'GT': 0, 'GC': 0, 'GG': 0,
                               '**': 0}

    amino_acid_frequency = {'L': 0, 'R': 0, 'S': 0, 'A': 0, 'G': 0, 'P': 0, 'T': 0, 'V': 0, 'I': 0, 'C': 0,
                            'D': 0, 'E': 0, 'F': 0, 'H': 0, 'K': 0, 'N': 0, 'Q': 0, 'Y': 0, 'M': 0, 'W': 0,
                            '*': 0, 'X': 0}  # stop codon frequency

    amino_acid_degeneracy_frequency = {'high': 0, 'medium': 0, 'low': 0}

    # endregion

    # additive method, which analyzes sequence and adds stats to the object variables
    def analyze_sequence_stats(self, sequence: str, frame, debug_value=""):
        assert frame != '-'
        frame = int(frame)  # make sure it is converted integer, if it's not

        sequence = sequence.upper()
        unframed_seq = sequence[frame: len(sequence)]
        unframed_seq = unframed_seq[0: len(unframed_seq) - len(unframed_seq) % 3]

        self.__analyze_sequence_by_nucleotides(sequence)
        self.__analyze_sequence_by_amino_acids(unframed_seq)

    # updates nucleotide_frequency and di_nucleotide_frequency dictionary accordingly
    def __analyze_sequence_by_nucleotides(self, seq):
        for i in range(0, len(seq)):
            base = seq[i] if self.nucleotide_frequency.__contains__(seq[i]) else '*'
            self.nucleotide_frequency[base] += 1

            if i > 0:
                di_nucleotide = seq[i - 1:i + 1] if self.di_nucleotide_frequency.__contains__(
                    seq[i - 1:i + 1]) else '**'
                self.di_nucleotide_frequency[di_nucleotide] += 1

    # updates nucleotide_frequency and di_nucleotide_frequency dictionary accordingly
    def __analyze_sequence_by_amino_acids(self, seq):
        my_seq = Seq(seq)
        amino_acids_seq = my_seq.translate(table="Standard")
        for i in range(0, len(amino_acids_seq)):
            self.amino_acid_frequency[amino_acids_seq[i]] += 1
            if self.CODONS_DEGENERACY.__contains__(amino_acids_seq[i]):
                level = self.degeneracy2level(self.CODONS_DEGENERACY[amino_acids_seq[i]])
                self.amino_acid_degeneracy_frequency[level] += 1

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
        amino_acid_degeneracy_frequency = self.beauty_dict(self.amino_acid_degeneracy_frequency, use_percents)

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
