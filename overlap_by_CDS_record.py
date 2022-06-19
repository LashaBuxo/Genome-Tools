from worker_genome import *
from worker_analyzer import AnalyzerData
from os.path import exists
import requests, sys


class OverlapRecord:
    FLANK_LENGTH = 50

    def __init__(self, transcript1, transcript2, chr_id, genome_ref: GenomeWorker, overlapped_segments):
        # MAKE SURE 1st is positive strand if there is diff strand overlap

        self.transcript1 = transcript1
        self.transcript2 = transcript2
        self.gene1 = genome_ref.get_transcript_parent(transcript1.id)
        self.gene2 = genome_ref.get_transcript_parent(transcript2.id)

        self.chr_id = chr_id
        self.genome_ref = genome_ref
        self.overlapped_segments = overlapped_segments

        # calculate initial stats without waiting 'analyze' command
        self.gene1_symbol = genome_ref.get_gene_symbol(self.gene1)
        self.gene2_symbol = genome_ref.get_gene_symbol(self.gene2)

        fragments1 = self.genome_ref.get_fragments_from_transcript(self.transcript1.id)
        fragments2 = self.genome_ref.get_fragments_from_transcript(self.transcript2.id)

        self.transcriptic_overlap_type = genome_ref.get_overlap_type(self.transcript1, self.transcript2)

        exonic_overlap_types = []
        for frag1 in fragments1:
            for frag2 in fragments2:
                if frag1.featuretype != 'CDS' or frag2.featuretype != 'CDS': continue
                if genome_ref.get_overlap_type(frag1, frag2) != OVERLAP_TYPE.NONE:
                    exonic_overlap_types.append(genome_ref.get_overlap_type(frag1, frag2))

        flag = False
        for index in range(1, len(exonic_overlap_types)):
            if exonic_overlap_types[index] != exonic_overlap_types[0]: flag = True

        self.exonic_overlap_type = OVERLAP_TYPE.MULTI if flag else exonic_overlap_types[0]

        self.genomic_overlap_type = genome_ref.get_overlap_type(self.gene1, self.gene2)
        self.is_fully_analyzed = False

        # declare parameters for later analyses
        self.exons_outline1 = ''
        self.exons_outline2 = ''
        self.flanks_positive = []
        self.flanks_negative = []
        self.seq1_stats = AnalyzerData()
        self.seq2_stats = AnalyzerData()

        self.max_record_length_in_siblings = -1

    def get_overlapped_sequence(self):
        if not self.is_fully_analyzed: self.analyze_record()
        if self.gene1.strand != self.gene2.strand:
            return self.seq1_stats.analyzed_sequences if self.gene1.strand == '+' else self.seq2_stats.analyzed_sequences
        return self.seq1_stats.analyzed_sequences if self.gene1.end - self.gene1.start > self.gene2.end - self.gene2.start else self.seq2_stats.analyzed_sequences

    def get_transcript_conservation_score(self, transcript_id):
        return self.genome_ref.get_transcript_conservation_info(transcript_id)

    def get_record_conservation_score(self):
        if not self.is_fully_analyzed: self.analyze_record()
        return min(self.get_transcript_conservation_score(self.transcript1.id),
                   self.get_transcript_conservation_score(self.transcript2.id))

    def get_record_length(self, max_in_siblings=False):
        if not self.is_fully_analyzed: self.analyze_record()

        if max_in_siblings:
            if self.max_record_length_in_siblings != -1:
                return self.max_record_length_in_siblings
            else:
                return self.seq1_stats.total_analyzed_sequence_length

        return self.seq1_stats.total_analyzed_sequence_length

    def get_record_GC_content(self):
        if not self.is_fully_analyzed: self.analyze_record()
        return self.seq1_stats.get_gc_content()

    def is_record_diff_strand(self):
        return self.transcriptic_overlap_type == OVERLAP_TYPE.DIVERGENT or \
               self.transcriptic_overlap_type == OVERLAP_TYPE.CONVERGENT or self.transcriptic_overlap_type == OVERLAP_TYPE.DIFF_NESTED

    def is_record_about_gene(self, g_id):
        return g_id == self.gene1.id or g_id == self.gene2.id

    def analyze_record(self):
        fragments1 = self.genome_ref.get_fragments_from_transcript(self.transcript1.id)
        fragments2 = self.genome_ref.get_fragments_from_transcript(self.transcript2.id)

        self.flanks_positive = self._calculate_flanks('+')
        self.flanks_negative = self._calculate_flanks('-')

        self.exons_outline1 = self._calculate_exons_outline(fragments1, self.transcript1.strand)
        self.exons_outline2 = self._calculate_exons_outline(fragments2, self.transcript2.strand)

        self.seq1_stats = self._calculate_sequence_stats(self.transcript1, 0)
        self.seq2_stats = self._calculate_sequence_stats(self.transcript2, 1)

        self.is_fully_analyzed = True

    def _calculate_flanks(self, strand):
        flanks = []
        for interval, frames in self.overlapped_segments:
            l, r = interval[0], interval[1]
            flank_pos_1 = self.genome_ref.retrieve_interval_sequence(self.chr_id, l - self.FLANK_LENGTH, l - 1, strand)
            flank_pos_2 = self.genome_ref.retrieve_interval_sequence(self.chr_id, r + 1, r + self.FLANK_LENGTH, strand)
            if strand == '-':
                flanks.append((flank_pos_1[::-1], flank_pos_2[::-1]))
            else:
                flanks.append((flank_pos_1, flank_pos_2))
        return flanks

    def _calculate_exons_outline(self, frags, strand):
        outline = ''
        for frag in frags:
            if frag.featuretype == 'three_prime_UTR' or frag.featuretype == 'five_prime_UTR':
                outline += 'u'
            if frag.featuretype == 'CDS':
                is_in_overlap = False
                for interval, _ in self.overlapped_segments:
                    if self.genome_ref.are_segments_overlapped(interval, (frag.start, frag.end)):
                        is_in_overlap = True
                outline += 'X' if is_in_overlap else 'x'
        return outline if strand == '+' else outline[::-1]

    def _calculate_sequence_stats(self, transcript, frame_index):
        # sequence/peptide compositions,   functional hits (swissprot)
        stats = AnalyzerData()
        for interval, frames in self.overlapped_segments:
            l, r = interval[0], interval[1]
            frame = frames[frame_index]
            seq = self.genome_ref.retrieve_interval_sequence(self.chr_id, l, r, transcript.strand)

            peptide_start, peptide_end = self.__get_peptide_range(transcript, (l, r))
            stats.analyze_sequence_stats(seq, frame, peptide_indexes=(peptide_start, peptide_end))
        return stats

    def __get_peptide_range(self, transcript, cds_interval):
        frags = self.genome_ref.get_fragments_from_transcript(transcript.id)
        l, r = cds_interval[0], cds_interval[1]
        cds = []
        for i in range(0, len(frags)):
            frag = frags[i]
            if frag.featuretype == 'CDS':
                for j in range(frag.start, frag.end + 1):
                    cds.append(j)
        if transcript.strand == '+':
            return cds.index(l) // 3 + 1, cds.index(r) // 3 + 1
        else:
            cds.reverse()
            return cds.index(r) // 3 + 1, cds.index(l) // 3 + 1

    def get_formatted_stats(self, indent_tabs=0):
        if not self.is_fully_analyzed: self.analyze_record()

        indent = '\t' * indent_tabs

        results = f'{indent}{self.transcript1.id} ({self.gene1_symbol}) - {self.transcript2.id} ({self.gene2_symbol}):\n'

        results += f"{indent}\tspecifics: transcriptic overlap type - {self.transcriptic_overlap_type.short_name()} | exonic overlap type - {self.exonic_overlap_type.short_name()} | len - {self.seq1_stats.total_analyzed_sequence_length} " \
                   f"| exons - [5'-{self.exons_outline1}-3'], [5'-{self.exons_outline2}-3']\n"
        results += f'{indent}\tflanks(+): '
        for flank in self.flanks_positive:
            results += f'[{flank[0]}]-[overlap]-[{flank[1]}] '
        results += f'\n{indent}\tflanks(-): '
        for flank in self.flanks_negative:
            results += f'[{flank[0]}]-[overlap]-[{flank[1]}] '

        cons_score1 = self.get_transcript_conservation_score(self.transcript1.id)
        cons_score2 = self.get_transcript_conservation_score(self.transcript2.id)
        results += f"\n{indent}\trecord conservation score: {self.get_record_conservation_score()}  [{self.gene1_symbol} - {cons_score1}] [{self.gene2_symbol} - {cons_score2}]\n"

        results += f"{indent}\tseq 1: 5'...{self.seq1_stats.analyzed_sequences}...3'\n"
        results += f"{indent}\t\t\t{self.seq1_stats.get_short_sequence_stats()}\n"

        results += f"{indent}\tseq 2: 5'...{self.seq2_stats.analyzed_sequences}...3'\n"
        results += f"{indent}\t\t\t{self.seq2_stats.get_short_sequence_stats()}\n"

        ###
        results += f"{indent}\tpeptide 1: N...{self.seq1_stats.analyzed_peptide}...C\n"
        results += f"{indent}\tpeptide 2: N...{self.seq2_stats.analyzed_peptide}...C\n"

        return results
