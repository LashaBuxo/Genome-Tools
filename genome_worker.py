import random
from os.path import exists

import gffutils
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from gffutils import Feature

from genome_worker_enums import *
from genome_worker_values import *


# Author: "Lasha Bukhnikashvili"
#
# Common Methods for working with genome data:
#   1) Creating database based on NCBI, Ensembl genome annotation
#   2) Import annotation features from created database
#   3) Filter specific features
#   4) calculate overlapping segments
#   and so on...
#
# Other scripts, which does specific calculations on genome use this script
# for working with Ensembl or NCBI annotation (with corresponding sequence data if necessary)

# Specific tools for working with Ensembl annotation and with sequence data
class GenomeWorker:
    def __init__(self, species: SPECIES, annotation_source: ANNOTATIONS, annotation_load_type: ANNOTATION_LOAD,
                 sequence_load_type: SEQUENCE_LOAD):

        assert annotation_source == ANNOTATIONS.ENSEMBL or annotation_source == ANNOTATIONS.NCBI

        print("loading data for species: " + str(species))

        self.species = species
        self.__genes_on_chr = [None] * (self.chromosomes_count() + 1)
        self.__sequences = [None] * (self.chromosomes_count() + 1)

        self.annotation_source = annotation_source
        self.annotation_load_type = annotation_load_type
        self.sequence_load_type = sequence_load_type

        self.__load_requested_sequences()
        self.__load_requested_features()

    # region Variables and Retrievers

    # gene features clustered by chromosome, they are on.
    # Example: genes_on_chr[22] = ['gene' features on chromosome 22]
    __genes_on_chr = []

    # {chr_id,DNA sequence of chromosome with id chr_id}
    __sequences = []

    #   {gene_id, [list of mRNA features whose parent is gene_id]}
    __gene_transcripts = {}

    #   {transcript_id, [list of 'CDS' and/or 'UTR(3'/5')' and/or 'exons' features whose parent is mRNA_transcript_id]}
    __transcript_fragments = {}

    def chromosomes_count(self):
        return NUMBER_OF_CHROMOSOMES[self.species.value]

    def genes_count_on_chr(self, chr_id):
        assert chr_id <= self.chromosomes_count()
        return len(self.__genes_on_chr[chr_id])

    def gene_by_ind(self, chr_id, index) -> Feature:
        assert chr_id <= self.chromosomes_count()
        return self.__genes_on_chr[chr_id][index]

    def __get_annotation_file_path(self):
        if self.annotation_source != ANNOTATIONS.ENSEMBL:
            return NCBI_ANNOTATIONS[self.species.value]
        return ENSEMBL_ANNOTATIONS[self.species.value]

    def __get_annotation_db_path(self):
        added_name = self.__get_annotation_file_path().replace('/', '_')
        generated_path = WORKING_DATABASES_DIRECTORY + self.species.name + "_" + self.annotation_source.name + "_" + added_name

        return generated_path + '.db'

    def __get_sequence_file_path(self):
        if self.annotation_source == ANNOTATIONS.ENSEMBL:
            return ENSEMBL_SEQUENCES[self.species.value]
        return NCBI_SEQUENCES[self.species.value]

    # endregion

    # region Features (from annotation) and Sequences (from assembly) load methods

    def __load_requested_sequences(self):
        if self.sequence_load_type == SEQUENCE_LOAD.LOAD:
            sequence_file_path = self.__get_sequence_file_path()
            for record in SeqIO.parse(sequence_file_path, 'fasta'):
                chr_id = self.chr_id_from_seq_id(self.annotation_source, record.id)
                if chr_id == -1: continue
                self.__sequences[chr_id] = record
            print("Using " + str(self.annotation_source) + ": sequence data loaded successfully!")

    # if database does not exists for specific chromosome, then
    # builds database from annotation file and fills/loads all
    # necessary features from the database.
    #
    # if database exists, then fills/loads all necessary data structures

    def __load_requested_features(self):
        if not exists(self.__get_annotation_db_path()):
            print("Creating database for " + str(self.annotation_source) + " only for first RUN it takes that long!")
            gffutils.create_db(self.__get_annotation_file_path(),
                               dbfn=self.__get_annotation_db_path(),
                               verbose=True, force=False, keep_order=False, merge_strategy='create_unique',
                               sort_attribute_values=False)

        features_db = gffutils.FeatureDB(self.__get_annotation_db_path(), keep_order=False)

        features_generator = features_db.features_of_type('gene')
        feature_genes = list(features_generator)

        # choose genes who has protein_coding attribute
        for gene in feature_genes:
            if self.annotation_source == ANNOTATIONS.NCBI:
                attribute_filter = 'gene_biotype'
            else:
                attribute_filter = 'biotype'
            attribute_value = 'protein_coding'
            if gene.attributes.__contains__(attribute_filter):
                if len(gene.attributes[attribute_filter]) != 1:
                    print("Lasha some issue found!")
                else:
                    # bio type must be 1 from these list
                    #  'miRNA', 'lncRNA', 'protein_coding', 'snoRNA', 'antisense_RNA', 'snRNA', 'ncRNA', 'tRNA',
                    #  'V_segment', 'misc_RNA', 'rRNA', 'other', 'C_region', 'J_segment', 'telomerase_RNA', 'vault_RNA'
                    #  'D_segment', 'Y_RNA', 'RNase_MRP_RNA', 'scRNA', 'RNase_P_RNA'
                    if gene.attributes[attribute_filter][0] == attribute_value:
                        chr_id = self.chr_id_from_seq_id(self.annotation_source, gene.chrom)
                        if chr_id != -1:
                            if self.__genes_on_chr[chr_id] is None:
                                self.__genes_on_chr[chr_id] = []
                            self.__genes_on_chr[chr_id].append(gene)

        print("Using " + str(self.annotation_source) + ": genes loaded successfully!")

        if self.annotation_load_type == ANNOTATION_LOAD.GENES:
            return

        # load preferred transcripts and link them to genes in dictionaries
        # preprocess_feature_by_type(__gene_transcripts, feature, 'unconfirmed_transcript')
        # preprocess_feature_by_type(__gene_transcripts, feature, 'V_gene_segment')
        # preprocess_feature_by_type(__gene_transcripts, feature, 'J_gene_segment')
        # preprocess_feature_by_type(__gene_transcripts, feature, 'C_gene_segment')
        features_generator = features_db.features_of_type('mRNA')
        features_mRNA = list(features_generator)
        for mRNA in features_mRNA:
            self.__load_feature_by_type(self.__gene_transcripts, mRNA, 'mRNA')

        print("Using " + str(self.annotation_source) + ": mRNAs loaded successfully!")

        if self.annotation_load_type == ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS:
            return

        # load preferred fragments and link them to transcripts in dictionaries
        features_generator = features_db.features_of_type('CDS')
        features_CDSs = list(features_generator)
        for CDS in features_CDSs:
            self.__load_feature_by_type(self.__transcript_fragments, CDS, 'CDS')

        print("Using " + str(self.annotation_source) + ": CDSs loaded successfully!")

        if self.annotation_load_type == ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_CDS:
            return

        features_generator = features_db.features_of_type('exon')
        features_exons = list(features_generator)
        for exon in features_exons:
            self.__load_feature_by_type(self.__transcript_fragments, exon, 'exon')
        print("Using " + str(self.annotation_source) + ": EXONs loaded successfully!")

        features_generator = features_db.features_of_type('three_prime_UTR')
        features_UTR3s = list(features_generator)
        for utr3 in features_UTR3s:
            self.__load_feature_by_type(self.__transcript_fragments, utr3, 'three_prime_UTR')
        print("Using " + str(self.annotation_source) + ": UTR3s loaded successfully!")

        features_generator = features_db.features_of_type('five_prime_UTR')
        features_UTR5s = list(features_generator)
        for utr5 in features_UTR5s:
            self.__load_feature_by_type(self.__transcript_fragments, utr5, 'five_prime_UTR')
        print("Using " + str(self.annotation_source) + ": UTR5s loaded successfully!")

    def __load_feature_by_type(self, dict_to_fill, feature: Feature, feature_type):
        if feature.featuretype != feature_type: return

        if not feature.attributes.__contains__('Parent'):
            print("Lasha some issue found!")
        if len(feature.attributes['Parent']) != 1:
            print("Lasha some issue found!")

        # parent_gene_id for mRNA and parent_transcript_id for CDS/UTR5'/UTR3'/Exon
        parent_id = feature.attributes['Parent'][0]

        if not dict_to_fill.__contains__(parent_id):
            dict_to_fill[parent_id] = []

        dict_to_fill[parent_id].append(feature)

    def chr_id_from_seq_id(self, annotation_source, id):
        if annotation_source == ANNOTATIONS.ENSEMBL:
            if id == 'MT': return 25
            if id == 'X': return 23
            if id == 'Y': return 24
            try:
                id = int(id)
            except ValueError:
                return -1
            id = id if 1 <= id <= 22 else -1
            return id

        # NCBI 1st column rules: accession number instead of chromosome number
        if not id.__contains__('.'): return -1
        parts = id.split('.')
        if not parts[0][0:3] == 'NC_': return -1
        num = parts[0][3:len(parts[0])]
        x = int(num)

        # let mitochondrion be chrom 25
        if x == 12920: x = 25
        if x < 1 or x > 25:
            return -1
        return x

    # endregion

    # region Find necessary features (sub-features) from a Gene
    def __find_best_gene_transcript(self, chr_id, gene_feature: Feature) -> Feature:
        if not self.__gene_transcripts.__contains__(gene_feature.id):
            # there is no valid mRNA transcript fot this specific gene
            return None

        mRNA_transcripts = self.__gene_transcripts[gene_feature.id]
        best_mRNA_transcript = None

        for transcript in mRNA_transcripts:
            if not transcript.attributes.__contains__('Parent'):
                print("Lasha some issue found!")
            if not len(transcript.attributes) != 1:
                print("Lasha some issue found!")
            if transcript.attributes['Parent'][0] != gene_feature.id:
                print("Lasha some issue found!")
            best_mRNA_transcript = self.__better_annotated_transcript(best_mRNA_transcript, transcript)

        # print("mRNA variants for Gene [" + gene_feature.id + "]: " + str(mRNA_variants))
        # print("well annotated mRNA: " + best_mRNA_transcript.id)

        return best_mRNA_transcript

    def __find_best_gene_fragments(self, chr_id, gene_feature: Feature, force_sorted) -> list[Feature]:
        # mRNA_transcript (start,end) always would fit in parent's (start,end) interval
        mRNA_transcript = self.__find_best_gene_transcript(chr_id, gene_feature)
        if mRNA_transcript is None:
            # it seems there is no valid mRNA for the gene
            # print('no valid mRNA and exons for a gene')
            return []

        fragments = [] if mRNA_transcript.id not in self.__transcript_fragments else self.__transcript_fragments[
            mRNA_transcript.id]

        if not force_sorted:
            return fragments

        fragments.sort(key=lambda x: x.start)

        return fragments

    # endregion

    # region Get sequence of a Gene

    # if sequence for specific chromosome is not loaded then loads it.
    # if its already loaded retrieves it.
    def retrieve_sequence_record(self, chr_id) -> SeqRecord:
        if self.__sequences[chr_id] is None:
            print("Sequence must be loaded during processing! Error...")
        return self.__sequences[chr_id]

    def retrieve_feature_sequence(self, chr_id, feature: Feature) -> str:
        seq_record = self.retrieve_sequence_record(chr_id)
        feature_record = seq_record[feature.start - 1:feature.end]
        if feature.strand == '-':
            feature_record = feature_record.reverse_complement()
        return str(feature_record.seq)

    def retrieve_interval_sequence(self, chr_id, start, end, strand) -> str:
        seq_record = self.retrieve_sequence_record(chr_id)
        interval_record = seq_record[max(0, start - 1):min(end, len(seq_record))]
        if strand == '-': interval_record = interval_record.reverse_complement()
        return str(interval_record.seq)

    # endregion

    # region Additional Methods

    def __better_annotated_transcript(self, A: Feature, B: Feature) -> Feature:
        if A == None: return B
        if B == None: return A

        # no definitive criteria, choose random transcript for a gene to late build 'average' gene model

        return A if random.randint(0, 1000) % 2 == 0 else B

    # between each pairs of transcripts (mRNA) from each gene, finds overlapped fragments (CDS or exons)
    # and returns overlapped segments which has maximum total overlapped length
    def get_fragments_overlap(self, gene_A: Feature, gene_B: Feature, ORF_similarity):
        mRNA_transcripts_A = self.__gene_transcripts[gene_A.id]
        mRNA_transcripts_B = self.__gene_transcripts[gene_B.id]

        max_overlap_segments = []
        max_overlap_length = 0
        for transcript_a in mRNA_transcripts_A:
            fragments_A = [] if transcript_a.id not in self.__transcript_fragments else self.__transcript_fragments[
                transcript_a.id]
            if len(fragments_A) == 0: continue

            for transcript_b in mRNA_transcripts_B:

                fragments_B = [] if transcript_b.id not in self.__transcript_fragments else self.__transcript_fragments[
                    transcript_b.id]
                if len(fragments_B) == 0: continue

                overlap_length = 0
                overlaps = []
                for fragment_a in fragments_A:
                    for fragment_b in fragments_B:
                        if ORF_similarity == 'diff_ORF' and self.are_features_same_framed(fragment_a,
                                                                                          fragment_b): continue
                        if ORF_similarity == 'same_ORF' and not self.are_features_same_framed(fragment_a,
                                                                                              fragment_b): continue

                        if fragment_b.end < fragment_a.start or fragment_b.start > fragment_a.end: continue
                        if fragment_b.end <= fragment_a.end:
                            if fragment_b.start >= fragment_a.start:
                                overlap = (fragment_b.start, fragment_b.end)
                            else:
                                overlap = (fragment_a.start, fragment_b.end)
                        else:
                            if fragment_b.start <= fragment_a.start:
                                overlap = (fragment_a.start, fragment_a.end)
                            else:
                                overlap = (fragment_b.start, fragment_a.end)
                        overlaps.append(overlap)
                        overlap_length += overlap[1] - overlap[0] + 1
                if overlap_length > max_overlap_length:
                    max_overlap_length = overlap_length
                    max_overlap_segments = overlaps

        return max_overlap_segments

    # endregion

    # region Static Methods

    @staticmethod
    def are_genes_overlapped(gene_feature_A: Feature, gene_feature_B: Feature) -> bool:
        l_a, r_a = gene_feature_A.start, gene_feature_A.end
        l_b, r_b = gene_feature_B.start, gene_feature_B.end
        if l_a <= l_b <= r_a: return True
        if l_a <= r_b <= r_a: return True
        if l_b <= l_a <= r_b: return True
        return False

    @staticmethod
    def get_genes_overlap(gene_feature_A: Feature, gene_feature_B: Feature):
        l_a, r_a = gene_feature_A.start, gene_feature_A.end
        l_b, r_b = gene_feature_B.start, gene_feature_B.end
        if l_a <= l_b <= r_a: return l_b, min(r_a, r_b)
        if l_a <= r_b <= r_a: return max(l_a, l_b), r_b
        if l_b <= l_a <= r_b: return l_a, l_b
        assert False

    @staticmethod
    def are_features_same_framed(fragment_a: Feature, fragment_b: Feature) -> bool:
        l = fragment_a.start + int(fragment_a.frame) if fragment_a.strand == '+' else fragment_a.end - int(
            fragment_a.frame)
        r = fragment_b.start + int(fragment_b.frame) if fragment_b.strand == '+' else fragment_b.end - int(
            fragment_b.frame)
        return l % 3 == r % 3

    # retrieves nucleotide composition of sequence
    # output: (C_count,G_count,A_count,T_count)
    @staticmethod
    def sequence_composition(sequence):
        stats = [0, 0, 0, 0]
        for char in sequence:
            if char == 'C': stats[0] += 1
            if char == 'G': stats[1] += 1
            if char == 'A': stats[2] += 1
            if char == 'T': stats[3] += 1
        # if stats[0] + stats[1] + stats[2] + stats[3] != len(sequence):
        #     print("Fuck you! unknown symbol in fragment")

        return stats

    @staticmethod
    def sequence_composition_by_parts(sequence, k):
        part_length = len(sequence) // k

        compositions_by_parts = [None] * k
        for i in range(0, k):
            l = i * part_length
            r = l + part_length - 1
            compositions_by_parts[i] = GenomeWorker.sequence_composition(sequence[l:(r + 1)])

        return compositions_by_parts

    # endregion

    # region Analyze Gene nucleotide composition(stranded) by fragments

    # gene structure: (UTR'5-CDS1)=EXON1 - INTRON - CDS2=EXON2 - ... - CDS(n-1)=EXON(n-1) - INTRON - (CDSn-UTR'3)=EXONn
    # subregions of interest: UTR'5_Procession, UTR'5, inner CDSs, inner Introns, UTR'3, UTR'3_Procession
    # for each subregion calculate (C_count,G_count,A_count,T_count)

    def __get_regional_merged_sequences_from_gene(self, chr_id, gene, procession_length):
        # for NCBI database
        fragments = self.__find_best_gene_fragments(chr_id, gene, force_sorted=True)

        utr5_sequence = ""
        utr3_sequence = ""
        cds_sequence = ""
        intron_sequence = ""
        utr5_procession_sequence = ""
        utr3_procession_sequence = ""

        if len(fragments) == 0: return [utr5_procession_sequence, utr5_sequence, cds_sequence, intron_sequence,
                                        utr3_sequence,
                                        utr3_procession_sequence]

        for fragment in fragments:
            if fragment.featuretype == 'five_prime_UTR':
                if gene.strand == '+':
                    utr5_sequence += self.retrieve_feature_sequence(chr_id, fragment)
                else:
                    utr5_sequence = self.retrieve_feature_sequence(chr_id, fragment) + utr5_sequence
            if fragment.featuretype == 'three_prime_UTR':
                if gene.strand == '+':
                    utr3_sequence += self.retrieve_feature_sequence(chr_id, fragment)
                else:
                    utr3_sequence = self.retrieve_feature_sequence(chr_id, fragment) + utr3_sequence

            if fragment.featuretype == 'CDS':
                if gene.strand == '+':
                    cds_sequence += self.retrieve_feature_sequence(chr_id, fragment)
                else:
                    cds_sequence = self.retrieve_feature_sequence(chr_id, fragment) + cds_sequence

        last_exon = None
        introns_len = 0

        left_fragment = None
        right_fragment = None

        for fragment in fragments:
            if fragment.featuretype != 'exon': continue

            if left_fragment is None or left_fragment.start > fragment.start:
                left_fragment = fragment
            if right_fragment is None or right_fragment.end < fragment.end:
                right_fragment = fragment

            if last_exon is None:
                last_exon = fragment
                continue
            if gene.strand == '+':
                intron_sequence += self.retrieve_interval_sequence(chr_id, last_exon.end + 1, fragment.start - 1,
                                                                   fragment.strand)
            else:
                intron_sequence = self.retrieve_interval_sequence(chr_id, last_exon.end + 1, fragment.start - 1,
                                                                  fragment.strand) + intron_sequence
            introns_len += fragment.start - last_exon.end - 1
            last_exon = fragment

        if left_fragment.strand == '+':
            utr5_procession_sequence = self.retrieve_interval_sequence(chr_id, left_fragment.start - procession_length,
                                                                       left_fragment.start - 1, '+')

            utr3_procession_sequence = self.retrieve_interval_sequence(chr_id, right_fragment.end + 1,
                                                                       right_fragment.end + procession_length, '+')

        else:
            utr3_procession_sequence = self.retrieve_interval_sequence(chr_id, left_fragment.start - procession_length,
                                                                       left_fragment.start - 1, '-')

            utr5_procession_sequence = self.retrieve_interval_sequence(chr_id, right_fragment.end + 1,
                                                                       right_fragment.end + procession_length, '-')

        return [utr5_procession_sequence, utr5_sequence, cds_sequence, intron_sequence, utr3_sequence,
                utr3_procession_sequence]

    # output: occurrences[part][base] = [k][4]
    def analyze_gene_occurrences(self, chr_id, gene: Feature, procession_length):
        regional_sequences = self.__get_regional_merged_sequences_from_gene(chr_id, gene, procession_length)

        utr5_procession_sequence = regional_sequences[0]
        utr5_sequence = regional_sequences[1]
        cds_sequence = regional_sequences[2]
        intron_sequence = regional_sequences[3]
        utr3_sequence = regional_sequences[4]
        utr3_procession_sequence = regional_sequences[5]

        occurrences_UTR5_procession = GenomeWorker.sequence_composition(utr5_procession_sequence)
        occurrences_UTR5 = GenomeWorker.sequence_composition(utr5_sequence)
        occurrences_CDS = GenomeWorker.sequence_composition(cds_sequence)
        occurrences_Introns = GenomeWorker.sequence_composition(intron_sequence)
        occurrences_UTR3 = GenomeWorker.sequence_composition(utr3_sequence)
        occurrences_UTR3_procession = GenomeWorker.sequence_composition(utr3_procession_sequence)

        return [occurrences_UTR5_procession, occurrences_UTR5, occurrences_CDS, occurrences_Introns, occurrences_UTR3,
                occurrences_UTR3_procession]

    # output: occurrences[region][part][base] = [6][k][4]
    def analyze_gene_occurrences_by_parts(self, chr_id, gene: Feature, k, procession_length):
        regional_sequences = self.__get_regional_merged_sequences_from_gene(chr_id, gene, procession_length)

        utr5_procession_sequence = regional_sequences[0]
        utr5_sequence = regional_sequences[1]
        cds_sequence = regional_sequences[2]
        intron_sequence = regional_sequences[3]
        utr3_sequence = regional_sequences[4]
        utr3_procession_sequence = regional_sequences[5]

        # subregions of interest: UTR'5_Procession, UTR'5, inner CDSs, inner Introns, UTR'3, UTR'3_Procession
        # each region divided into k-part
        # for each part (a,c,g,t) is calculated
        # occurrences[region][part][base] = [6][k][4]
        occurrences = [self.sequence_composition_by_parts(utr5_procession_sequence, k),
                       self.sequence_composition_by_parts(utr5_sequence, k),
                       self.sequence_composition_by_parts(cds_sequence, k),
                       self.sequence_composition_by_parts(intron_sequence, k),
                       self.sequence_composition_by_parts(utr3_sequence, k),
                       self.sequence_composition_by_parts(utr3_procession_sequence, k)]
        return occurrences

    # endregion
