import gc
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

    def clear_precomputed_values(self):
        # only transcripts for a gene are saved from previous computations later,
        # criteria could be changed, so this method make sures new computing by new criteria
        self.__gene_transcript_by_criteria.clear()

    def release_memory(self):
        keys3 = []
        for key in self.__gene_transcript_by_criteria:
            keys3.append(key)
        for key in keys3:
            del self.__gene_transcript_by_criteria[key]

        keys = []
        for key in self.__transcript_fragments:
            keys.append(key)
        for key in keys:
            del self.__transcript_fragments[key]

        keys2 = []
        for key in self.__gene_transcripts:
            keys2.append(key)
        for key in keys2:
            del self.__gene_transcripts[key]

        del keys
        del keys2
        del keys3
        del self.__genes_on_chr
        del self.__sequences
        gc.collect()

    # region Variables and Retrievers

    imported_protein_coding_genes = 0
    ignored_protein_coding_genes = 0

    # gene features clustered by chromosome, they are on.
    # Example: genes_on_chr[22] = ['gene' features on chromosome 22]
    __genes_on_chr = []

    # {chr_id,DNA sequence of chromosome with id chr_id}
    __sequences = []

    #   {gene_id, [list of mRNA features whose parent is gene_id]}
    __gene_transcripts = {}

    #   {transcript_id, [list of 'CDS' and/or 'UTR(3'/5')' and/or 'exons' features whose parent is mRNA_transcript_id]}
    __transcript_fragments = {}

    #   {gene_id, [list of 'best' mRNA features whose parent is gene_id]}
    __gene_transcript_by_criteria = {}

    __gene_names_set = {}
    __gene_accessions_set = {}

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

    # region Load methods for Features (from annotation) and Sequences (from assembly)

    def __load_requested_sequences(self):
        if self.sequence_load_type == SEQUENCE_LOAD.LOAD:
            sequence_file_path = self.__get_sequence_file_path()
            for record in SeqIO.parse(sequence_file_path, 'fasta'):
                chr_id = self.chr_id_from_seq_id(self.annotation_source, record.id)
                if chr_id == -1: continue
                self.__sequences[chr_id] = record
            loaded_chromosomes = 0
            for chr_id in range(1, NUMBER_OF_CHROMOSOMES[self.species.value] + 1):
                loaded_chromosomes += 1 if self.__sequences[chr_id] is not None and len(
                    self.__sequences[chr_id]) > 0 else 0
            print(
                "    Using " + self.annotation_source.name + ": sequence data loaded from chromosomes " + str(
                    loaded_chromosomes) + "/" +
                str(NUMBER_OF_CHROMOSOMES[self.species.value]))

    # if database does not exists for specific chromosome, then
    # builds database from annotation file and fills/loads all
    # necessary features from the database.
    #
    # if database exists, then fills/loads all necessary data structures

    def __check_gene_for_filters(self, gene: Feature):
        # ignore genes without Name or description
        if not gene.attributes.__contains__('Name') or not gene.attributes.__contains__('description'):
            return False

        # record must contain only 1 Name
        assert len(gene.attributes['Name']) == 1

        # ignore 'readthrough' genes
        if gene.attributes['description'][0].__contains__('readthrough'):
            return False

        # ignore genes with duplicate Names
        gene_name = GenomeWorker.get_gene_name(gene)
        if self.__gene_names_set.__contains__(gene_name):
            return False

        # ignore genes if gene with same NCBI accession already imported
        gene_accession = GenomeWorker.get_gene_accession(gene)
        if gene_accession != 'no_acc' and self.__gene_accessions_set.__contains__(gene_accession):
            return False

        return True

    def __load_requested_features(self):
        if not exists(self.__get_annotation_db_path()):
            print("Creating database for " + self.annotation_source.name + " only for first RUN it takes that long!")
            gffutils.create_db(self.__get_annotation_file_path(),
                               dbfn=self.__get_annotation_db_path(),
                               verbose=True, force=False, keep_order=False, merge_strategy='create_unique',
                               sort_attribute_values=False)

        features_db = gffutils.FeatureDB(self.__get_annotation_db_path(), keep_order=False)

        features_generator = features_db.features_of_type('gene')
        feature_genes = list(features_generator)

        # choose genes who has protein_coding attribute and additional filter values
        for gene in feature_genes:
            attribute_filter = 'gene_biotype' if self.annotation_source == ANNOTATIONS.NCBI else 'biotype'
            attribute_value = 'protein_coding'

            if gene.attributes.__contains__(attribute_filter):
                assert len(gene.attributes[attribute_filter]) == 1
                if gene.attributes[attribute_filter][0] == attribute_value:
                    chr_id = self.chr_id_from_seq_id(self.annotation_source, gene.chrom)
                    if chr_id != -1:
                        if self.__check_gene_for_filters(gene):
                            if self.__genes_on_chr[chr_id] is None:
                                self.__genes_on_chr[chr_id] = []

                            self.imported_protein_coding_genes += 1
                            self.__genes_on_chr[chr_id].append(gene)

                            # store these values, for filtering upcoming genes
                            self.__gene_names_set[GenomeWorker.get_gene_name(gene)] = len(self.__genes_on_chr[chr_id])
                            self.__gene_accessions_set[GenomeWorker.get_gene_accession(gene)] = len(
                                self.__genes_on_chr[chr_id])
                        else:
                            self.ignored_protein_coding_genes += 1

        loaded_chromosomes = 0
        for chr_id in range(1, NUMBER_OF_CHROMOSOMES[self.species.value] + 1):
            if self.__genes_on_chr[chr_id] is not None and len(self.__genes_on_chr[chr_id]) > 0:
                loaded_chromosomes += 1
            else:
                print(chr_id)
        print("    Using " + self.annotation_source.name + ": " + str(
            self.imported_protein_coding_genes) + " genes loaded (" + str(
            self.ignored_protein_coding_genes) + " filtered out) from chromosomes " + str(
            loaded_chromosomes) + "/" + str(NUMBER_OF_CHROMOSOMES[self.species.value]))

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

        print("    Using " + self.annotation_source.name + ": mRNAs loaded successfully!")

        if self.annotation_load_type == ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS:
            return

        # load preferred fragments and link them to transcripts in dictionaries
        features_generator = features_db.features_of_type('CDS')
        features_CDSs = list(features_generator)
        for CDS in features_CDSs:
            self.__load_feature_by_type(self.__transcript_fragments, CDS, 'CDS')

        print("    Using " + self.annotation_source.name + ": CDSs loaded successfully!")

        if self.annotation_load_type == ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_CDS:
            return

        features_generator = features_db.features_of_type('exon')
        features_exons = list(features_generator)
        for exon in features_exons:
            self.__load_feature_by_type(self.__transcript_fragments, exon, 'exon')
        print("    Using " + self.annotation_source.name + ": EXONs loaded successfully!")

        features_generator = features_db.features_of_type('three_prime_UTR')
        features_UTR3s = list(features_generator)
        for utr3 in features_UTR3s:
            self.__load_feature_by_type(self.__transcript_fragments, utr3, 'three_prime_UTR')
        print("    Using " + self.annotation_source.name + ": UTR3s loaded successfully!")

        features_generator = features_db.features_of_type('five_prime_UTR')
        features_UTR5s = list(features_generator)
        for utr5 in features_UTR5s:
            self.__load_feature_by_type(self.__transcript_fragments, utr5, 'five_prime_UTR')
        print("    Using " + self.annotation_source.name + ": UTR5s loaded successfully!")

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
            if self.species == SPECIES.Homo_sapiens or self.species == SPECIES.Mus_musculus \
                    or self.species == SPECIES.Rattus_norvegicus or self.species == SPECIES.Danio_rerio:

                if id == 'MT': return -1  # ignore mito
                if id == 'Y': return NUMBER_OF_CHROMOSOMES[self.species.value]
                if id == 'X': return NUMBER_OF_CHROMOSOMES[self.species.value] - 1

                try:
                    index = int(id)
                except ValueError:
                    index = -1

                return index if 1 <= index <= NUMBER_OF_CHROMOSOMES[self.species.value] else -1
            else:
                if self.species == SPECIES.Caenorhabditis_elegans:
                    if not ENSEMBL_CHR_MAP_FOR_Caenorhabditis.__contains__(id): return -1
                    return ENSEMBL_CHR_MAP_FOR_Caenorhabditis[id]
                elif self.species == SPECIES.Drosophila_melanogaster:
                    if not ENSEMBL_CHR_MAP_FOR_DROSOPHILA.__contains__(id): return -1
                    return ENSEMBL_CHR_MAP_FOR_DROSOPHILA[id]
                else:
                    return -1
        # NCBI 1st column rules: accession number instead of chromosome number
        if not id.__contains__('.'): return -1
        parts = id.split('.')

        if self.species == SPECIES.Drosophila_melanogaster:
            if parts[0] == 'NC_004354': return 1
            if parts[0] == 'NT_033779': return 2
            if parts[0] == 'NT_033778': return 3
            if parts[0] == 'NT_037436': return 4
            if parts[0] == 'NT_033777': return 5
            if parts[0] == 'NC_004353': return 6
            # ignore mito
            if parts[0] == 'NC_024511': return -1
            return -1

        if not parts[0][0:3] == 'NC_': return -1
        num = parts[0][3:len(parts[0])]
        x = int(num)

        if self.species == SPECIES.Homo_sapiens:
            if x == 12920: return -1  # ignore mito
            if x < 1 or x > 24: return -1
            return x
        elif self.species == SPECIES.Rattus_norvegicus:
            if x == 1665: return -1  # ignore mito
            base = 51335
            x = x - base
            if x < 1 or x > 22: return -1
            return x
        elif self.species == SPECIES.Mus_musculus:
            if x == 5089: return -1  # ignore mito
            base = 66
            x = x - base
            if x < 1 or x > 21: return -1
            return x
        elif self.species == SPECIES.Caenorhabditis_elegans:
            if x == 1328: return -1  # ignore mito
            base = 3278
            x = x - base
            if x < 1 or x > 6: return -1
            return x
        elif self.species == SPECIES.Danio_rerio:
            if x == 2333: return -1  # ignore mito
            base = 7111
            x = x - base
            if x < 1 or x > 25: return -1
            return x
        else:
            print("NCBI chrid needs conversation to index!!!!!!!!!!!!!")
            return -1

    # endregion

    def __get_transcript_score_by_criteria(self, mRNA_transcript: Feature, criteria: TRANSCRIPT_CRITERIA):
        if criteria == TRANSCRIPT_CRITERIA.LONGEST: return mRNA_transcript.end - mRNA_transcript.start + 1

        fragments = self.__transcript_fragments[mRNA_transcript.id]

        score = 0
        for fragment in fragments:
            # for other criteria CDS must be counted: LONGEST_CDS and LONGEST_CDS_AND_UTRs
            if fragment.featuretype == 'CDS':
                score += fragment.end - fragment.start + 1

            # for LONGEST_CDS_AND_UTRs we must also count UTR sequences
            # same as exon
            if criteria == TRANSCRIPT_CRITERIA.LONGEST_CDS_AND_UTRs:
                if fragment.featuretype == 'three_prime_UTR' or fragment.featuretype == 'five_prime_UTR':
                    score += fragment.end - fragment.start + 1

        return score

    # region Find necessary feature (transcript/mRNA) from a Gene
    def find_gene_transcript_by_criteria(self, gene_feature: Feature, criteria) -> Feature:
        if self.__gene_transcript_by_criteria.__contains__(gene_feature.id):
            # if 'best' gene transcript for a gene already found
            return self.__gene_transcript_by_criteria[gene_feature.id]

        if not self.__gene_transcripts.__contains__(gene_feature.id):
            # there is no valid mRNA transcript fot this specific gene
            return None

        mRNA_transcripts = self.__gene_transcripts[gene_feature.id]
        assert len(mRNA_transcripts) != 0

        # if criteria is RANDOM, just choose if randomly from the lists of all transcripts
        if criteria == TRANSCRIPT_CRITERIA.RANDOM:
            return mRNA_transcripts[random.randint(0, len(mRNA_transcripts))]

        best_mRNA_transcript = None
        best_score = 0

        for transcript in mRNA_transcripts:
            assert transcript.attributes.__contains__('Parent')
            assert len(transcript.attributes) != 1
            assert transcript.attributes['Parent'][0] == gene_feature.id

            score = self.__get_transcript_score_by_criteria(transcript, criteria)

            if score > best_score:
                best_mRNA_transcript = transcript
                best_score = score

        # save in dict, so accessing later would not require additional findings/calculations
        self.__gene_transcript_by_criteria[gene_feature.id] = best_mRNA_transcript

        return best_mRNA_transcript

    def __find_gene_fragments_by_criteria(self, chr_id, gene_feature: Feature, criteria, force_sorted) -> list[Feature]:
        # mRNA_transcript (start,end) always would fit in parent's (start,end) interval
        mRNA_transcript = self.find_gene_transcript_by_criteria(gene_feature, criteria)
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
    def get_fragments_overlap(self, gene_A: Feature, gene_B: Feature, ORF_similarity, criteria: TRANSCRIPT_CRITERIA):
        mRNA_transcripts_a = self.find_gene_transcript_by_criteria(gene_A, criteria)
        mRNA_transcripts_b = self.find_gene_transcript_by_criteria(gene_B, criteria)

        max_overlap_segments = []
        max_overlap_length = 0

        if mRNA_transcripts_a is None or mRNA_transcripts_b is None: return max_overlap_segments
        fragments_A = [] if mRNA_transcripts_a.id not in self.__transcript_fragments else self.__transcript_fragments[
            mRNA_transcripts_a.id]
        fragments_B = [] if mRNA_transcripts_b.id not in self.__transcript_fragments else self.__transcript_fragments[
            mRNA_transcripts_b.id]
        if len(fragments_B) == 0 or len(fragments_B) == 0: return max_overlap_segments

        overlap_length = 0
        overlaps = []
        for fragment_a in fragments_A:
            if fragment_a.featuretype != 'CDS': continue
            for fragment_b in fragments_B:
                if fragment_b.featuretype != 'CDS': continue

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
    def get_gene_accession(gene: Feature):
        desc_parts = gene.attributes['description'][0].split("Acc:")
        if len(desc_parts) != 2: return "no_acc"  # it seems record does not have NCBI ID
        ncbi_id = desc_parts[1][0:(len(desc_parts[1]) - 1)]
        return ncbi_id

    @staticmethod
    def get_gene_name(gene: Feature):
        assert gene.attributes.__contains__('Name')  # loaded genes must be filtered
        assert len(gene.attributes['Name']) == 1
        return gene.attributes['Name'][0]

    @staticmethod
    def are_features_overlapped(gene_feature_A: Feature, gene_feature_B: Feature) -> bool:
        # if gene_feature_A is None or gene_feature_B is None: return False
        l_a = gene_feature_A.start
        r_a = gene_feature_A.end
        l_b = gene_feature_B.start
        r_b = gene_feature_B.end
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
        # sliding window with size k.

        if len(sequence) % k != 0:  # place equally distanced 'N's to make it k-dividable for later smooth analyze
            total_gaps_to_place = k - (len(sequence) % k)
            gaps_placed = 0
            new_seq = ""
            length = len(sequence)
            for index in range(length):
                new_seq += sequence[index]
                # index / len * gaps_to_place = gaps_placed   index*(gaps_to_place+1) // len= 15
                if round((index + 1) * total_gaps_to_place / length) > gaps_placed:
                    n = round((index + 1) * total_gaps_to_place / length) - gaps_placed
                    new_seq += 'N' * n  # unknown
                    gaps_placed += n
            sequence = new_seq

        part_length = len(sequence) // k  # it is now 100% dividable by k
        compositions_by_parts = []

        for i in range(0, k):
            l = i * part_length
            r = l + part_length - 1
            compositions_by_parts.append(GenomeWorker.sequence_composition(sequence[l:(r + 1)]))

        return compositions_by_parts

    # endregion

    # region Analyze Gene nucleotide composition(stranded) by fragments

    # gene structure: (UTR'5-CDS1)=EXON1 - INTRON - CDS2=EXON2 - ... - CDS(n-1)=EXON(n-1) - INTRON - (CDSn-UTR'3)=EXONn
    # subregions of interest: UTR'5_Procession, UTR'5, inner CDSs, inner Introns, UTR'3, UTR'3_Procession
    # for each subregion calculate (C_count,G_count,A_count,T_count)

    def __get_regional_merged_sequences_from_gene(self, chr_id, gene, procession_length, criteria):
        # for NCBI database
        fragments = self.__find_gene_fragments_by_criteria(chr_id, gene, criteria, force_sorted=True)

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
    def analyze_gene_occurrences(self, chr_id, gene: Feature, procession_length, criteria: TRANSCRIPT_CRITERIA):
        regional_sequences = self.__get_regional_merged_sequences_from_gene(chr_id, gene, procession_length, criteria)

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
    def analyze_gene_occurrences_by_parts(self, chr_id, gene: Feature, k, procession_length,
                                          criteria: TRANSCRIPT_CRITERIA):
        regional_sequences = self.__get_regional_merged_sequences_from_gene(chr_id, gene, procession_length, criteria)

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
