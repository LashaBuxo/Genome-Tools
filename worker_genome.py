from os.path import exists
import gffutils
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from gffutils import Feature
import random

# Specific tools for working with Ensembl annotation and with sequence data
from worker_genome_values import *
from worker_genome_enums import *
from worker_analyzer import *


class GenomeWorker:
    def __init__(self, species: SPECIES, annotation_load_type: ANNOTATION_LOAD,
                 sequence_load_type: SEQUENCE_LOAD, load_APPRIS=False, load_MitoCARTA=False):

        print(f"loading data for {species}:")

        self.species = species
        self.annotation_load_type = annotation_load_type
        self.sequence_load_type = sequence_load_type
        self.load_APPRIS = load_APPRIS
        self.load_MitoCARTA = load_MitoCARTA

        self.imported_protein_coding_genes = 0
        self.ignored_protein_coding_genes = 0
        self.ignored_genes_by_types = {}

        # storage for gene or transcript features to access by id
        self.__loaded_feature_by_id = {}

        # gene features clustered by chromosome, they are on.
        self.chromosome_name2index = {}
        self.__genes_on_chr = [None]
        self.__gene_symbols_set = {}
        self.__gene_accessions_set = {}

        self.__gene_transcripts = {}  # {gene_id, [list of mRNA features/isoforms whose parent is gene_id]}

        self.__transcript_fragments = {}
        self.__transcript_fragments_is_sorted = {}  # optimizes speed, if we are not sorting it if we don't need

        self.__gene_transcript_by_criteria = {}  # {gene_id, representative mRNA features whose parent is gene_id}

        # {chr_id,DNA sequence of chromosome with id chr_id}
        self.__sequences = [None]

        self.__load_requested_features()
        self.__load_requested_sequences()

    # region easy retriever methods

    def get_chromosome_length(self, chr_id):
        return len(self.retrieve_sequence_record(chr_id))

    def get_feature_chromosomal_position(self, feature_id):
        feature = self.feature_by_id(feature_id)
        chr_id = self.seq_id2_chr_index(feature.chrom)

        if self.sequence_load_type == SEQUENCE_LOAD.NOT_LOAD:
            return chr_id, ''
        chr_len = len(self.retrieve_sequence_record(chr_id))
        return chr_id, f"{chr_id} ({int(feature.start / chr_len * 100)}%)"

    def chromosomes_count(self):
        return len(self.chromosome_name2index)

    def chromosome_length(self, chr_id):
        return len(self.__sequences[chr_id])

    def genes_count_on_chr(self, chr_id):
        assert chr_id <= self.chromosomes_count()
        return 0 if self.__genes_on_chr[chr_id] is None else len(self.__genes_on_chr[chr_id])

    def gene_by_indexes(self, chr_id, index) -> Feature:
        assert chr_id <= self.chromosomes_count()
        return self.__genes_on_chr[chr_id][index]

    def gene_by_symbol(self, symbol) -> Feature:
        if self.__gene_symbols_set.__contains__(symbol):
            return self.feature_by_id(self.__gene_symbols_set[symbol])
        return None

    def get_transcript_parent(self, transcript_id):
        feature = self.feature_by_id(transcript_id)
        return self.feature_by_id(feature.attributes['Parent'][0])

    def feature_by_id(self, feature_id) -> Feature:
        return self.__loaded_feature_by_id[feature_id] if self.__loaded_feature_by_id.__contains__(feature_id) else None

    def __get_annotation_file_path(self):
        return ENSEMBL_ANNOTATIONS[self.species.value]

    def __get_annotation_db_path(self):
        added_name = self.__get_annotation_file_path().replace('/', '_')
        annotation_name = "ENSEMBL"
        generated_path = GENOME_DATABASES_DIRECTORY + self.species.name + "_" + annotation_name + "_" + added_name
        return generated_path + '.db'

    def __get_sequence_file_path(self):
        return ENSEMBL_SEQUENCES[self.species.value]

    def get_transcripts_from_gene(self, gene_id) -> list[Feature]:
        return self.__gene_transcripts[gene_id] if self.__gene_transcripts.__contains__(gene_id) else []

    def get_gene_incomplete_level(self, gene_id):
        transcripts = self.get_transcripts_from_gene(gene_id)
        complete_cnt = 0
        protein_coding_isoforms = 0
        for transcript in transcripts:
            frags = self.get_fragments_from_transcript(transcript.id)
            utr5_found = False
            utr3_found = False
            CDS_found = False

            for frag in frags:
                if frag.featuretype == 'three_prime_UTR':
                    utr3_found = True
                if frag.featuretype == 'CDS':
                    CDS_found = True
                if frag.featuretype == 'five_prime_UTR':
                    utr5_found = True
            if CDS_found:
                protein_coding_isoforms += 1
                complete_cnt += 1 if utr3_found and utr5_found else 0
        return (protein_coding_isoforms - complete_cnt) / protein_coding_isoforms

    def get_transcript_CDS_length(self, transcript_id):
        frags = self.__transcript_fragments[transcript_id]
        CDS_length = 0
        for frag in frags:
            if frag.featuretype == 'CDS':
                CDS_length += frag.end - frag.start + 1
        return CDS_length

    def get_transcript_first_CDS_len(self, transcript_id):
        frags = self.get_fragments_from_transcript(transcript_id)
        feature = self.feature_by_id(transcript_id)

        for index in range(0, len(frags)):
            frag = frags[index] if feature.strand == '+' else frags[len(frags) - 1 - index]
            if frag.featuretype == 'CDS': return frag.end - frag.start + 1

        assert False

    def get_fragments_from_transcript(self, transcript_id) -> list[Feature]:
        if transcript_id not in self.__transcript_fragments: return []

        fragments = self.__transcript_fragments[transcript_id]
        if self.__transcript_fragments_is_sorted.__contains__(transcript_id):
            return fragments

        fragments.sort(key=lambda x: x.start)

        self.__transcript_fragments_is_sorted[transcript_id] = True
        self.__transcript_fragments[transcript_id] = fragments

        return fragments

    # endregion

    # region load methods for Annotation, Assembly(sequence)

    def __load_requested_sequences(self):
        if self.sequence_load_type != SEQUENCE_LOAD.LOAD:
            return

        sequence_file_path = self.__get_sequence_file_path()
        for record in SeqIO.parse(sequence_file_path, 'fasta'):
            chr_id = self.seq_id2_chr_index(record.id)
            if chr_id == -1: continue
            self.__sequences[chr_id] = record

        loaded_chromosomes = 0
        for _, chr_index in self.chromosome_name2index.items():
            loaded_chromosomes += 1 if self.__sequences[chr_index] is not None and len(
                self.__sequences[chr_index]) > 0 else 0

        print(f"\tDNA chromosomes loaded {loaded_chromosomes}/{self.chromosomes_count()}!")

    def __check_gene_for_Ensembl_filters(self, gene: Feature, feature_ids, feature_syms, feature_accs):
        # https: // www.biostars.org / p / 5304 /  # 9521037
        # ignore genes with miscellaneous chromosome/scaffold names
        chr_id = self.seq_id2_chr_index(gene.chrom)
        if chr_id == -1: return False, "located_on_miscellaneous_scaffold"

        gene_symbol = GenomeWorker.get_gene_symbol(gene)
        gene_accession = GenomeWorker.get_gene_accession(gene)
        description = GenomeWorker.get_gene_description(gene)

        # ignore genes without description or gene_symbol
        # if description == "no_desc": return False, "no_desc"
        # if gene_symbol == "no_sym": return False, "no_symbol"
        # if gene_accession == "no_acc": return False, "no_accession"

        # ignore genes if they are pseudogene, novel or predicted, readthrough
        if description.__contains__('readthrough'):
            return False, "is_readthrough"
        if description.__contains__('pseudogene') or description.__contains__('pseudogene'):
            return False, "is_pseudogene"
        if description.__contains__('novel') or description.__contains__('Novel'):
            return False, "is_novel"
        if description.__contains__('predicted') or description.__contains__('Predicted'):
            return False, "is_predicted"

        # ignore genes with duplicate Names or accessions
        if gene_symbol != "no_sym" and feature_syms.__contains__(gene_symbol):
            return False, "name_duplicated"
        # ignore genes if gene with same NCBI accession already imported
        if gene_accession != 'no_acc' and feature_accs.__contains__(gene_accession):
            return False, "acc_duplicated"

        # for any annotation, there must not be 2 gene record with same id
        assert not feature_ids.__contains__(gene.id)

        return True, "passed_filter"

    def filter_by_ensembl_attributes(self):
        feature_ids = {}
        feature_syms = {}
        feature_accs = {}
        for chr_index in range(1, self.chromosomes_count() + 1):
            genes_on_chr = self.__genes_on_chr[chr_index]
            new_genes_on_chr = []
            cnt = len(genes_on_chr)
            for i in range(0, cnt):
                gene = genes_on_chr[i]
                is_valid, status = self.__check_gene_for_Ensembl_filters(gene, feature_ids, feature_syms, feature_accs)
                if is_valid:
                    new_genes_on_chr.append(gene)
                    feature_ids[gene.id] = True
                    feature_syms[GenomeWorker.get_gene_symbol(gene)] = True
                    feature_accs[GenomeWorker.get_gene_accession(gene)] = True
                else:
                    self.filter_gene_by_reason(status)

            self.__genes_on_chr[chr_index] = new_genes_on_chr

    def filter_gene_by_reason(self, status):
        if not self.ignored_genes_by_types.__contains__(status):
            self.ignored_genes_by_types[status] = 0
        self.ignored_genes_by_types[status] += 1
        self.ignored_protein_coding_genes += 1
        self.imported_protein_coding_genes -= 1

    def filter_nested_genes(self):
        self.ignored_genes_by_types["nested"] = 0
        for chr_index in range(1, self.chromosomes_count() + 1):
            genes_on_chr = self.__genes_on_chr[chr_index]
            new_genes_on_chr = []
            cnt = len(genes_on_chr)
            for i in range(0, cnt):
                gene1 = genes_on_chr[i]
                flag = False
                for j in range(0, cnt):
                    if i == j: continue
                    gene2 = genes_on_chr[j]
                    if self._is_gene_semi_nested_on_same_strand(gene1, gene2):
                        flag = True
                        break

                if not flag:
                    new_genes_on_chr.append(gene1)
                else:
                    self.filter_gene_by_reason("nested")

            self.__genes_on_chr[chr_index] = new_genes_on_chr

    threshold = 0.9

    def _is_gene_semi_nested_on_same_strand(self, nested_gene: Feature, parent_gene: Feature):
        if nested_gene.strand != parent_gene.strand: return False
        len = self.get_overlap_length(nested_gene, parent_gene)
        if len > (nested_gene.end - nested_gene.start + 1) * self.threshold: return True
        return False

    # if database does not exists for specific chromosome, then
    # builds database from annotation file and fills/loads all
    # necessary features from the database.
    #
    # if database exists, then fills/loads all necessary data structures
    def __load_requested_features(self):
        if not exists(self.__get_annotation_db_path()):
            print("Creating database for " + self.species.name + " only for first RUN it takes that long!")
            gffutils.create_db(self.__get_annotation_file_path(),
                               dbfn=self.__get_annotation_db_path(),
                               verbose=True, force=False, keep_order=False, merge_strategy='create_unique',
                               sort_attribute_values=False)

        features_db = gffutils.FeatureDB(self.__get_annotation_db_path(), keep_order=False)
        features_generator = features_db.features_of_type('gene')
        feature_genes = list(features_generator)

        file = open(self.__get_annotation_file_path(), 'r')

        lines = file.readlines()
        chr_index = 0
        for line in lines:
            if line.__contains__("#!genome-build"): break
            if line.__contains__("##sequence-region"):
                arr = line.split(' ')
                chr_index += 1
                self.chromosome_name2index[arr[3]] = chr_index

        self.__genes_on_chr = [None] * (self.chromosomes_count() + 1)
        self.__sequences = [None] * (self.chromosomes_count() + 1)

        # choose genes who has protein_coding attribute and additional filter values
        for gene in feature_genes:
            attribute_filter = 'biotype'
            attribute_value = 'protein_coding'
            if gene.attributes.__contains__(attribute_filter):
                assert len(gene.attributes[attribute_filter]) == 1
                if gene.attributes[attribute_filter][0] == attribute_value:
                    chr_id = self.seq_id2_chr_index(gene.chrom)
                    if self.__genes_on_chr[chr_id] is None:
                        self.__genes_on_chr[chr_id] = []

                    self.imported_protein_coding_genes += 1
                    self.__genes_on_chr[chr_id].append(gene)  # store gene feature on chromosome list

        self.filter_nested_genes()
        self.filter_by_ensembl_attributes()

        # store gene feature->id, id->symbol, id->accession
        for chr_index in range(1, self.chromosomes_count() + 1):
            genes_on_chr = self.__genes_on_chr[chr_index]
            cnt = len(genes_on_chr)
            for i in range(0, cnt):
                gene = genes_on_chr[i]
                self.__loaded_feature_by_id[gene.id] = gene
                self.__gene_symbols_set[GenomeWorker.get_gene_symbol(gene)] = gene.id
                self.__gene_accessions_set[GenomeWorker.get_gene_accession(gene)] = gene.id

        loaded_chromosomes = 0
        for _, chr_index in self.chromosome_name2index.items():
            if self.__genes_on_chr[chr_index] is not None and len(self.__genes_on_chr[chr_index]) > 0:
                loaded_chromosomes += 1
            else:
                print(f'WARNING: Chromosome {chr_index} is not loaded!')

        print(f"\t{self.imported_protein_coding_genes} genes loaded "
              f"({self.ignored_protein_coding_genes} filtered out) "
              f"from chromosomes {loaded_chromosomes}/{self.chromosomes_count()}")
        print(f"\t\tFiltered out Genes: {str(self.ignored_genes_by_types)}")

        if self.annotation_load_type == ANNOTATION_LOAD.GENES:
            return

        features_generator = features_db.features_of_type('mRNA')
        features_mRNA = list(features_generator)
        loaded_transcripts = 0
        for mRNA in features_mRNA:
            if self.__load_feature_by_type(self.__gene_transcripts, mRNA, 'mRNA'): loaded_transcripts += 1
        print(f"\t{loaded_transcripts} mRNAs loaded successfully!")

        if self.annotation_load_type == ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS:
            return

        # load preferred fragments and link them to transcripts in dictionaries
        features_generator = features_db.features_of_type('CDS')
        features_CDSs = list(features_generator)
        loaded_CDSs = 0
        for CDS in features_CDSs:
            if self.__load_feature_by_type(self.__transcript_fragments, CDS, 'CDS'): loaded_CDSs += 1
        print(f"\t{loaded_CDSs} CDSs loaded successfully!")

        if self.annotation_load_type == ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_CDS:
            return

        features_generator = features_db.features_of_type('exon')
        features_exons = list(features_generator)
        loaded_exons = 0
        for exon in features_exons:
            if self.__load_feature_by_type(self.__transcript_fragments, exon, 'exon'): loaded_exons += 1
        print(f"\t{loaded_exons} EXONs loaded successfully!")

        features_generator = features_db.features_of_type('three_prime_UTR')
        features_UTR3s = list(features_generator)
        loaded_utr3s = 0
        for utr3 in features_UTR3s:
            if self.__load_feature_by_type(self.__transcript_fragments, utr3, 'three_prime_UTR'): loaded_utr3s += 1
        print(f"\t{loaded_utr3s} UTR3s loaded successfully!")

        features_generator = features_db.features_of_type('five_prime_UTR')
        features_UTR5s = list(features_generator)
        loaded_utr5s = 0
        for utr5 in features_UTR5s:
            if self.__load_feature_by_type(self.__transcript_fragments, utr5, 'five_prime_UTR'): loaded_utr5s += 1
        print(f"\t{loaded_utr5s} UTR5s loaded successfully!")

    def __load_feature_by_type(self, dict_to_fill, feature: Feature, feature_type):
        if feature.featuretype != feature_type: return False

        assert feature.attributes.__contains__('Parent') and len(feature.attributes['Parent']) == 1
        assert len(feature.attributes['Parent']) == 1

        # parent_gene_id for mRNA and parent_transcript_id for CDS/UTR5'/UTR3'/Exon
        parent_id = feature.attributes['Parent'][0]

        # for mRNA parent is gene, and for CDS,Exon,UTR's it is transcript.
        parent_is_loaded = self.__loaded_feature_by_id.__contains__(parent_id)

        if not parent_is_loaded:
            return False  # it seems feature belongs to transcript or gene, which was filtered out

        if not dict_to_fill.__contains__(parent_id):
            dict_to_fill[parent_id] = []

        dict_to_fill[parent_id].append(feature)

        assert not self.__loaded_feature_by_id.__contains__(feature.id)

        self.__loaded_feature_by_id[feature.id] = feature

        return True

    def seq_id2_chr_index(self, chr_name):
        if not self.chromosome_name2index.__contains__(chr_name):
            return -1
        assert self.chromosome_name2index.__contains__(chr_name)
        return self.chromosome_name2index[chr_name]

    def chr_index2_seq_id(self, index):
        for chr_label, ind in self.chromosome_name2index.items():
            if ind == index: return chr_label
        assert False

    # endregion

    # region choose representative transcript by criteria

    def __get_transcript_score_by_criteria(self, mRNA_transcript: Feature, criteria: TRANSCRIPT_CRITERIA):
        if criteria == TRANSCRIPT_CRITERIA.LONGEST:
            return mRNA_transcript.end - mRNA_transcript.start + 1

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

    def get_transcript_from_gene_by_criteria(self, gene_id, criteria, tie_breaker_criteria) -> Feature:
        if self.__gene_transcript_by_criteria.__contains__(gene_id):
            # Note: criteria must be same for any call made on single object
            return self.__gene_transcript_by_criteria[gene_id]

        if not self.__gene_transcripts.__contains__(gene_id):
            # there is no valid mRNA transcript fot this specific gene
            return None

        mRNA_transcripts = self.__gene_transcripts[gene_id]
        assert len(mRNA_transcripts) != 0

        # if criteria is RANDOM, just choose if randomly from the lists of all transcripts
        if criteria == TRANSCRIPT_CRITERIA.RANDOM:
            return random.choice(mRNA_transcripts)

        best_isoform = None
        best_score = 0

        for transcript in mRNA_transcripts:
            score = self.__get_transcript_score_by_criteria(transcript, criteria)
            if best_isoform is None or score > best_score:
                best_isoform = transcript
                best_score = score
            elif score == best_score:
                score_2nd = self.__get_transcript_score_by_criteria(transcript, tie_breaker_criteria)
                score_old = self.__get_transcript_score_by_criteria(best_isoform, tie_breaker_criteria)
                if score_2nd > score_old:
                    best_isoform = transcript
                    best_score = score

        # save in dict, so accessing later would not require additional findings/calculations
        self.__gene_transcript_by_criteria[gene_id] = best_isoform

        return best_isoform

    # endregion

    # region sequence retriever methods

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

    # region overlap findings between transcript coding sequences (CDS)

    def __calculate_frame_of_fragment_interval(self, fragment: Feature, interval):
        if fragment.strand == '+':
            return (3 + int(fragment.frame) - (interval[0] - fragment.start) % 3) % 3
        return (3 + int(fragment.frame) - (fragment.end - interval[1]) % 3) % 3

    def __attach_frames_to_interval(self, fragment_a: Feature, fragment_b: Feature, interval):
        frame1 = self.__calculate_frame_of_fragment_interval(fragment_a, interval)
        frame2 = self.__calculate_frame_of_fragment_interval(fragment_b, interval)
        return interval, (frame1, frame2)

        # between each pairs of fragments from each transcript (mRNA) finds fragments overlapped by CDS

    def _are_features_same_framed(self, fragment_a: Feature, fragment_b: Feature) -> bool:
        l = fragment_a.start + int(fragment_a.frame) if fragment_a.strand == '+' else fragment_a.end - int(
            fragment_a.frame)
        r = fragment_b.start + int(fragment_b.frame) if fragment_b.strand == '+' else fragment_b.end - int(
            fragment_b.frame)
        return l % 3 == r % 3

    def get_overlaps_between_transcripts(self, transcript1_id, transcript2_id):
        non_ati_overlaps = []

        fragments_a = self.get_fragments_from_transcript(transcript1_id)
        fragments_b = self.get_fragments_from_transcript(transcript2_id)

        if len(fragments_b) == 0 or len(fragments_b) == 0:
            return non_ati_overlaps

        contains_in_phase_overlaps = False

        ov_length = 0
        for fragment_a in fragments_a:
            if fragment_a.featuretype != 'CDS': continue
            for fragment_b in fragments_b:
                if fragment_b.featuretype != 'CDS': continue
                if fragment_b.end < fragment_a.start or fragment_b.start > fragment_a.end: continue
                if fragment_b.end <= fragment_a.end:
                    if fragment_b.start >= fragment_a.start:
                        overlap = (fragment_b.start, fragment_b.end)
                    else:
                        overlap = (fragment_a.start, fragment_b.end,)
                else:
                    if fragment_b.start <= fragment_a.start:
                        overlap = (fragment_a.start, fragment_a.end)
                    else:
                        overlap = (fragment_b.start, fragment_a.end)

                if fragment_a.strand == fragment_b.strand and self._are_features_same_framed(fragment_a, fragment_b):
                    contains_in_phase_overlaps = True
                    continue

                overlap = self.__attach_frames_to_interval(fragment_a, fragment_b, overlap)
                ov_length += overlap[0][1] - overlap[0][0] + 1
                non_ati_overlaps.append(overlap)

        return non_ati_overlaps, ov_length, contains_in_phase_overlaps

    # endregion

    # region Analyze Gene nucleotide composition(stranded) by fragments

    # for a transcript returns merged sequences of UTR's, CDS's, introns and tails.
    def __get_regional_merged_sequences_from_transcript(self, chr_id, transcript, tail_length):
        fragments = self.get_fragments_from_transcript(transcript.id)

        utr5_sequence = ""
        utr3_sequence = ""
        cds_sequence = ""
        intron_sequence = ""
        upstream_sequence = ""
        downstream_sequence = ""

        if len(fragments) == 0: return [upstream_sequence, utr5_sequence, cds_sequence, intron_sequence,
                                        utr3_sequence,
                                        downstream_sequence]

        for fragment in fragments:
            if fragment.featuretype == 'five_prime_UTR':
                if transcript.strand == '+':
                    utr5_sequence += self.retrieve_feature_sequence(chr_id, fragment)
                else:
                    utr5_sequence = self.retrieve_feature_sequence(chr_id, fragment) + utr5_sequence
            if fragment.featuretype == 'three_prime_UTR':
                if transcript.strand == '+':
                    utr3_sequence += self.retrieve_feature_sequence(chr_id, fragment)
                else:
                    utr3_sequence = self.retrieve_feature_sequence(chr_id, fragment) + utr3_sequence

            if fragment.featuretype == 'CDS':
                if transcript.strand == '+':
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
            if transcript.strand == '+':
                intron_sequence += self.retrieve_interval_sequence(chr_id, last_exon.end + 1, fragment.start - 1,
                                                                   fragment.strand)
            else:
                intron_sequence = self.retrieve_interval_sequence(chr_id, last_exon.end + 1, fragment.start - 1,
                                                                  fragment.strand) + intron_sequence
            introns_len += fragment.start - last_exon.end - 1
            last_exon = fragment

        if left_fragment.strand == '+':
            upstream_sequence = self.retrieve_interval_sequence(chr_id, left_fragment.start - tail_length,
                                                                left_fragment.start - 1, '+')

            downstream_sequence = self.retrieve_interval_sequence(chr_id, right_fragment.end + 1,
                                                                  right_fragment.end + tail_length, '+')

        else:
            downstream_sequence = self.retrieve_interval_sequence(chr_id, left_fragment.start - tail_length,
                                                                  left_fragment.start - 1, '-')

            upstream_sequence = self.retrieve_interval_sequence(chr_id, right_fragment.end + 1,
                                                                right_fragment.end + tail_length, '-')

        return [upstream_sequence, utr5_sequence, cds_sequence, intron_sequence, utr3_sequence,
                downstream_sequence]

    # subregions of interest: upstream tail, UTR'5, inner CDSs, inner Introns, UTR'3, downstream tail
    # each region divided into k-part
    # for each part (a,c,g,t) is calculated
    # occurrences[region][part][base] = [6][k][4]
    def __get_transcript_occurrences_by_parts(self, chr_id, transcript: Feature, k, tail_length):
        regional_sequences = self.__get_regional_merged_sequences_from_transcript(chr_id, transcript, tail_length)

        upstream_sequence = regional_sequences[0]
        utr5_sequence = regional_sequences[1]
        cds_sequence = regional_sequences[2]
        intron_sequence = regional_sequences[3]
        utr3_sequence = regional_sequences[4]
        downstream_sequence = regional_sequences[5]

        occurrences = [self.sequence_composition_by_parts(upstream_sequence, k),
                       self.sequence_composition_by_parts(utr5_sequence, k),
                       self.sequence_composition_by_parts(cds_sequence, k),
                       self.sequence_composition_by_parts(intron_sequence, k),
                       self.sequence_composition_by_parts(utr3_sequence, k),
                       self.sequence_composition_by_parts(downstream_sequence, k)]
        return occurrences

    def add_matrix(self, base, new, parts):
        if len(base) == 0: return new
        for i in range(0, 6):
            for j in range(0, parts):
                for k in range(0, 4):
                    base[i][j][k] += new[i][j][k]
        return base

    # output: occurrences[region][part][base] = [6][k][4]
    def analyze_gene_occurrences_by_parts(self, chr_id, gene: Feature, k, tail_length,
                                          criteria=TRANSCRIPT_CRITERIA.NONE):
        if criteria == TRANSCRIPT_CRITERIA.NONE:
            mRNA_transcripts = self.__gene_transcripts[gene.id] if self.__gene_transcripts.__contains__(gene.id) else []
        else:
            transcript_by_criteria = self.get_gene_transcript_by_critereia(gene, criteria)
            mRNA_transcripts = [] if transcript_by_criteria is None else [transcript_by_criteria]

        all_transcript_occurrences = []

        for transcript in mRNA_transcripts:
            occurrences = self.__get_transcript_occurrences_by_parts(chr_id, transcript, k, tail_length)
            all_transcript_occurrences = self.add_matrix(all_transcript_occurrences, occurrences, k)

        return all_transcript_occurrences

    # endregion

    # region static Methods

    @staticmethod
    def get_gene_accession(gene: Feature):
        if not gene.attributes.__contains__('description'): return "no_acc"
        desc_parts = gene.attributes['description'][0].split("Acc:")
        if len(desc_parts) != 2: return "no_acc"  # it seems record does not have NCBI ID
        ncbi_id = desc_parts[1][0:(len(desc_parts[1]) - 1)]
        return ncbi_id

    @staticmethod
    def get_gene_symbol(gene: Feature):
        if not gene.attributes.__contains__('Name'): return "no_sym"  # loaded genes must be filtered
        assert len(gene.attributes['Name']) == 1
        return gene.attributes['Name'][0].upper()

    @staticmethod
    def get_gene_description(gene: Feature):
        if not gene.attributes.__contains__('description'): return "no_desc"  # loaded genes must be filtered
        return gene.attributes['description'][0]

    @staticmethod
    def are_segments_overlapped(segment1, segment2) -> bool:
        l_a, r_a, l_b, r_b = segment1[0], segment1[1], segment2[0], segment2[1]
        return l_a <= l_b <= r_a or l_a <= r_b <= r_a or l_b <= l_a <= r_b

    def are_genes_presented(self, sym1, sym2):
        return self.__gene_symbols_set.__contains__(sym1) and self.__gene_symbols_set.__contains__(sym2)

    def are_genes_overlapped(self, sym1, sym2):
        if not self.are_genes_presented(sym1, sym2): return False
        if not self.__gene_symbols_set.__contains__(sym1) or not self.__gene_symbols_set.__contains__(sym2):
            return False
        gene1 = self.feature_by_id(self.__gene_symbols_set[sym1])
        gene2 = self.feature_by_id(self.__gene_symbols_set[sym2])
        if gene1.chrom != gene2.chrom: return False
        return self.are_segments_overlapped((gene1.start, gene1.end), (gene2.start, gene2.end))

    @staticmethod
    def get_overlap_length(gene_feature_A: Feature, gene_feature_B: Feature):
        if gene_feature_A.chrom != gene_feature_B.chrom: return 0
        l_a, r_a = gene_feature_A.start, gene_feature_A.end
        l_b, r_b = gene_feature_B.start, gene_feature_B.end
        if l_a <= l_b <= r_a: return min(r_a, r_b) - l_b + 1
        if l_a <= r_b <= r_a: return r_b - max(l_a, l_b) + 1
        if l_b <= l_a <= r_b: return r_a - l_a + 1
        return 0

    @staticmethod
    def get_features_distance(segment1: Feature, segment2: Feature):
        if GenomeWorker.get_overlap_type(segment1, segment2) != OVERLAP_TYPE.NONE: return 0
        if segment1.chrom != segment2.chrom: return 1000
        if segment1.end > segment2.end:
            return segment1.start - segment2.end
        else:
            return segment2.start - segment1.end

    @staticmethod
    def get_overlap_type(segment1: Feature, segment2: Feature) -> OVERLAP_TYPE:
        assert segment1 is not None and segment2 is not None

        l_a = segment1.start
        r_a = segment1.end
        l_b = segment2.start
        r_b = segment2.end

        if (l_a <= l_b <= r_a and l_a <= r_b <= r_a) or (l_b <= l_a <= r_b and l_b <= r_a <= r_b):
            return OVERLAP_TYPE.SAME_NESTED if segment1.strand == segment2.strand else OVERLAP_TYPE.DIFF_NESTED

        if l_a <= l_b <= r_a or l_a <= r_b <= r_a:
            if segment1.strand == segment2.strand:
                return OVERLAP_TYPE.TANDEM
            if l_a <= l_b <= r_a and segment1.strand == '+':
                return OVERLAP_TYPE.CONVERGENT
            if l_a <= l_b <= r_a and segment1.strand == '-':
                return OVERLAP_TYPE.DIVERGENT
            if segment1.strand == '+': return OVERLAP_TYPE.DIVERGENT
            return OVERLAP_TYPE.CONVERGENT

        return OVERLAP_TYPE.NONE

    # retrieves nucleotide composition of sequence
    # output: (C_count,G_count,A_count,T_count)
    @staticmethod
    def sequence_composition(sequence):
        sequence = sequence.upper()
        stats = [0, 0, 0, 0]
        for char in sequence:
            if char == 'C':
                stats[0] += 1
            elif char == 'G':
                stats[1] += 1
            elif char == 'A':
                stats[2] += 1
            elif char == 'T':
                stats[3] += 1
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
