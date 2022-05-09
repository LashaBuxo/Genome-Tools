import gc
import random
from os.path import exists
import time
import gffutils
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from gffutils import Feature
import random

import xlrd
from openpyxl import load_workbook
import urllib.request, json
# Author: "Lasha Bukhnikashvili"
#
# Common Methods for working with Genome Worker data:
#   1) Creating database based on NCBI, Ensembl Genome Worker annotation
#   2) Import annotation features from created database
#   3) Filter specific features
#   4) calculate overlapping intervals
#   and so on...
#
# Other codes, which does specific calculations on Genome Worker use this script
# for working with Ensembl or NCBI annotation (with corresponding sequence data if necessary)

# Specific tools for working with Ensembl annotation and with sequence data
from genome_worker_values import *
from genome_worker_enums import *


class GenomeWorker:
    def __init__(self, species: SPECIES, annotation_source: ANNOTATIONS, annotation_load_type: ANNOTATION_LOAD,
                 sequence_load_type: SEQUENCE_LOAD, expression_load_type: EXPRESSIONS_LOAD):

        assert annotation_source == ANNOTATIONS.ENSEMBL or annotation_source == ANNOTATIONS.NCBI

        print(f"loading data for {species} ({annotation_source.name}):")

        self.species = species
        self.annotation_source = annotation_source
        self.annotation_load_type = annotation_load_type
        self.sequence_load_type = sequence_load_type
        self.expression_load_type = expression_load_type

        self.imported_protein_coding_genes = 0
        self.ignored_protein_coding_genes = 0

        self.__loaded_feature_by_id = {}  # storage for gene or transcript features to access by id

        # gene features clustered by chromosome, they are on.
        self.__genes_on_chr = [None] * (self.chromosomes_count() + 1)
        self.__gene_symbols_set = {}
        self.__gene_accessions_set = {}

        self.__sequences = [None] * (self.chromosomes_count() + 1)  # {chr_id,DNA sequence of chromosome with id chr_id}

        self.__gene_transcripts = {}  # {gene_id, [list of mRNA features/isoforms whose parent is gene_id]}

        self.__transcript_tissue_expressions = {}  # {transcript_id, [list of 'CDS' and/or 'UTR(3'/5')' and/or exons]]

        self.__transcript_fragments = {}
        self.__transcript_fragments_is_sorted = {}  # optimizes speed, if we are not sorting it if we don't need

        self.__gene_transcript_by_criteria = {}  # {gene_id, representative mRNA features whose parent is gene_id}

        self.__load_requested_features()
        self.__load_requested_expressions()
        self.__load_requested_sequences()

    # region easy retriever methods

    def get_transcript_median_expression(self, transcript_id):
        expression_data = self.get_transcript_expression_values(transcript_id)

    def get_transcript_expression_values(self, transcript_id):
        assert self.annotation_source == ANNOTATIONS.ENSEMBL  # gene expression data is only available for Ensembl
        assert self.species == SPECIES.Homo_sapiens  # gene expression data is only available for Homo Sapiens

        if not self.__transcript_tissue_expressions.__contains__(transcript_id):
            return {}

        return self.__transcript_tissue_expressions[transcript_id]

    def chromosomes_count(self):
        return NUMBER_OF_CHROMOSOMES[self.species.value]

    def genes_count_on_chr(self, chr_id):
        assert chr_id <= self.chromosomes_count()
        return len(self.__genes_on_chr[chr_id])

    def gene_by_indexes(self, chr_id, index) -> Feature:
        assert chr_id <= self.chromosomes_count()
        return self.__genes_on_chr[chr_id][index]

    def feature_by_id(self, feature_id) -> Feature:
        return self.__loaded_feature_by_id[feature_id]

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

    def get_transcripts_from_gene(self, gene_id) -> list[Feature]:
        return self.__gene_transcripts[gene_id] if self.__gene_transcripts.__contains__(gene_id) else []

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

    # region load methods for Annotation, Assembly(sequence) and Expression values

    sample_tissues = []

    def __load_requested_expressions(self):
        if self.expression_load_type != EXPRESSIONS_LOAD.LOAD:
            return

        assert exists(EXPRESSION_DATA_PATH)

        file = open(EXPRESSION_DATA_PATH, 'r')
        headers = file.readline().split('\t')

        for header in headers:
            if header != 'ensgid' and header != 'enstid' and header != 'gene_name' and header != 'transcript_biotype':
                sample_name = header.split('.')[1]
                self.sample_tissues.append(sample_name)

        expression_records = file.readlines()

        for record_line in expression_records:
            record_data = record_line.replace('\n', '').split('\t')
            if record_data[3] != 'protein_coding': continue

            transcript_id = f"transcript:{record_data[0]}"  # convert transcript to Ensembl format id
            if not self.__transcript_fragments.__contains__(transcript_id): continue

            transcript_expression = {}
            for i in range(0, len(self.sample_tissues)):
                tissue_name = self.sample_tissues[i].split('_')[0]
                if not transcript_expression.__contains__(tissue_name):
                    transcript_expression[tissue_name] = 0
                value = float(record_data[i + 4]) if record_data[i + 4] != 'NA' else 0
                transcript_expression[tissue_name] += value
            self.__transcript_tissue_expressions[transcript_id] = transcript_expression

        print(f"\t{len(self.__transcript_tissue_expressions)} mRNA's expression values loaded successfully!")

    def __load_requested_sequences(self):
        if self.sequence_load_type != SEQUENCE_LOAD.LOAD:
            return

        sequence_file_path = self.__get_sequence_file_path()
        for record in SeqIO.parse(sequence_file_path, 'fasta'):
            chr_id = self.chr_id_from_seq_id(self.annotation_source, record.id)
            if chr_id == -1: continue
            self.__sequences[chr_id] = record
        loaded_chromosomes = 0
        for chr_id in range(1, NUMBER_OF_CHROMOSOMES[self.species.value] + 1):
            loaded_chromosomes += 1 if self.__sequences[chr_id] is not None and len(self.__sequences[chr_id]) > 0 else 0

        print(f"\tDNA chromosomes loaded {loaded_chromosomes}/{NUMBER_OF_CHROMOSOMES[self.species.value]}!")

    def __check_gene_for_Ensembl_filters(self, gene: Feature):
        # ignore genes without description
        if not gene.attributes.__contains__('description'): return (False, "no_desc")

        # ignore genes if they are novel, predicted or readthrough
        description = gene.attributes['description'][0]
        if description.__contains__('novel') or description.__contains__('Novel'): return (False, "is_novel")
        if description.__contains__('predicted') or description.__contains__('Predicted'): return (
            False, "is_predicted")
        if description.__contains__('readthrough'): return (False, "is_readthrough")

        # ignore genes with duplicate Names if it has Name
        gene_symbol = GenomeWorker.get_gene_symbol(gene)
        if gene_symbol != 'no_name' and self.__gene_symbols_set.__contains__(gene_symbol):
            return (False, "name_duplicated")

        # ignore genes if gene with same NCBI accession already imported
        gene_accession = GenomeWorker.get_gene_accession(gene)
        if gene_accession != 'no_acc' and self.__gene_accessions_set.__contains__(gene_accession):
            return (False, "acc_duplicated")

        return (True, "passed_filter")

    def __check_gene_for_NCBI_filters(self, gene: Feature):
        # ignore genes without Name
        if not gene.attributes.__contains__('Name'): return (False, "no_name")

        # ignore genes without description
        if not gene.attributes.__contains__('description'): return (False, "no_desc")

        # ignore genes if they are novel, predicted or readthrough
        description = gene.attributes['description'][0]
        if description.__contains__('novel') or description.__contains__('Novel'): return (False, "is_novel")
        if description.__contains__('predicted') or description.__contains__('Predicted'): return (
            False, "is_predicted")
        if description.__contains__('readthrough'): return (False, "is_readthrough")

        # ignore genes with duplicate Names if it has Name
        gene_name = GenomeWorker.get_gene_symbol(gene)
        if gene_name != 'no_name' and self.__gene_symbols_set.__contains__(gene_name):
            return (False, "name_duplicated")

        # ignore genes if gene with same NCBI accession already imported
        gene_accession = GenomeWorker.get_gene_accession(gene)
        if gene_accession != 'no_acc' and self.__gene_accessions_set.__contains__(gene_accession):
            return (False, "acc_duplicated")

        return (True, "passed_filter")

    def __check_gene_for_filters(self, gene: Feature):
        # https: // www.biostars.org / p / 5304 /  # 9521037
        if self.annotation_source == ANNOTATIONS.NCBI:
            return self.__check_gene_for_NCBI_filters(gene)
        else:
            return self.__check_gene_for_Ensembl_filters(gene)

    # if database does not exists for specific chromosome, then
    # builds database from annotation file and fills/loads all
    # necessary features from the database.
    #
    # if database exists, then fills/loads all necessary data structures
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

        self.imported_protein_coding_genes = 0
        self.ignored_protein_coding_genes = 0

        ignored_genes_by_types = {}
        # choose genes who has protein_coding attribute and additional filter values
        for gene in feature_genes:
            attribute_filter = 'gene_biotype' if self.annotation_source == ANNOTATIONS.NCBI else 'biotype'
            attribute_value = 'protein_coding'

            if gene.attributes.__contains__(attribute_filter):
                assert len(gene.attributes[attribute_filter]) == 1
                if gene.attributes[attribute_filter][0] == attribute_value:
                    chr_id = self.chr_id_from_seq_id(self.annotation_source, gene.chrom)
                    if chr_id != -1:
                        is_valid, status = self.__check_gene_for_filters(gene)
                        if is_valid:
                            if self.__genes_on_chr[chr_id] is None:
                                self.__genes_on_chr[chr_id] = []

                            self.imported_protein_coding_genes += 1

                            self.__genes_on_chr[chr_id].append(gene)  # store gene feature on chromosome list
                            self.__loaded_feature_by_id[gene.id] = gene  # store gene feature to access by id

                            # store these values, for filtering upcoming genes
                            self.__gene_symbols_set[GenomeWorker.get_gene_symbol(gene)] = gene.id
                            self.__gene_accessions_set[GenomeWorker.get_gene_accession(gene)] = gene.id
                        else:
                            if not ignored_genes_by_types.__contains__(status):
                                ignored_genes_by_types[status] = 0
                            ignored_genes_by_types[status] += 1
                            self.ignored_protein_coding_genes += 1

        loaded_chromosomes = 0
        for chr_id in range(1, NUMBER_OF_CHROMOSOMES[self.species.value] + 1):
            if self.__genes_on_chr[chr_id] is not None and len(self.__genes_on_chr[chr_id]) > 0:
                loaded_chromosomes += 1
            else:
                print(f'WARNING: Chromosome {chr_id} is not loaded!')

        print(f"\t{self.imported_protein_coding_genes} genes loaded "
              f"({self.ignored_protein_coding_genes} filtered out) "
              f"from chromosomes {loaded_chromosomes}/{NUMBER_OF_CHROMOSOMES[self.species.value]}")
        print(f"\t\tFiltered out Genes: {str(ignored_genes_by_types)}")

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
        for utr3 in features_UTR3s:
            self.__load_feature_by_type(self.__transcript_fragments, utr3, 'three_prime_UTR')
        print(f"\tUTR3s loaded successfully!")

        features_generator = features_db.features_of_type('five_prime_UTR')
        features_UTR5s = list(features_generator)
        for utr5 in features_UTR5s:
            self.__load_feature_by_type(self.__transcript_fragments, utr5, 'five_prime_UTR')
        print(f"\tUTR5s loaded successfully!")

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
        self.__loaded_feature_by_id[feature.id] = feature

        return True

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
                if self.species == SPECIES.Drosophila_melanogaster:
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

    # region choose representative transcript by criteria

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
        return (3 + int(fragment.frame) - (interval[0] - fragment.start) % 3) % 3

    def __attach_frames_to_interval(self, fragment_a: Feature, fragment_b: Feature, interval):
        frame1 = self.__calculate_frame_of_fragment_interval(fragment_a, interval)
        frame2 = self.__calculate_frame_of_fragment_interval(fragment_b, interval)
        return interval, (frame1, frame2)

        # between each pairs of fragments from each transcript (mRNA) finds fragments overlapped by CDS

    def get_fragments_overlap_between_transcripts(self, transcript1_id, transcript2_id):
        overlaps = []

        fragments_a = self.get_fragments_from_transcript(transcript1_id)
        fragments_b = self.get_fragments_from_transcript(transcript2_id)

        if len(fragments_b) == 0 or len(fragments_b) == 0:
            return overlaps

        for fragment_a in fragments_a:
            assert fragment_a.featuretype == 'CDS'
            for fragment_b in fragments_b:
                assert fragment_b.featuretype == 'CDS'

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

                overlap = self.__attach_frames_to_interval(fragment_a, fragment_b, overlap)
                overlaps.append(overlap)

        return overlaps

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
            transcript_by_criteria = self.find_gene_transcript_by_criteria(gene, criteria)
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
        if not gene.attributes.__contains__('description'): return "***"
        desc_parts = gene.attributes['description'][0].split("Acc:")
        if len(desc_parts) != 2: return "no_acc"  # it seems record does not have NCBI ID
        ncbi_id = desc_parts[1][0:(len(desc_parts[1]) - 1)]
        return ncbi_id

    @staticmethod
    def get_gene_symbol(gene: Feature):
        if not gene.attributes.__contains__('Name'): return "***"  # loaded genes must be filtered
        assert len(gene.attributes['Name']) == 1
        return gene.attributes['Name'][0]

    @staticmethod
    def get_gene_description(gene: Feature):
        if not gene.attributes.__contains__('description'): return "***"  # loaded genes must be filtered
        return gene.attributes['description'][0]

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
                return OVERLAP_TYPE.SAME_TANDEM
            if l_a <= l_b <= r_a and segment1.strand == '+':
                return OVERLAP_TYPE.DIFF_CONVERGENT
            if l_a <= l_b <= r_a and segment1.strand == '-':
                return OVERLAP_TYPE.DIFF_DIVERGENT
            if segment1.strand == '+': return OVERLAP_TYPE.DIFF_DIVERGENT
            return OVERLAP_TYPE.DIFF_CONVERGENT

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
