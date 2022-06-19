import gc
import random
from math import log2
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
from worker_genome_values import *
from worker_genome_enums import *
from worker_analyzer import *


class GenomeWorker:
    def __init__(self, species: SPECIES, annotation_source: ANNOTATIONS, annotation_load_type: ANNOTATION_LOAD,
                 sequence_load_type: SEQUENCE_LOAD):

        assert annotation_source == ANNOTATIONS.ENSEMBL or annotation_source == ANNOTATIONS.NCBI

        print(f"loading data for {species} ({annotation_source.name}):")

        self.species = species
        self.annotation_source = annotation_source
        self.annotation_load_type = annotation_load_type
        self.sequence_load_type = sequence_load_type

        self.imported_protein_coding_genes = 0
        self.ignored_protein_coding_genes = 0
        self.ignored_genes_by_types = {}

        # storage for gene or transcript features to access by id
        self.__loaded_feature_by_id = {}

        # gene features clustered by chromosome, they are on.
        self.__genes_on_chr = [None] * (self.chromosomes_count() + 1)
        self.__gene_symbols_set = {}
        self.__gene_accessions_set = {}

        self.__gene_transcripts = {}  # {gene_id, [list of mRNA features/isoforms whose parent is gene_id]}

        self.__transcript_APPRIS_data = {}
        self.__gene_mitoCarta_data = {}

        self.__transcript_fragments = {}
        self.__transcript_fragments_is_sorted = {}  # optimizes speed, if we are not sorting it if we don't need

        self.__gene_transcript_by_criteria = {}  # {gene_id, representative mRNA features whose parent is gene_id}

        # {chr_id,DNA sequence of chromosome with id chr_id}
        self.__sequences = [None] * (self.chromosomes_count() + 1)

        self.__load_requested_features()
        self.__load_mitoCarta_data()

        # we only have APPRIS data for Ensembl
        if annotation_source == ANNOTATIONS.ENSEMBL: self.__load_APPRIS_data()

        self.__load_requested_sequences()

    # region easy retriever methods

    def get_chromosome_length(self, chr_id):
        return len(self.retrieve_sequence_record(chr_id))

    def get_feature_chromosomal_position(self, feature_id):
        feature = self.feature_by_id(feature_id)
        chr_id = self.chr_id_from_seq_id(self.annotation_source, feature.chrom)

        if self.sequence_load_type == SEQUENCE_LOAD.NOT_LOAD:
            return chr_id, ''
        chr_len = len(self.retrieve_sequence_record(chr_id))
        return chr_id, f"{chr_id} ({int(feature.start / chr_len * 100)}%)"

    def chromosomes_count(self):
        return NUMBER_OF_CHROMOSOMES[self.species.value]

    def chromosome_length(self, chr_id):
        return len(self.__sequences[chr_id])

    def genes_count_on_chr(self, chr_id):
        assert chr_id <= self.chromosomes_count()
        return 0 if self.__genes_on_chr[chr_id] is None else len(self.__genes_on_chr[chr_id])

    def gene_by_indexes(self, chr_id, index) -> Feature:
        assert chr_id <= self.chromosomes_count()
        return self.__genes_on_chr[chr_id][index]

    def get_transcript_parent(self, transcript_id):
        feature = self.feature_by_id(transcript_id)
        return self.feature_by_id(feature.attributes['Parent'][0])

    def feature_by_id(self, feature_id) -> Feature:
        return self.__loaded_feature_by_id[feature_id] if self.__loaded_feature_by_id.__contains__(feature_id) else None

    def __get_annotation_file_path(self):
        if self.annotation_source != ANNOTATIONS.ENSEMBL:
            return NCBI_ANNOTATIONS[self.species.value]
        return ENSEMBL_ANNOTATIONS[self.species.value]

    def __get_annotation_db_path(self):
        added_name = self.__get_annotation_file_path().replace('/', '_')
        generated_path = GENOME_DATABASES_DIRECTORY + self.species.name + "_" + self.annotation_source.name + "_" + added_name
        return generated_path + '.db'

    def __get_sequence_file_path(self):
        if self.annotation_source == ANNOTATIONS.ENSEMBL:
            return ENSEMBL_SEQUENCES[self.species.value]
        return NCBI_SEQUENCES[self.species.value]

    def get_transcripts_from_gene(self, gene_id) -> list[Feature]:
        return self.__gene_transcripts[gene_id] if self.__gene_transcripts.__contains__(gene_id) else []

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

    def get_transcript_conservation_info(self, transcript_id):
        assert self.annotation_source == ANNOTATIONS.ENSEMBL  # we only have APPRIS data for Ensembl
        transcript_id = transcript_id.replace('transcript:', '')
        if not self.__transcript_APPRIS_data.__contains__(transcript_id): return 0
        score = float(self.__transcript_APPRIS_data[transcript_id]['conservation_score'])
        n_score = score
        return n_score

    def get_transcript_homologue_species(self, transcript_id):
        assert self.annotation_source == ANNOTATIONS.ENSEMBL  # we only have APPRIS data for Ensembl
        transcript_id = transcript_id.replace('transcript:', '')
        if not self.__transcript_APPRIS_data.__contains__(transcript_id): return 0
        if not self.__transcript_APPRIS_data[transcript_id].__contains__('homologue'): return []
        return self.__transcript_APPRIS_data[transcript_id]['homologue']

    def get_gene_max_conserved_homologue_species(self, gene_id):
        assert self.annotation_source == ANNOTATIONS.ENSEMBL  # we only have APPRIS data for Ensembl
        transcripts = self.get_transcripts_from_gene(gene_id)
        arrs = []
        for transcript in transcripts:
            arrs.append((transcript.id, self.get_transcript_conservation_info(transcript.id)))
        arrs.sort(key=lambda x: x[1], reverse=True)
        transcript_id = arrs[0][0]
        transcript_id = transcript_id.replace('transcript:', '')
        if not self.__transcript_APPRIS_data.__contains__(transcript_id): return 0
        if not self.__transcript_APPRIS_data[transcript_id].__contains__('homologue'): return []
        return self.__transcript_APPRIS_data[transcript_id]['homologue']

    def is_gene_MITO(self, gene_id):
        feature = self.feature_by_id(gene_id)
        gene_sym = self.get_gene_symbol(feature)
        return self.__gene_mitoCarta_data.__contains__(gene_sym)

    # endregion

    # region load methods for Annotation, Assembly(sequence), APPRIS and mitoCarta

    def __load_mitoCarta_data(self):
        file = open(MITOCARTA_DATA[self.species.value], 'r')
        lines = file.readlines()
        total = 0
        for line in lines:
            if not line.__contains__('<tr><td>'): continue
            total += 1
            gene_name = line.split('<tr><td>')[1].split('</td><td')[0].upper()
            if not self.__gene_symbols_set.__contains__(gene_name):
                continue
            self.__gene_mitoCarta_data[gene_name] = True
        print(f"\t{len(self.__gene_mitoCarta_data)}/{total} genes status loaded from mitoCarta!")

    def __load_APPRIS_data(self):
        scores_file = APPRIS_DATA_DIRECTORY[self.species.value] + 'appris_data.appris.txt'
        assert exists(scores_file)
        x = {}
        file = open(scores_file, 'r')
        lines = file.readlines()
        for index in range(1, len(lines)):
            line_data = lines[index].split('\t')
            transcript_id = line_data[2]

            if self.feature_by_id('transcript:' + transcript_id) is None: continue

            # if self.get_transcript_APPRIS_status(transcript_id) is None: continue
            residues, structure, conservation, domains, helices, signals, trifid_score, mapped_peptides, appris_score, appris_annotation = \
                line_data[10], line_data[11], line_data[12], line_data[13], line_data[14], line_data[15], line_data[16], \
                line_data[17], line_data[18], line_data[19]

            not_found_tag = line_data[6]

            trifid_score = 0 if trifid_score == '' else float(trifid_score)
            if not self.__transcript_APPRIS_data.__contains__(transcript_id):
                self.__transcript_APPRIS_data[transcript_id] = {}
            self.__transcript_APPRIS_data[transcript_id]['functional_residues'] = residues
            self.__transcript_APPRIS_data[transcript_id]['structure_score'] = structure
            self.__transcript_APPRIS_data[transcript_id]['domain_residues'] = domains
            self.__transcript_APPRIS_data[transcript_id]['conservation_score'] = conservation
            self.__transcript_APPRIS_data[transcript_id]['Trifid_Score'] = trifid_score
            self.__transcript_APPRIS_data[transcript_id]['APPRIS_score'] = appris_score
            self.__transcript_APPRIS_data[transcript_id]['APPRIS_annotation'] = appris_annotation
            self.__transcript_APPRIS_data[transcript_id]['transcript_incomplete'] = not_found_tag

            x[not_found_tag] = 1
        print(f"\t{len(self.__transcript_APPRIS_data)}/{len(lines) - 1} transcripts scores loaded from APPRIS!")

        if self.species != SPECIES.Homo_sapiens: return
        corsair_file = APPRIS_DATA_DIRECTORY[self.species.value] + 'corsair.txt'
        assert exists(corsair_file)
        file = open(corsair_file, 'r')
        lines = file.readlines()
        cnt = 0
        transcript_id = ''
        for index in range(0, len(lines)):
            line = lines[index].replace('"', '')
            if line.startswith('>'):
                transcript_id = line.replace('>', '').split('\t')[0]
                continue
            if len(line) < 3: continue
            # if transcript_id is filtered, then ignore record
            if not self.__transcript_APPRIS_data.__contains__(transcript_id):
                continue
            if not self.__transcript_APPRIS_data[transcript_id].__contains__('homologue'):
                self.__transcript_APPRIS_data[transcript_id]['homologue'] = []
            species = line.split('\t')[0]
            self.__transcript_APPRIS_data[transcript_id]['homologue'].append(species)
        print(f"\t{cnt}/{len(lines)} homologue (CORSAIR) cross-species loaded from APPRIS!")

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
        # ignore genes with miscellaneous chromosome/scaffold names
        chr_id = self.chr_id_from_seq_id(self.annotation_source, gene.chrom)
        if chr_id == -1: return False, "located_on_miscellaneous_scaffold"

        # ignore genes without description or gene_symbol
        if not gene.attributes.__contains__('description'): return False, "no_desc"
        gene_symbol = GenomeWorker.get_gene_symbol(gene)
        if gene_symbol == "***": return False, "no_symbol"

        description = gene.attributes['description'][0]

        # ignore genes if they are readthrough
        if description.__contains__('readthrough'): return False, "is_readthrough"

        # ignore genes with duplicate Names or accessions
        if self.__gene_symbols_set.__contains__(gene_symbol):
            return False, "name_duplicated"
        # ignore genes if gene with same NCBI accession already imported
        gene_accession = GenomeWorker.get_gene_accession(gene)
        if gene_accession != 'no_acc' and self.__gene_accessions_set.__contains__(gene_accession):
            return False, "acc_duplicated"

        # ignore genes if they are pseudogene, novel or predicted
        if description.__contains__('pseudogene') or description.__contains__('pseudogene'):
            # print(f"{description} {gene.id}\n")
            return (False, "is_pseudogene")
        if description.__contains__('novel') or description.__contains__('Novel'):
            # print(f"{description} {gene.id}\n")
            return False, "is_novel"
        if description.__contains__('predicted') or description.__contains__('Predicted'):
            # print(f"{description} {gene.id}\n")
            return (False, "is_predicted")

        # for any annotation, there must not be 2 gene record with same id
        assert not self.__loaded_feature_by_id.__contains__(gene.id)

        return True, "passed_filter"

    def __check_gene_for_NCBI_filters(self, gene: Feature):
        # ignore genes with miscellaneous chromosome/scaffold names
        chr_id = self.chr_id_from_seq_id(self.annotation_source, gene.chrom)
        if chr_id == -1: return False, "located_on_miscellaneous_scaffold"

        # ignore genes without Name
        if not gene.attributes.__contains__('Name'): return False, "no_name"

        if chr_id != NUMBER_OF_CHROMOSOMES[self.species.value]:
            # ignore genes without description
            if not gene.attributes.__contains__('description'): return False, "no_desc"

            # ignore genes if they are novel, predicted or readthrough
            description = gene.attributes['description'][0]
            if description.__contains__('pseudogene') or description.__contains__('pseudogene'): return (
                False, "is_pseudogene")
            if description.__contains__('novel') or description.__contains__('Novel'): return False, "is_novel"
            if description.__contains__('predicted') or description.__contains__('Predicted'): return (
                False, "is_predicted")
            if description.__contains__('readthrough'): return False, "is_readthrough"

        # for any annotation, there must not be 2 gene record with same id
        assert not self.__loaded_feature_by_id.__contains__(gene.id)

        # ignore genes with duplicate Names if it has Name
        gene_name = GenomeWorker.get_gene_symbol(gene)
        if gene_name != 'no_name' and self.__gene_symbols_set.__contains__(gene_name):
            return False, "name_duplicated"

        if chr_id != NUMBER_OF_CHROMOSOMES[self.species.value]:
            # ignore genes if gene with same NCBI accession already imported
            gene_accession = GenomeWorker.get_gene_accession(gene)
            if gene_accession != 'no_acc' and self.__gene_accessions_set.__contains__(gene_accession):
                return False, "acc_duplicated"

        return True, "passed_filter"

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
        # choose genes who has protein_coding attribute and additional filter values

        for gene in feature_genes:
            attribute_filter = 'gene_biotype' if self.annotation_source == ANNOTATIONS.NCBI else 'biotype'
            attribute_value = 'protein_coding'
            if gene.attributes.__contains__(attribute_filter):
                assert len(gene.attributes[attribute_filter]) == 1
                if gene.attributes[attribute_filter][0] == attribute_value:
                    is_valid, status = self.__check_gene_for_filters(gene)
                    if is_valid:
                        chr_id = self.chr_id_from_seq_id(self.annotation_source, gene.chrom)
                        if self.__genes_on_chr[chr_id] is None:
                            self.__genes_on_chr[chr_id] = []

                        self.imported_protein_coding_genes += 1

                        self.__genes_on_chr[chr_id].append(gene)  # store gene feature on chromosome list
                        self.__loaded_feature_by_id[gene.id] = gene  # store gene feature to access by id

                        # store these values, for filtering upcoming genes
                        self.__gene_symbols_set[GenomeWorker.get_gene_symbol(gene)] = gene.id
                        self.__gene_accessions_set[GenomeWorker.get_gene_accession(gene)] = gene.id
                    else:
                        if not self.ignored_genes_by_types.__contains__(status):
                            self.ignored_genes_by_types[status] = 0
                        self.ignored_genes_by_types[status] += 1
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

    def chr_id_from_seq_id(self, annotation_source, id):
        if annotation_source == ANNOTATIONS.ENSEMBL:
            if id == 'MT': return NUMBER_OF_CHROMOSOMES[self.species.value]
            if id == 'Y': return NUMBER_OF_CHROMOSOMES[self.species.value] - 1
            if id == 'X': return NUMBER_OF_CHROMOSOMES[self.species.value] - 2

            try:
                index = int(id)
            except ValueError:
                index = -1

            return index if 1 <= index <= NUMBER_OF_CHROMOSOMES[self.species.value] else -1

        # NCBI 1st column rules: accession number instead of chromosome number
        if not id.__contains__('.'): return -1
        parts = id.split('.')

        if not parts[0][0:3] == 'NC_': return -1
        num = parts[0][3:len(parts[0])]
        x = int(num)

        if self.species == SPECIES.Homo_sapiens:
            if x == 12920: return 25
            if x < 1 or x > 24: return -1
            return x
        elif self.species == SPECIES.Mus_musculus:
            if x == 5089: return 22
            base = 66
            x = x - base
            if x < 1 or x > 21: return -1
            return x
        else:
            print("NCBI chrid needs conversation to index!!!!!!!!!!!!!")
            return -1

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
        if not gene.attributes.__contains__('description'): return "***"
        desc_parts = gene.attributes['description'][0].split("Acc:")
        if len(desc_parts) != 2: return "no_acc"  # it seems record does not have NCBI ID
        ncbi_id = desc_parts[1][0:(len(desc_parts[1]) - 1)]
        return ncbi_id

    @staticmethod
    def get_gene_symbol(gene: Feature):
        if not gene.attributes.__contains__('Name'): return "***"  # loaded genes must be filtered
        assert len(gene.attributes['Name']) == 1
        return gene.attributes['Name'][0].upper()

    @staticmethod
    def get_gene_description(gene: Feature):
        if not gene.attributes.__contains__('description'): return "***"  # loaded genes must be filtered
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
