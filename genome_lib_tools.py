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
from os.path import exists
import gffutils
from gffutils import FeatureDB, Feature
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import genome_lib_values as CONSTANTS
import enum


# region ENUM classes
class ANNOTATION(enum.Enum):
    NCBI = 1
    ENSEMBL = 2
    NONE = 3


class ANNOTATION_LOAD(enum.Enum):
    GENES = 1
    GENES_AND_TRANSCRIPTS = 2
    GENES_AND_TRANSCRIPTS_AND_CDS = 3
    GENES_AND_TRANSCRIPTS_AND_FRAGMENTS = 4


class SEQUENCE_LOAD(enum.Enum):
    LOAD = 1
    NOT_LOAD = 2


# endregion

# region Variables and Retrievers


# gene features clustered by chromosome, they are on.
# Example: genes_on_chr[22] = ['gene' features on chromosome 22]
__genes_on_chr = [None] * (CONSTANTS.NUMBER_OF_CHROMOSOMES + 1)

# {chr_id,DNA sequence of chromosome with id chr_id}
__sequences = [None] * (CONSTANTS.NUMBER_OF_CHROMOSOMES + 1)

#   {gene_id, [list of mRNA features whose parent is gene_id]}
__gene_transcripts = {}

#   {mRNA_transcript_id, [list of 'CDS' and/or 'UTR(3'/5')' and/or 'exons' features whose parent is mRNA_transcript_id]}
__transcript_fragments = {}


def chromosomes_count():
    return CONSTANTS.NUMBER_OF_CHROMOSOMES


def genes_count_on_chr(chr_id):
    assert chr_id <= CONSTANTS.NUMBER_OF_CHROMOSOMES
    return len(__genes_on_chr[chr_id])


def gene_by_ind(chr_id, index) -> Feature:
    assert chr_id <= CONSTANTS.NUMBER_OF_CHROMOSOMES
    return __genes_on_chr[chr_id][index]


def __get_annotation_file_path(annotation_source):
    if annotation_source != ANNOTATION.ENSEMBL:
        return CONSTANTS.NCBI_WHOLE_ANNOTATION_FILE_PATH
    return CONSTANTS.ENSEMBL_WHOLE_ANNOTATION_FILE_PATH


def __get_annotation_db_path(annotation_source):
    if annotation_source == ANNOTATION.ENSEMBL:
        return CONSTANTS.ENSEMBL_ANNOTATION_DB_FILE_PATH
    return CONSTANTS.NCBI_ANNOTATION_DB_FILE_PATH


def __get_sequence_file_path(annotation_source):
    if annotation_source == ANNOTATION.ENSEMBL:
        return CONSTANTS.ENSEMBL_WHOLE_SEQUENCE_FILE_PATH
    return CONSTANTS.NCBI_WHOLE_SEQUENCE_FILE_PATH


# endregion

# region Additional Methods

def __better_annotated_transcript(A: Feature, B: Feature) -> Feature:
    if A == None: return B
    if B == None: return A

    #  compare by length?
    if A.end - A.start > B.end - B.start: return A
    if B.end - B.start > A.end - A.start: return B

    tagged_A = A.attributes.__contains__('tag') and A.attributes['tag'][0] == 'basic'
    tagged_B = B.attributes.__contains__('tag') and B.attributes['tag'][0] == 'basic'

    if tagged_A == True and tagged_B == False: return A
    if tagged_B == True and tagged_A == False: return B

    transcript_support_lvl_A = 1000
    transcript_support_lvl_B = 1000

    if A.attributes.__contains__('transcript_support_level') and (
            A.attributes['transcript_support_level'][0][0]).isdigit():
        transcript_support_lvl_A = int(A.attributes['transcript_support_level'][0][0])
    if B.attributes.__contains__('transcript_support_level') and (
            B.attributes['transcript_support_level'][0][0]).isdigit():
        transcript_support_lvl_B = int(B.attributes['transcript_support_level'][0][0])

    if transcript_support_lvl_A > transcript_support_lvl_B: return B
    return A


def are_genes_overlapped(gene_feature_A: Feature, gene_feature_B: Feature) -> bool:
    l_a, r_a = gene_feature_A.start, gene_feature_A.end
    l_b, r_b = gene_feature_B.start, gene_feature_B.end
    if l_a <= l_b <= r_a: return True
    if l_a <= r_b <= r_a: return True
    if l_b <= l_a <= r_b: return True
    return False


def get_genes_overlap(gene_feature_A: Feature, gene_feature_B: Feature):
    l_a, r_a = gene_feature_A.start, gene_feature_A.end
    l_b, r_b = gene_feature_B.start, gene_feature_B.end
    if l_a <= l_b <= r_a: return l_b, min(r_a, r_b)
    if l_a <= r_b <= r_a: return max(l_a, l_b), r_b
    if l_b <= l_a <= r_b: return l_a, l_b
    assert False


def is_same_frame_overlap(fragment_a: Feature, fragment_b: Feature) -> bool:
    l = fragment_a.start + int(fragment_a.frame) if fragment_a.strand == '+' else fragment_a.end - int(fragment_a.frame)
    r = fragment_b.start + int(fragment_b.frame) if fragment_b.strand == '+' else fragment_b.end - int(fragment_b.frame)
    return l % 3 == r % 3


def get_fragments_overlap_segments(fragments_A: list[Feature], fragments_B: list[Feature], ORF_similarity):
    if len(fragments_A) == 0 or len(fragments_B) == 0: return False

    overlaps = []
    for fragment_a in fragments_A:
        for fragment_b in fragments_B:
            if ORF_similarity == 'diff_ORF' and is_same_frame_overlap(fragment_a, fragment_b): continue
            if ORF_similarity == 'same_ORF' and not is_same_frame_overlap(fragment_a, fragment_b): continue

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
    return overlaps


# endregion

# region Annotation Preprocessing methods

# if database does not exists for specific chromosome, then
# builds database from annotation file and fills/loads all
# necessary features from the database.
#
# if database exists, then fills/loads all necessary data structures

def __preprocess_feature_by_type(dict_to_fill, feature: Feature, feature_type):
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


def chr_id_from_seq_id(annotation_source, id):
    if annotation_source == ANNOTATION.ENSEMBL:
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
        print("Lasha some issue found!")
    return x


def preprocess_annotation(annotation_source: ANNOTATION, annotation_load_type: ANNOTATION_LOAD,
                          sequence_load_type: SEQUENCE_LOAD):
    assert annotation_source == ANNOTATION.ENSEMBL or annotation_source == ANNOTATION.NCBI

    if sequence_load_type == SEQUENCE_LOAD.LOAD:
        sequence_file_path = __get_sequence_file_path(annotation_source)
        for record in SeqIO.parse(sequence_file_path, 'fasta'):
            chr_id = chr_id_from_seq_id(annotation_source, record.id)
            if chr_id == -1: continue
            __sequences[chr_id] = record
        print("Using " + str(annotation_source) + ": sequence data loaded successfully!")

    if not exists(__get_annotation_db_path(annotation_source)):
        print("Creating database for " + str(annotation_source) + " only for first RUN it takes that long!")
        gffutils.create_db(__get_annotation_file_path(annotation_source),
                           dbfn=__get_annotation_db_path(annotation_source),
                           verbose=True, force=False, keep_order=False, merge_strategy='create_unique',
                           sort_attribute_values=False)

    features_db = gffutils.FeatureDB(__get_annotation_db_path(annotation_source), keep_order=False)

    features_generator = features_db.features_of_type('gene')
    feature_genes = list(features_generator)

    # choose genes who has protein_coding attribute
    for gene in feature_genes:
        if annotation_source == ANNOTATION.NCBI:
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
                    chr_id = chr_id_from_seq_id(annotation_source, gene.chrom)
                    if chr_id != -1:
                        if __genes_on_chr[chr_id] is None:
                            __genes_on_chr[chr_id] = []
                        __genes_on_chr[chr_id].append(gene)

    print("Using " + str(annotation_source) + ": genes loaded successfully!")

    if annotation_load_type == ANNOTATION_LOAD.GENES:
        return

    # load preferred transcripts and link them to genes in dictionaries
    # preprocess_feature_by_type(__gene_transcripts, feature, 'unconfirmed_transcript')
    # preprocess_feature_by_type(__gene_transcripts, feature, 'V_gene_segment')
    # preprocess_feature_by_type(__gene_transcripts, feature, 'J_gene_segment')
    # preprocess_feature_by_type(__gene_transcripts, feature, 'C_gene_segment')
    features_generator = features_db.features_of_type('mRNA')
    features_mRNA = list(features_generator)
    for mRNA in features_mRNA:
        __preprocess_feature_by_type(__gene_transcripts, mRNA, 'mRNA')

    print("Using " + str(annotation_source) + ": mRNAs loaded successfully!")

    if annotation_load_type == ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS:
        return

    # load preferred fragments and link them to transcripts in dictionaries
    features_generator = features_db.features_of_type('CDS')
    features_CDSs = list(features_generator)
    for CDS in features_CDSs:
        __preprocess_feature_by_type(__transcript_fragments, CDS, 'CDS')

    print("Using " + str(annotation_source) + ": CDSs loaded successfully!")

    if annotation_load_type == ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_CDS:
        return

    features_generator = features_db.features_of_type('exon')
    features_exons = list(features_generator)
    for exon in features_exons:
        __preprocess_feature_by_type(__transcript_fragments, exon, 'exon')
    print("Using " + str(annotation_source) + ": EXONs loaded successfully!")

    features_generator = features_db.features_of_type('three_prime_UTR')
    features_UTR3s = list(features_generator)
    for utr3 in features_UTR3s:
        __preprocess_feature_by_type(__transcript_fragments, utr3, 'three_prime_UTR')
    print("Using " + str(annotation_source) + ": UTR3s loaded successfully!")

    features_generator = features_db.features_of_type('five_prime_UTR')
    features_UTR5s = list(features_generator)
    for utr5 in features_UTR5s:
        __preprocess_feature_by_type(__transcript_fragments, utr5, 'five_prime_UTR')
    print("Using " + str(annotation_source) + ": UTR5s loaded successfully!")


# endregion

# region Get sequence of a Gene

# if sequence for specific chromosome is not loaded then loads it.
# if its already loaded retrieves it.
def retrieve_sequence_record(chr_id) -> SeqRecord:
    if __sequences[chr_id] is None:
        print("Sequence must be loaded during processing! Error...")
    return __sequences[chr_id]


def retrieve_feature_sequence(chr_id, feature: Feature) -> str:
    seq_record = retrieve_sequence_record(chr_id)
    feature_record = seq_record[feature.start - 1:feature.end]
    if feature.strand == '-':
        feature_record = feature_record.reverse_complement()
    return str(feature_record.seq)


def retrieve_interval_sequence(chr_id, start, end, strand) -> str:
    seq_record = retrieve_sequence_record(chr_id)
    interval_record = seq_record[max(0, start - 1):min(end, len(seq_record))]
    if strand == '-': interval_record = interval_record.reverse_complement()
    return str(interval_record.seq)


# endregion

# region Get transcript of a gene


def __find_best_gene_transcript(chr_id, gene_feature: Feature) -> Feature:
    if not __gene_transcripts.__contains__(gene_feature.id):
        # there is no valid mRNA transcript fot this specific gene
        return None

    mRNA_transcripts = __gene_transcripts[gene_feature.id]

    best_mRNA_transcript = None
    mRNA_variants = 0

    for transcript in mRNA_transcripts:
        if not transcript.attributes.__contains__('Parent'):
            print("Lasha some issue found!")
        if not len(transcript.attributes) != 1:
            print("Lasha some issue found!")
        if transcript.attributes['Parent'][0] == gene_feature.id:
            mRNA_variants = mRNA_variants + 1
            best_mRNA_transcript = __better_annotated_transcript(best_mRNA_transcript, transcript)

    # print("mRNA variants for Gene [" + gene_feature.id + "]: " + str(mRNA_variants))
    # print("well annotated mRNA: " + best_mRNA_transcript.id)

    return best_mRNA_transcript


# endregion

# region Get fragments of a gene


def __find_transcript_fragments(chr_id, mRNA_transcript: Feature) -> list[Feature]:
    if mRNA_transcript.id not in __transcript_fragments:
        # there is no valid fragment transcript fot this specific transcript
        return []

    fragment_features = __transcript_fragments[mRNA_transcript.id]

    fragments = []
    for fragment in fragment_features:
        if fragment.attributes.__contains__('Parent'):
            if len(fragment.attributes['Parent']) != 1:
                print("Lasha some issue found!")
            if fragment.attributes['Parent'][0] == mRNA_transcript.id:
                fragments.append(fragment)

    return fragments


def get_fragments_on_gene(chr_id, gene_feature: Feature, force_sorted) -> list[Feature]:
    # mRNA_transcript (start,end) always would fit in parent's (start,end) interval
    mRNA_transcript = __find_best_gene_transcript(chr_id, gene_feature)
    if mRNA_transcript is None:
        # it seems there is no valid mRNA for the gene
        # print('no valid mRNA and exons for a gene')
        return []

    fragments = __find_transcript_fragments(chr_id, mRNA_transcript)
    if not force_sorted:
        return fragments
    fragments.sort(key=lambda x: x.start)
    return fragments


# endregion

# region Get nucleotide Composition(stranded)

# retrieves nucleotide composition of sequence
# output: (C_count,G_count,A_count,T_count)
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


def sequence_composition_by_parts(sequence, k):
    part_length = len(sequence) // k

    compositions_by_parts = [None] * k
    for i in range(0, k):
        l = i * part_length
        r = l + part_length - 1
        compositions_by_parts[i] = sequence_composition(sequence[l:(r + 1)])

    return compositions_by_parts


# endregion

# region Analyze Gene nucleotide composition(stranded) by fragments

# gene structure: (UTR'5-CDS1)=EXON1 - INTRON - CDS2=EXON2 - ... - CDS(n-1)=EXON(n-1) - INTRON - (CDSn-UTR'3)=EXONn
# subregions of interest: UTR'5_Procession, UTR'5, inner CDSs, inner Introns, UTR'3, UTR'3_Procession
# for each subregion calculate (C_count,G_count,A_count,T_count)


def __get_regional_merged_sequences_from_gene(chr_id, gene, procession_length):
    # for NCBI database
    fragments = get_fragments_on_gene(chr_id, gene, force_sorted=True)

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
                utr5_sequence += retrieve_feature_sequence(chr_id, fragment)
            else:
                utr5_sequence = retrieve_feature_sequence(chr_id, fragment) + utr5_sequence
        if fragment.featuretype == 'three_prime_UTR':
            if gene.strand == '+':
                utr3_sequence += retrieve_feature_sequence(chr_id, fragment)
            else:
                utr3_sequence = retrieve_feature_sequence(chr_id, fragment) + utr3_sequence

        if fragment.featuretype == 'CDS':
            if gene.strand == '+':
                cds_sequence += retrieve_feature_sequence(chr_id, fragment)
            else:
                cds_sequence = retrieve_feature_sequence(chr_id, fragment) + cds_sequence

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
            intron_sequence += retrieve_interval_sequence(chr_id, last_exon.end + 1, fragment.start - 1,
                                                          fragment.strand)
        else:
            intron_sequence = retrieve_interval_sequence(chr_id, last_exon.end + 1, fragment.start - 1,
                                                         fragment.strand) + intron_sequence
        introns_len += fragment.start - last_exon.end - 1
        last_exon = fragment

    if left_fragment.strand == '+':
        utr5_procession_sequence = retrieve_interval_sequence(chr_id, left_fragment.start - procession_length,
                                                              left_fragment.start - 1, '+')

        utr3_procession_sequence = retrieve_interval_sequence(chr_id, right_fragment.end + 1,
                                                              right_fragment.end + procession_length, '+')

    else:
        utr3_procession_sequence = retrieve_interval_sequence(chr_id, left_fragment.start - procession_length,
                                                              left_fragment.start - 1, '-')

        utr5_procession_sequence = retrieve_interval_sequence(chr_id, right_fragment.end + 1,
                                                              right_fragment.end + procession_length, '-')

    return [utr5_procession_sequence, utr5_sequence, cds_sequence, intron_sequence, utr3_sequence,
            utr3_procession_sequence]


# output: occurrences[part][base] = [k][4]
def analyze_gene_occurrences(chr_id, gene: Feature, procession_length):
    regional_sequences = __get_regional_merged_sequences_from_gene(chr_id, gene, procession_length)

    utr5_procession_sequence = regional_sequences[0]
    utr5_sequence = regional_sequences[1]
    cds_sequence = regional_sequences[2]
    intron_sequence = regional_sequences[3]
    utr3_sequence = regional_sequences[4]
    utr3_procession_sequence = regional_sequences[5]

    occurrences_UTR5_procession = sequence_composition(utr5_procession_sequence)
    occurrences_UTR5 = sequence_composition(utr5_sequence)
    occurrences_CDS = sequence_composition(cds_sequence)
    occurrences_Introns = sequence_composition(intron_sequence)
    occurrences_UTR3 = sequence_composition(utr3_sequence)
    occurrences_UTR3_procession = sequence_composition(utr3_procession_sequence)

    return [occurrences_UTR5_procession, occurrences_UTR5, occurrences_CDS, occurrences_Introns, occurrences_UTR3,
            occurrences_UTR3_procession]


# output: occurrences[region][part][base] = [6][k][4]
def analyze_gene_occurrences_by_parts(chr_id, gene: Feature, k, procession_length):
    regional_sequences = __get_regional_merged_sequences_from_gene(chr_id, gene, procession_length)

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
    occurrences = [sequence_composition_by_parts(utr5_procession_sequence, k),
                   sequence_composition_by_parts(utr5_sequence, k), sequence_composition_by_parts(cds_sequence, k),
                   sequence_composition_by_parts(intron_sequence, k), sequence_composition_by_parts(utr3_sequence, k),
                   sequence_composition_by_parts(utr3_procession_sequence, k)]
    return occurrences

# endregion
