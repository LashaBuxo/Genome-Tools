# Lasha
# Specific tools for working with Ensembl annotation and with sequence data
from os.path import exists
import gffutils
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from gffutils import FeatureDB, Feature

# region Variables and Retrievers

#  22 autosome + 2 sex + 1 mitochondrion chromosome
number_of_chromosomes = 25

# genome files: annotations & DNA sequences
annotation_files_prefix = './genome_data/Ensembl/annotations/'
annotation_db_files_prefix = './databases/Ensembl/'
sequence_files_prefix = './genome_data/Ensembl/sequences/'

# gene features clustered by chromosome, they are on.
# Example: genes_on_chr[22] = ['gene' features on chromosome 22]
__genes_on_chr = [None] * (number_of_chromosomes + 1)

# {chr_id,DNA sequence of chromosome with id chr_id}
__sequences = [None] * (number_of_chromosomes + 1)

#   {gene_id, [list of mRNA features whose parent is gene_id]}
__gene_transcripts = {}

#   {mRNA_transcript_id, [list of 'CDS' and/or 'UTR(3'/5')' and/or 'exons' features whose parent is mRNA_transcript_id]}
__transcript_fragments = {}


def genes_count_on_chr(chr_id):
    assert chr_id <= number_of_chromosomes
    return len(__genes_on_chr[chr_id])


def gene_by_ind(chr_id, index) -> Feature:
    assert chr_id <= number_of_chromosomes
    return __genes_on_chr[chr_id][index]


def get_annotation_file_path(chr_id):
    return annotation_files_prefix + 'chr' + str(chr_id) + ".gff3"


def get_annotation_db_path(chr_id):
    return annotation_db_files_prefix + 'chr' + str(chr_id) + '.db'


def get_sequence_file_path(chr_id):
    return sequence_files_prefix + 'chr' + str(chr_id) + ".fa"


# endregion

# region Additional Methods

def better_annotated_transcript(A: Feature, B: Feature) -> Feature:
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
    return False


def is_same_frame_overlap(fragment_a: Feature, fragment_b: Feature) -> bool:
    if fragment_a.strand != fragment_b.strand: return False
    l = 0
    r = 0
    if fragment_a.strand == '+':
        l = fragment_a.start + int(fragment_a.frame)
        r = fragment_b.start + int(fragment_b.frame)
    else:
        l = fragment_a.end - int(fragment_a.frame)
        r = fragment_b.end - int(fragment_b.frame)
    return l % 3 == r % 3


def max_fragments_overlap_length(fragments_A: list[Feature], fragments_B: list[Feature],
                                 exclude_same_stranded_and_same_frame_overlaps) -> int:
    if len(fragments_A) == 0 or len(fragments_B) == 0: return False

    same_stranded = fragments_A[0].strand == fragments_B[0].strand

    overlaps = []
    max_overlap = 0
    for fragment_a in fragments_A:
        for fragment_b in fragments_B:
            if exclude_same_stranded_and_same_frame_overlaps and is_same_frame_overlap(fragment_a, fragment_b): continue
            overlap_len = 0
            overlap = (0, 0)
            if fragment_b.end < fragment_a.start or fragment_b.start > fragment_a.end: continue
            if fragment_b.end <= fragment_a.end:
                if fragment_b.start >= fragment_a.start:
                    overlap_len = fragment_b.end - fragment_b.start + 1
                    overlap = (fragment_b.start, fragment_b.end)
                else:
                    overlap_len = fragment_b.end - fragment_a.start + 1
                    overlap = (fragment_a.start, fragment_b.end)
            else:
                if fragment_b.start <= fragment_a.start:
                    overlap_len = fragment_a.end - fragment_a.start + 1
                    overlap = (fragment_a.start, fragment_a.end)
                else:
                    overlap_len = fragment_a.end - fragment_b.start + 1
                    overlap = (fragment_b.start, fragment_a.end)
            max_overlap = max(max_overlap, overlap_len)
            overlaps.append(overlap)
    return max_overlap, overlaps


# endregion

# region Annotation Preprocessing methods

# if database does not exists for specific chromosome, then
# builds database from annotation file and fills/loads all
# necessary features from the database.
#
# if database exists, then fills/loads all necessary data structures

def preprocess_feature_by_type(dict_to_fill, feature: Feature, feature_type):
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


def preprocess_annotation_for_chr(chr_id, use_only_cds=False):
    if not exists(get_annotation_db_path(chr_id)):
        gffutils.create_db(get_annotation_file_path(chr_id), dbfn=get_annotation_db_path(chr_id), force=False,
                           keep_order=False,
                           merge_strategy='merge', sort_attribute_values=False)

    features_db = gffutils.FeatureDB(get_annotation_db_path(chr_id), keep_order=False)

    features_generator = features_db.all_features()
    features = list(features_generator)

    for feature in features:

        if feature.featuretype == 'gene':
            # choose genes who has protein_coding attribute
            if feature.attributes.__contains__('biotype') and feature.attributes['biotype'][0] == 'protein_coding':
                if __genes_on_chr[chr_id] is None:
                    __genes_on_chr[chr_id] = []
                __genes_on_chr[chr_id].append(feature)

        # load best transcripts and link them to genes in dictionaries

        preprocess_feature_by_type(__gene_transcripts, feature, 'mRNA')
        # preprocess_feature_by_type(__gene_transcripts, feature, 'unconfirmed_transcript')
        # preprocess_feature_by_type(__gene_transcripts, feature, 'V_gene_segment')
        # preprocess_feature_by_type(__gene_transcripts, feature, 'J_gene_segment')
        # preprocess_feature_by_type(__gene_transcripts, feature, 'C_gene_segment')

        # load best fragments and link them to transcripts in dictionaries

        preprocess_feature_by_type(__transcript_fragments, feature, 'CDS')
        if not use_only_cds:
            preprocess_feature_by_type(__transcript_fragments, feature, 'exon')
            preprocess_feature_by_type(__transcript_fragments, feature, 'three_prime_UTR')
            preprocess_feature_by_type(__transcript_fragments, feature, 'five_prime_UTR')


# endregion

# region Get sequence of a Gene

# if sequence for specific chromosome is not loaded then loads it.
# if its already loaded retrieves it.
def retrieve_sequence_record(chr_id) -> SeqRecord:
    if __sequences[chr_id] is None:
        sequence_path = get_sequence_file_path(chr_id)
        __sequences[chr_id] = next(SeqIO.parse(sequence_path, 'fasta'))

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


def find_best_gene_transcript(chr_id, gene_feature: Feature) -> Feature:
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
            best_mRNA_transcript = better_annotated_transcript(best_mRNA_transcript, transcript)

    # print("mRNA variants for Gene [" + gene_feature.id + "]: " + str(mRNA_variants))
    # print("well annotated mRNA: " + best_mRNA_transcript.id)

    return best_mRNA_transcript


# endregion

# region Get fragments of a gene


def find_transcript_fragments(chr_id, mRNA_transcript: Feature) -> list[Feature]:
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


def get_fragments_on_gene(chr_id, gene_feature: Feature):
    # mRNA_transcript (start,end) always would fit in parent's (start,end) interval
    mRNA_transcript = find_best_gene_transcript(chr_id, gene_feature)
    if mRNA_transcript is None:
        # it seems there is no valid mRNA for the gene
        # print('no valid mRNA and exons for a gene')
        return []

    fragments = find_transcript_fragments(chr_id, mRNA_transcript)

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

procession_seq_len = 1000


def get_regional_merged_sequences_from_gene(chr_id, gene):
    fragments = get_fragments_on_gene(chr_id, gene)

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
        if fragment.featuretype == 'five_prime_UTR': utr5_sequence = retrieve_feature_sequence(chr_id, fragment)
        if fragment.featuretype == 'three_prime_UTR': utr3_sequence = retrieve_feature_sequence(chr_id, fragment)
        if fragment.featuretype == 'CDS': cds_sequence += retrieve_feature_sequence(chr_id, fragment)

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

        intron_sequence += retrieve_interval_sequence(chr_id, last_exon.end + 1, fragment.start - 1, fragment.strand)
        introns_len += fragment.start - last_exon.end - 1
        last_exon = fragment

    if left_fragment.strand == '+':
        utr5_procession_sequence = retrieve_interval_sequence(chr_id, left_fragment.start - procession_seq_len,
                                                              left_fragment.start - 1, '+')

        utr3_procession_sequence = retrieve_interval_sequence(chr_id, right_fragment.end + 1,
                                                              right_fragment.end + procession_seq_len, '+')

    else:
        utr3_procession_sequence = retrieve_interval_sequence(chr_id, left_fragment.start - procession_seq_len,
                                                              left_fragment.start - 1, '-')

        utr5_procession_sequence = retrieve_interval_sequence(chr_id, right_fragment.end + 1,
                                                              right_fragment.end + procession_seq_len, '-')

    return [utr5_procession_sequence, utr5_sequence, cds_sequence, intron_sequence, utr3_sequence,
            utr3_procession_sequence]


# output: occurrences[part][base] = [k][4]
def analyze_gene_occurrences(chr_id, gene: Feature):
    regional_sequences = get_regional_merged_sequences_from_gene(chr_id, gene)

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
def analyze_gene_occurrences_by_parts(chr_id, gene: Feature, k):
    regional_sequences = get_regional_merged_sequences_from_gene(chr_id, gene)

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
