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
sequence_files_prefix = './genome_data/Ensembl/chromosomes/'

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


def max_fragments_overlap_length(fragments_A: list[Feature], fragments_B: list[Feature]) -> int:
    if len(fragments_A) == 0 or len(fragments_B) == 0: return False

    overlaps = []
    max_overlap = 0
    for fragment_a in fragments_A:
        for fragment_b in fragments_B:
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
        print("Fuck You!")
    if len(feature.attributes['Parent']) != 1:
        print("Fuck You!")

    # parent_gene_id for mRNA and parent_transcript_id for CDS/UTR5'/UTR3'/Exon
    parent_id = feature.attributes['Parent'][0]

    if not dict_to_fill.__contains__(parent_id):
        dict_to_fill[parent_id] = []

    dict_to_fill[parent_id].append(feature)


def preprocess_annotation_for_chr(chr_id, use_only_cds=False):
    if not exists(get_annotation_db_path(chr_id)):
        gffutils.create_db(get_annotation_file_path(chr_id), dbfn=get_annotation_db_path(chr_id), force=True,
                           keep_order=True,
                           merge_strategy='merge', sort_attribute_values=False)

    features_db = gffutils.FeatureDB(get_annotation_db_path(chr_id), keep_order=True)
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
def retrieve_sequence(chr_id) -> SeqRecord:
    if __sequences[chr_id] is None:
        sequence_path = get_sequence_file_path(chr_id)
        __sequences[chr_id] = next(SeqIO.parse(sequence_path, 'fasta'))

    return __sequences[chr_id]


def retrieve_feature_sequence(chr_id, feature: Feature) -> SeqRecord:
    seq_record = retrieve_sequence(chr_id)
    return seq_record[feature.start - 1:feature.end]


def retrieve_interval_sequence(chr_id, start, end) -> SeqRecord:
    seq_record = retrieve_sequence(chr_id)
    return seq_record[max(0, start - 1):min(end, len(seq_record))]

    # endregion

    # region Get transcript of a gene


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
                print("Fuck you!")
            if fragment.attributes['Parent'][0] == mRNA_transcript.id:
                fragments.append(fragment)

    return fragments


def find_best_gene_transcript(chr_id, gene_feature: Feature) -> Feature:
    if not __gene_transcripts.__contains__(gene_feature.id):
        # there is no valid mRNA transcript fot this specific gene
        return None

    mRNA_transcripts = __gene_transcripts[gene_feature.id]

    best_mRNA_transcript = None
    mRNA_variants = 0

    for transcript in mRNA_transcripts:
        if not transcript.attributes.__contains__('Parent'):
            print('Fuck You!')
        if not len(transcript.attributes) != 1:
            print('Fuck You!')
        if transcript.attributes['Parent'][0] == gene_feature.id:
            mRNA_variants = mRNA_variants + 1
            best_mRNA_transcript = better_annotated_transcript(best_mRNA_transcript, transcript)

    # print("mRNA variants for Gene [" + gene_feature.id + "]: " + str(mRNA_variants))
    # print("well annotated mRNA: " + best_mRNA_transcript.id)

    return best_mRNA_transcript


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


# retrieves nucleotide composition of sequence of given fragment/feature
# output: (C_count,G_count,A_count,T_count)
def fragment_composition(chr_id, fragment: Feature):
    sequence = retrieve_feature_sequence(chr_id, fragment)
    if fragment.strand == '-': sequence = sequence.reverse_complement()

    return sequence_composition(sequence)


# retrieves nucleotide composition of sequence of given interval (start,end,strand)
# output: (C_count,G_count,A_count,T_count)
def interval_composition(chr_id, start, end, strand: str):
    sequence = retrieve_interval_sequence(chr_id, start, end)
    if strand == '-': sequence = sequence.reverse_complement()

    return sequence_composition(sequence)


# endregion

# region Analyze Gene nucleotide composition(stranded) by fragments

# gene structure: (UTR'5-CDS1)=EXON1 - INTRON - CDS2=EXON2 - ... - CDS(n-1)=EXON(n-1) - INTRON - (CDSn-UTR'3)=EXONn
# subregions of interest: UTR'5_Procession, UTR'5, inner CDSs, inner Introns, UTR'3, UTR'3_Procession
# for each subregion calculate (C_count,G_count,A_count,T_count)

procession_seq_len = 2000


def analyze_gene_composition(chr_id, gene: Feature):
    fragments = get_fragments_on_gene(chr_id, gene)

    stats_UTR5_procession = [0, 0, 0, 0]
    stats_UTR5 = [0, 0, 0, 0]
    stats_CDS = [0, 0, 0, 0]
    stats_Introns = [0, 0, 0, 0]
    stats_UTR3 = [0, 0, 0, 0]
    stats_UTR3_procession = [0, 0, 0, 0]

    if len(fragments) == 0: return [stats_UTR5, stats_CDS, stats_Introns, stats_UTR3]

    for fragment in fragments:
        if fragment.featuretype == 'five_prime_UTR':
            stats = fragment_composition(chr_id, fragment)
            stats_UTR5 = [stats_UTR5[x] + stats[x] for x in range(len(stats))]
        if fragment.featuretype == 'three_prime_UTR':
            stats = fragment_composition(chr_id, fragment)
            stats_UTR3 = [stats_UTR3[x] + stats[x] for x in range(len(stats))]
        if fragment.featuretype == 'CDS':
            stats = fragment_composition(chr_id, fragment)
            stats_CDS = [stats_CDS[x] + stats[x] for x in range(len(stats))]

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

        # fragment.strand must be same for all fragments/exons
        stats = interval_composition(chr_id, last_exon.end + 1, fragment.start - 1, fragment.strand)
        stats_Introns = [stats_Introns[x] + stats[x] for x in range(len(stats))]

        introns_len += fragment.start - last_exon.end - 1
        last_exon = fragment

    if left_fragment.strand == '+':
        stats = interval_composition(chr_id, left_fragment.start - procession_seq_len, left_fragment.start - 1, '+')
        stats_UTR5_procession = [stats_UTR5_procession[x] + stats[x] for x in range(len(stats))]

        stats = interval_composition(chr_id, right_fragment.end + 1, right_fragment.end + procession_seq_len, '+')
        stats_UTR3_procession = [stats_UTR3_procession[x] + stats[x] for x in range(len(stats))]
    else:
        stats = interval_composition(chr_id, left_fragment.start - procession_seq_len, left_fragment.start - 1, '-')
        stats_UTR3_procession = [stats_UTR3_procession[x] + stats[x] for x in range(len(stats))]

        stats = interval_composition(chr_id, right_fragment.end + 1, right_fragment.end + procession_seq_len, '-')
        stats_UTR5_procession = [stats_UTR5_procession[x] + stats[x] for x in range(len(stats))]

    return [stats_UTR5_procession, stats_UTR5, stats_CDS, stats_Introns, stats_UTR3, stats_UTR3_procession]
# endregion
