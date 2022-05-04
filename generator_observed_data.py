# Author: "Lasha Bukhnikashvili"
#
# Description:
#   Calculates overlapping gene pairs by
#   specific parts: CDS (coding sequence), EXON or whole gene
#
# Usage:
#   generator_observed_data.py <species> <annotation> <strand_similarity> <ORF_similarity>
#
# Params (possible) to run:
#   annotation: NCBI / Ensembl
#   fragment_type: CDS / exon / gene
#   strand_similarity: same_stranded / diff_stranded / both_stranded
#   ORF_similarity: same_ORF / diff_ORF / both_ORF
#
# Example:
#   python ./generator_observed_data.py Drosophila_melanogaster Ensembl diff_stranded both_ORF
#   Finds all overlapping genes by CDS fragments, which are on different strands
#   and which have same ORF or different ORF
#
# Output:
#   creates .txt file in ./results folder

from analyzer_worker import *
from genome_worker import *
from genome_worker_enums import *
import time
import sys

start_time = time.time()

# region clusters (by size) counting method

cluster_indexes = []
cluster_sizes = []
clusters_count_by_size = []


def initialize_clusters(cnt):
    global cluster_indexes, cluster_sizes, clusters_count_by_size
    cluster_indexes = [0] * cnt
    cluster_sizes = [0] * cnt
    clusters_count_by_size = [0] * cnt
    for i in range(0, cnt):
        cluster_indexes[i] = i
        cluster_sizes[i] = 1
        clusters_count_by_size[1] += 1


def decrease_cluster_size(cluster_index):
    # x-sized clusters count will decrease, if x was size for cluster_index
    clusters_count_by_size[cluster_sizes[cluster_index]] -= 1
    # cluster size for cluster_index reduce
    cluster_sizes[cluster_index] -= 1
    # x-sized clusters count will increase, if x is now size for cluster_index
    clusters_count_by_size[cluster_sizes[cluster_index]] += 1


def increase_cluster_size(cluster_index):
    # x-sized clusters count will decrease, if x was size for cluster_index
    clusters_count_by_size[cluster_sizes[cluster_index]] -= 1
    # cluster size for cluster_index increase
    cluster_sizes[cluster_index] += 1
    # x-sized clusters count will increase, if x is now size for cluster_index
    clusters_count_by_size[cluster_sizes[cluster_index]] += 1


def union_clusters(x, y, cnt):
    global cluster_indexes, cluster_sizes, clusters_count_by_size
    # let cluster-x merge to cluster-y, so new cluster index for all genes will be y
    old_cluster_index = cluster_indexes[x]
    new_cluster_index = cluster_indexes[y]
    for element_index in range(0, cnt):
        if cluster_indexes[element_index] == old_cluster_index:
            cluster_indexes[element_index] = new_cluster_index  # change cluster index for this gene, with new
            decrease_cluster_size(old_cluster_index)  # decrease cluster size of old_cluster_index
            increase_cluster_size(new_cluster_index)  # increase cluster size of new_cluster_index


# gives clusters by size, where size >1
def get_clusters_by_size(cnt):
    result = {}
    for size in range(2, cnt):
        if clusters_count_by_size[size] > 0:
            result[size] = clusters_count_by_size[size]
    return result


# endregion

# region variables and input management
assert len(sys.argv) == 5
assert sys.argv[2] == 'NCBI' or sys.argv[2] == 'Ensembl'
assert sys.argv[3] == 'same_stranded' or sys.argv[3] == 'diff_stranded' or sys.argv[3] == 'both_stranded'
assert sys.argv[4] == 'same_ORF' or sys.argv[4] == 'diff_ORF' or sys.argv[4] == 'both_ORF'

species = SPECIES.from_string(sys.argv[1])
annotation = ANNOTATIONS.NCBI if sys.argv[2] == 'NCBI' else ANNOTATIONS.ENSEMBL
strand_similarity, ORF_similarity = sys.argv[3], sys.argv[4]

overlapped_genes = 0
overlapped_genes_by_cluster_size = {}
cds_overlaps = []

# endregion

# region observed data ( OGs by CDS) collector

# Load necessary annotation data
genome = GenomeWorker(species, annotation, ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_CDS, SEQUENCE_LOAD.LOAD)
stats = AnalyzerData()

total_max_overlapped_length = 0

for chr_id in range(1, genome.chromosomes_count() + 1):
    genes_cnt = genome.genes_count_on_chr(chr_id)

    # initialise clusters, for later cluster counting
    initialize_clusters(genes_cnt)

    for i in range(0, genes_cnt):
        for j in range(i + 1, genes_cnt):
            gene_A = genome.gene_by_ind(chr_id, i)
            gene_B = genome.gene_by_ind(chr_id, j)

            # check gene pairs by strand_similarity param
            if strand_similarity == 'same_stranded' and gene_A.strand != gene_B.strand: continue
            if strand_similarity == 'diff_stranded' and gene_A.strand == gene_B.strand: continue

            # if genes don't overlap, there is no point in searching overlaps by fragments (CDS in this case)
            if not genome.are_segments_overlapped((gene_A.start, gene_A.end), (gene_B.start, gene_B.end)): continue

            # if criteria != TRANSCRIPT_CRITERIA.NONE:
            #     mRNA_transcript_a = self.find_gene_transcript_by_criteria(gene_A, criteria)
            #     mRNA_transcript_b = self.find_gene_transcript_by_criteria(gene_B, criteria)
            #     return self.__get_fragments_overlap_between_transcripts(mRNA_transcript_a, mRNA_transcript_b,
            #                                                             ORF_similarity)[0]
            max_overlapped_segments = []
            max_overlapped_length = 0

            overlapped_segments = []
            transcripts_A = genome.get_gene_all_transcripts(gene_A)
            transcripts_B = genome.get_gene_all_transcripts(gene_B)
            for transcript_a in transcripts_A:
                for transcript_b in transcripts_B:
                    segments = genome.get_fragments_overlap_between_transcripts(transcript_a, transcript_b,
                                                                                ORF_similarity)
                    overlapped_segments += segments
                    length = 0
                    for segment in segments:
                        assert segment[1] >= segment[0]
                        length += segment[1] - segment[0]+1
                    if length > max_overlapped_length:
                        max_overlapped_length = length
                        max_overlapped_segments = segments
            overlapped_segments = max_overlapped_segments

            # maybe genes overlapped, but none of the fragments/transcripts (cds) overlapped
            if len(overlapped_segments) > 0:
                union_clusters(i, j, genes_cnt)
                cds_overlaps.append((max_overlapped_length, "", chr_id, gene_A, gene_B))
                total_max_overlapped_length += max_overlapped_length

            for interval in overlapped_segments:
                assert interval[1] >= interval[0]
                if gene_A.strand == gene_B.strand:
                    seq = genome.retrieve_interval_sequence(chr_id, interval[0], interval[1], gene_A.strand)
                    stats.analyze_sequence_stats(seq, interval[2])
                    stats.analyze_sequence_stats(seq, interval[3])
                else:
                    sense_seq = genome.retrieve_interval_sequence(chr_id, interval[0], interval[1], '+')
                    anti_sense_seq = genome.retrieve_interval_sequence(chr_id, interval[0], interval[1], '-')
                    stats.analyze_sequence_stats(sense_seq, interval[2])
                    stats.analyze_sequence_stats(anti_sense_seq, interval[3])

    # merge clusters size data from this chromosome, to overall data
    clusters_by_size_on_chr = get_clusters_by_size(genes_cnt)  # it only gives clusters with size >1

    for cluster_size in clusters_by_size_on_chr.keys():
        if not overlapped_genes_by_cluster_size.__contains__(cluster_size):
            overlapped_genes_by_cluster_size[cluster_size] = 0
        overlapped_genes_by_cluster_size[cluster_size] += clusters_by_size_on_chr[cluster_size]

        # increase genes count which is any of the cluster
        overlapped_genes += cluster_size * clusters_by_size_on_chr[cluster_size]

# endregion

# region export generated data

# sort overlapped gene pairs by max overlapped length between their fragments
cds_overlaps.sort(key=lambda tup: tup[0], reverse=True)

# example path './results/overlaps - CDS [diff_stranded, SAME_ORF] (NCBI) .txt'
output_data_path = f'generated_data/observed/OGs by CDS [{species.short_name()}, {annotation.short_name()}, {strand_similarity}, {ORF_similarity}].txt'

with open(output_data_path, 'w') as f:
    f.write(f'Total overlapped (by CDS) genes: {overlapped_genes} =')
    is_first = True
    keys = list(overlapped_genes_by_cluster_size.keys())
    keys.sort()
    for i in range(0, len(keys)):
        if i > 0: f.write(' +')
        f.write(f' [{keys[i]}]x{overlapped_genes_by_cluster_size[keys[i]]}')
    if len(keys) == 0: f.write(' 0')
    f.write('\n')
    f.write(f'Total overlapped (by CDS) sequence length: {total_max_overlapped_length}\n\n')
    f.write(f'-----------------------------Overlapped Sequence Composition-----------------------------\n\n')
    f.write(stats.get_formatted_results(indent_tabs=1, use_percents=True))
    f.write(f'\n------------------------------------Overlapped Gene------------------------------------\n\n')

    for item in cds_overlaps:
        f.write(f'Total overlap length: {item[0]}{item[1]}  |  chr: {item[2]} |  genes: [{item[3].id};{item[4].id}]\n')
        gene1_name = GenomeWorker.get_gene_name(item[3])
        gene1_desc = GenomeWorker.get_gene_description(item[3])

        gene2_name = GenomeWorker.get_gene_name(item[4])
        gene2_desc = GenomeWorker.get_gene_description(item[4])

        f.write(f'\tGene 1: {gene1_name}; {gene1_desc}\n')
        f.write(f'\tGene 2: {gene2_name}; {gene2_desc}\n')

        f.write('\n')

# endregion

print("--- script was running for %s seconds. ---" % (time.time() - start_time))
