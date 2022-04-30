# Author: "Lasha Bukhnikashvili"
#
# Description:
#   Calculates overlapping gene pairs by
#   specific parts: CDS (coding sequence), EXON or whole gene
#
# Usage:
#   analyzer_OGs_by_CDS.py <species> <annotation> <strand_similarity> <ORF_similarity>
#
# Params (possible) to run:
#   annotation: NCBI / Ensembl
#   fragment_type: CDS / exon / gene
#   strand_similarity: same_stranded / diff_stranded / both_stranded
#   ORF_similarity: same_ORF / diff_ORF / both_ORF
#
# Example:
#   python ./analyzer_OGs_by_CDS.py Drosophila_melanogaster Ensembl diff_stranded both_ORF
#   Finds all overlapping genes by CDS fragments, which are on different strands
#   and which have same ORF or different ORF
#
# Output:
#   creates .txt file in ./results folder

from genome_worker import *
from analyzer_data import *

import time
import sys

start_time = time.time()

assert len(sys.argv) == 5
assert sys.argv[2] == 'NCBI' or sys.argv[2] == 'Ensembl'
assert sys.argv[3] == 'same_stranded' or sys.argv[3] == 'diff_stranded' or sys.argv[3] == 'both_stranded'
assert sys.argv[4] == 'same_ORF' or sys.argv[4] == 'diff_ORF' or sys.argv[4] == 'both_ORF'

species = SPECIES.from_string(sys.argv[1])
annotation = ANNOTATIONS.NCBI if sys.argv[2] == 'NCBI' else ANNOTATIONS.ENSEMBL
strand_similarity, ORF_similarity = sys.argv[3], sys.argv[4]

cds_overlaps = []
overlapped_genes = 0
overlapped_gene_pairs = 0
total_overlap_sum = 0

# Load necessary annotation data
genome = GenomeWorker(species, annotation, ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_CDS, SEQUENCE_LOAD.LOAD)
stats = AnalyzerData()

for chr_id in range(1, genome.chromosomes_count() + 1):
    genes_cnt = genome.genes_count_on_chr(chr_id)
    overlapped_genes_fix = [False] * genes_cnt

    for i in range(0, genes_cnt):
        for j in range(i + 1, genes_cnt):
            gene_A = genome.gene_by_ind(chr_id, i)
            gene_B = genome.gene_by_ind(chr_id, j)

            # check gene pairs by strand_similarity param
            if strand_similarity == 'same_stranded' and gene_A.strand != gene_B.strand: continue
            if strand_similarity == 'diff_stranded' and gene_A.strand == gene_B.strand: continue

            # if genes don't overlap, there is no point in searching overlaps by fragments (CDS in this case)
            if not genome.are_features_overlapped(gene_A, gene_B): continue

            #  we are searching overlaps within fragments of gene pairs
            overlap_intervals = genome.get_fragments_overlap_between_genes(gene_A, gene_B, ORF_similarity,
                                                                           TRANSCRIPT_CRITERIA.NONE)

            total_overlap = 0
            overlaps_description = ""
            for interval in overlap_intervals:
                assert interval[1] >= interval[0]

                if len(overlaps_description) == 0:
                    overlaps_description = f" = {interval[1] - interval[0] + 1}"
                else:
                    overlaps_description = f'{overlaps_description} + {interval[1] - interval[0] + 1}'

                total_overlap += interval[1] - interval[0] + 1

                if gene_A.strand == gene_B.strand:
                    seq = genome.retrieve_interval_sequence(chr_id, interval[0], interval[1], gene_A.strand)

                    stats.analyze_sequence_stats(seq, interval[2])
                    stats.analyze_sequence_stats(seq, interval[3])
                else:
                    sense_seq = genome.retrieve_interval_sequence(chr_id, interval[0], interval[1], '+')
                    anti_sense_seq = genome.retrieve_interval_sequence(chr_id, interval[0], interval[1], '-')

                    stats.analyze_sequence_stats(sense_seq, interval[2])
                    stats.analyze_sequence_stats(anti_sense_seq, interval[3])

            # maybe genes overlapped, but none of the fragments (cds or exon) overlapped
            if len(overlap_intervals) > 0:
                overlapped_genes += 1 if overlapped_genes_fix[i] is False else 0
                overlapped_genes += 1 if overlapped_genes_fix[j] is False else 0

                overlapped_genes_fix[i] = True
                overlapped_genes_fix[j] = True
                overlapped_gene_pairs += 1

                overlaps_description = overlaps_description if len(overlap_intervals) > 1 else ""

                total_overlap_sum += total_overlap

                cds_overlaps.append((total_overlap, overlaps_description, chr_id, gene_A, gene_B))

# sort overlapped gene pairs by max overlapped length between their fragments
sorted_by_second = sorted(cds_overlaps, key=lambda tup: tup[0], reverse=True)

# example path './results/overlaps - CDS [diff_stranded, SAME_ORF] (NCBI) .txt'
output_data_path = f'./results/OGs by CDS [{strand_similarity}, {ORF_similarity}, {species.short_name()}, {annotation.short_name()}].txt'

with open(output_data_path, 'w') as f:
    f.write(f'Total overlapped (by CDS) genes: {overlapped_genes}\n')
    f.write(f'Total overlapped (by CDS) gene pairs: {overlapped_gene_pairs}\n')
    f.write(f'Total overlapped (by CDS) sequence length: {total_overlap_sum}\n\n')
    f.write(f'-----------------------------Overlapped Sequence Composition-----------------------------\n\n')
    f.write(stats.get_formatted_results(indent_tabs=1, use_percents=True))
    f.write(f'\n------------------------------------Overlapped Gene------------------------------------\n\n')

    for item in sorted_by_second:
        f.write(f'Total overlap length: {item[0]}{item[1]}  |  chr: {item[2]} |  genes: [{item[3].id};{item[4].id}]\n')
        gene1_name = GenomeWorker.get_gene_name(item[3])
        gene1_desc = GenomeWorker.get_gene_description(item[3])

        gene2_name = GenomeWorker.get_gene_name(item[4])
        gene2_desc = GenomeWorker.get_gene_description(item[4])

        f.write(f'\tGene 1: {gene1_name}; {gene1_desc}\n')
        f.write(f'\tGene 2: {gene2_name}; {gene2_desc}\n')

        f.write('\n')

print("--- script was running for %s seconds. ---" % (time.time() - start_time))
