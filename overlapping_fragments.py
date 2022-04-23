# Author: "Lasha Bukhnikashvili"
#
# Description:
#   Calculates overlapping gene pairs by
#   specific parts: CDS (coding sequence), EXON or whole gene
#
# Usage:
#   overlapping_fragments.py <annotation> <fragment_type> <strand_similarity> <ORF_similarity> <output_type>
#
# Params (possible) to run:
#   annotation: NCBI / Ensembl
#   fragment_type: CDS / exon / gene
#   strand_similarity: same_stranded / diff_stranded / both_stranded
#   ORF_similarity: same_ORF / diff_ORF / both_ORF
#   output_type: print_gene_pairs / print_gene_pairs_with_descriptions
#
# Example:
#   python overlapping_fragments.py NCBI CDS diff_stranded both_framed print_gene_pairs_with_descriptions
#   Finds all overlapping genes by CDS fragments, which are on different strands
#   and which have same ORF or different ORF
#
# Output:
#   creates .txt file in ./results folder

import genome_worker as genome
from genome_worker import ANNOTATIONS
from genome_worker import ANNOTATION_LOAD
from genome_worker import SEQUENCE_LOAD

import time
import sys

start_time = time.time()

assert len(sys.argv) == 6

annotation = ANNOTATIONS.NCBI if sys.argv[1] == 'NCBI' else ANNOTATIONS.ENSEMBL
annotation_load_type = ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_CDS if sys.argv[2] == 'CDS' \
    else ANNOTATION_LOAD.GENES if sys.argv[2] == 'gene' else ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_FRAGMENTS
strand_similarity, ORF_similarity, output_type = sys.argv[3], sys.argv[4], sys.argv[5]

coding_overlapped_genes = []
overlapped_gene_pairs = 0
total_overlap_sum = 0
overlapped_merged_sequence = ""

# Load necessary annotation data
genome.preprocess_annotation(annotation, annotation_load_type, SEQUENCE_LOAD.LOAD)

for chr_id in range(1, genome.chromosomes_count() + 1):
    genes_cnt = genome.genes_count_on_chr(chr_id)

    genes = 0
    maxOverlap = -1
    desired_pair = (-1, -1)

    for i in range(0, genes_cnt):
        for j in range(i + 1, genes_cnt):
            gene_A = genome.gene_by_ind(chr_id, i)
            gene_B = genome.gene_by_ind(chr_id, j)

            # check gene pairs by strand_similarity param
            if strand_similarity == 'same_stranded' and gene_A.strand != gene_B.strand: continue
            if strand_similarity == 'diff_stranded' and gene_A.strand == gene_B.strand: continue

            # if genes don't overlap, there is no point for searching fragments overlap
            if not genome.are_genes_overlapped(gene_A, gene_B): continue

            # if we are searching genes overlap
            if annotation_load_type == ANNOTATION_LOAD.GENES:
                l, r = genome.get_genes_overlap(gene_A, gene_B)
                overlap_intervals = [(l, r)]
            else:  # if we are searching overlaps within fragments of gene pairs
                overlap_intervals = genome.get_fragments_overlap(gene_A, gene_B, ORF_similarity)

            total_overlap = 0
            max_overlap = 0
            for interval in overlap_intervals:
                assert interval[1] >= interval[0]
                max_overlap = max(max_overlap, interval[1] - interval[0] + 1)
                total_overlap += interval[1] - interval[0] + 1
                overlapped_merged_sequence += genome.retrieve_interval_sequence(chr_id, interval[0], interval[1],
                                                                                gene_A.strand)

            # maybe genes overlapped, but none of the fragments (cds or exon) overlapped
            if max_overlap > 0:
                overlapped_gene_pairs += 1
                total_overlap_sum += total_overlap
                coding_overlapped_genes.append(
                    (max_overlap, total_overlap, chr_id, gene_A, gene_B))

    print("Total overlaps found: " + str(len(coding_overlapped_genes)))

# sort overlapped gene pairs by max overlapped length between their fragments
sorted_by_second = sorted(coding_overlapped_genes, key=lambda tup: tup[0], reverse=True)
total_overlapping_composition = genome.sequence_composition(overlapped_merged_sequence)

# output data
# example path './results/overlaps - CDS [diff_stranded, SAME_ORF] (NCBI) .txt'
output_data_path = './results/OGs' + ('' if sys.argv[5] == 'print_gene_pairs' else ' (+descriptions)') \
                   + ' - ' + sys.argv[2] + " [" + strand_similarity + ", " + ORF_similarity + "] (" + \
                   sys.argv[1] + ").txt"
with open(output_data_path, 'w') as f:
    for item in sorted_by_second:
        data = "Max overlap: " + str(item[0]) + "  |  Total overlap: " + str(item[1]) + "  |  chr: " + str(
            item[2]) + "  |  genes: [" + item[3].id + ";" + item[4].id + "]"

        # print genes description?
        if output_type == 'print_gene_pairs_with_descriptions':
            data += "\n     Gene 1 Description: "
            if item[3].attributes.__contains__('Name') and item[3].attributes.__contains__('description'):
                data += item[3].attributes['Name'][0] + "; " + item[3].attributes['description'][0]
            data += "\n     Gene 2 Description: "
            if item[4].attributes.__contains__('Name') and item[4].attributes.__contains__('description'):
                data += item[4].attributes['Name'][0] + "; " + item[4].attributes['description'][0]
            data += "\n"

            f.write("%s\n" % data)

    f.write("\n%s\n" % "Stats:")
    data2 = "   Total overlapped gene pairs: " + str(overlapped_gene_pairs)
    f.write("%s\n" % data2)
    data2 = "   Total overlapped length: " + str(total_overlap_sum)
    f.write("%s\n" % data2)
    data3 = "   Overlapped segments composition: C-" + str(total_overlapping_composition[0]) + " G-" + str(
        total_overlapping_composition[1]) + " A-" + str(total_overlapping_composition[2]) + " T-" + str(
        total_overlapping_composition[3])
    f.write("%s\n" % data3)

print("--- script was running for %s seconds. ---" % (time.time() - start_time))
