from Bio.SeqRecord import SeqRecord

import genome_tools as genome
import time

start_time = time.time()

coding_overlapped_genes = []

max_overlap_sums = 0
total_overlap_sum = 0

total_overlapping_composition = [0, 0, 0, 0]

for chr_id in range(22, 23):
    genome.preprocess_annotation_for_chr(chr_id, True)
    genes_cnt = genome.genes_count_on_chr(chr_id)

    cds=""
    for i in range(0,genome.genes_count_on_chr(chr_id)):
        gene=genome.gene_by_ind(chr_id,i)
        if gene.id=="gene:ENSG00000099901":
            fragments=genome.get_fragments_on_gene(chr_id,gene)
            for fragment in fragments:
                if fragment.featuretype!="CDS": continue
                print("")
                gene_record=genome.retrieve_feature_sequence(chr_id,fragment)
                cds+=gene_record

                gene_record.translate()
                print(gene_record.translate().format("fasta"))
    x=SeqRecord();
#
#     genes = 0
#     maxOverlap = -1
#     desired_pair = (-1, -1)
#     gene_fragments = []
#     for i in range(0, genes_cnt):
#         gene_fragments.append([])
#
#     for i in range(0, genes_cnt):
#         gene_fragments[i] = genome.get_fragments_on_gene(chr_id, genome.gene_by_ind(chr_id, i))
#
#     for i in range(0, genes_cnt):
#         for j in range(i + 1, genes_cnt):
#             if not genome.are_genes_overlapped(genome.gene_by_ind(chr_id, i), genome.gene_by_ind(chr_id, j)): continue
#             if genome.gene_by_ind(chr_id, i).strand == genome.gene_by_ind(chr_id, j).strand: continue
#             max_overlap, overlap_intervals = genome.max_fragments_overlap_length(gene_fragments[i], gene_fragments[j])
#
#             total_overlap = 0
#             if len(overlap_intervals) > 1:
#                 print(genome.gene_by_ind(chr_id, i))
#
#             for interval in overlap_intervals:
#                 assert interval[1] >= interval[0]
#                 total_overlap += interval[1] - interval[0] + 1
#                 comp = genome.interval_composition(chr_id, interval[0], interval[1], '-')
#                 total_overlapping_composition = [total_overlapping_composition[x] + comp[x] for x in range(0, 4)]
#
#             if max_overlap > 0:
#                 max_overlap_sums += max_overlap
#                 total_overlap_sum += total_overlap
#                 coding_overlapped_genes.append(
#                     (chr_id, max_overlap, total_overlap, genome.gene_by_ind(chr_id, i), genome.gene_by_ind(chr_id, j)))
#
#     print("Overlaps Found: " + str(len(coding_overlapped_genes)))
#
# sorted_by_second = sorted(coding_overlapped_genes, key=lambda tup: tup[1])
#
# with open('cds_overlapped_genes_test.txt', 'w') as f:
#     for item in sorted_by_second:
#         data = "chr-" + str(item[0]) + "  overlap_max-" + str(item[1]) + "  overlap_total-" + str(
#             item[2]) + "  genes-[" + item[3].id + ";" + item[
#                    4].id + "]";
#
#         # data += "\n"
#         # if item[2].attributes.__contains__('Name') and item[2].attributes.__contains__('description'):
#         #     data += item[2].attributes['Name'][0] + "; " + item[2].attributes['description'][0]
#         # data += "\n"
#         # if item[3].attributes.__contains__('Name') and item[3].attributes.__contains__('description'):
#         #     data += item[3].attributes['Name'][0] + "; " + item[3].attributes['description'][0]
#         # data += "\n"
#         #
#         # data += "\n"
#
#         f.write("%s\n" % data)
#     data2 = "max_overlap_sums: " + str(max_overlap_sums) + " total_overlap_sums: " + str(total_overlap_sum)
#     f.write("%s\n" % data2)
#     data3 = "Total Overlap Composition: C-" + str(total_overlapping_composition[0]) + " G-" + str(
#         total_overlapping_composition[1]) + " A-" + str(total_overlapping_composition[2]) + " T-" + str(
#         total_overlapping_composition[3])
#     f.write("%s\n" % data3)
#     print("--- %s seconds ---" % (time.time() - start_time))
