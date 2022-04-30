# python ./analyzer_CDS_control_stats.py Drosophila_melanogaster Ensembl

from genome_worker import *
from analyzer_data import *
import sys

species = SPECIES.from_string(sys.argv[1])
annotation = ANNOTATIONS.NCBI if sys.argv[2] == 'NCBI' else ANNOTATIONS.ENSEMBL

genome = GenomeWorker(species, annotation, ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_CDS, SEQUENCE_LOAD.LOAD)
stats = AnalyzerData()

for chr_id in range(1, genome.chromosomes_count() + 1):
    genes_cnt = genome.genes_count_on_chr(chr_id)
    for i in range(0, genes_cnt):
        gene = genome.gene_by_ind(chr_id, i)
        CDSs = genome.get_gene_all_CDSs(gene)
        for CDS in CDSs:
            seq = genome.retrieve_interval_sequence(chr_id, CDS.start, CDS.end, CDS.strand)
            stats.analyze_sequence_stats(seq, CDS.frame, gene.id + " " + CDS.id)

output_data_path = f"./results/CDS control stats [{species.short_name()}, {annotation.short_name()}].txt"
with open(output_data_path, 'w') as f:
    f.write(stats.get_formatted_results(indent_tabs=0, use_percents=True))
