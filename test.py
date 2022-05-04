import plotly.graph_objects as go
from plotly.subplots import make_subplots
import sys
import math

from analyzer_worker import *
from genome_worker import *

genome = GenomeWorker(SPECIES.Homo_sapiens, ANNOTATIONS.ENSEMBL, ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_FRAGMENTS,
                      SEQUENCE_LOAD.LOAD)
stats = AnalyzerData()

for chr_id in range(1, genome.chromosomes_count() + 1):
    genes_cnt = genome.genes_count_on_chr(chr_id)
    for i in range(genes_cnt):
        gene = genome.gene_by_ind(chr_id, i)
        transcripts = genome.get_gene_all_transcripts(gene)
        for transcript in transcripts:
            fragments = genome.get_transcript_all_fragments_sorted(transcript)
            for fragment in fragments:
                if fragment.featuretype == 'five_prime_UTR':
                    seq = genome.retrieve_interval_sequence(chr_id, fragment.start, fragment.end, fragment.strand)
                    stats.analyze_sequence_stats(seq, 0)

print(stats.get_formatted_results(indent_tabs=0, use_percents=True))
