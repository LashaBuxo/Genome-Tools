from worker_genome import *
from worker_analyzer import *

genome = GenomeWorker(SPECIES.Homo_sapiens, ANNOTATIONS.ENSEMBL,
                      ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_CDS, SEQUENCE_LOAD.LOAD)

count = 0
for chr_id in range(1, genome.chromosomes_count() + 1):
    genes_cnt = genome.genes_count_on_chr(chr_id)
    for i in range(0, genes_cnt):
        gene = genome.gene_by_indexes(chr_id, i)
        transcript = genome.get_transcript_from_gene_by_criteria(gene.id, criteria=TRANSCRIPT_CRITERIA.LONGEST_CDS,
                                                                 tie_breaker_criteria=TRANSCRIPT_CRITERIA.RANDOM)
        seq = genome.retrieve_feature_sequence(chr_id, transcript)
        count += 1 if seq.startswith("ATGGGG") else 0

print(f'{count} transcripts starts with {"ATGGGG"}')
