import numpy

from worker_genome import *
from GO_graph import *
import scipy.stats

specie_list = [
    SPECIES.Homo_sapiens,
    SPECIES.Mus_musculus,
    SPECIES.Drosophila_melanogaster,
    SPECIES.Arabidopsis_thaliana
]

ov_types = [REVIEWED_OVERLAP_TYPE.CONVERGENT, REVIEWED_OVERLAP_TYPE.DIFF_NESTED, REVIEWED_OVERLAP_TYPE.DIVERGENT,
            REVIEWED_OVERLAP_TYPE.PROMOTOR]

for species in specie_list:
    # Load necessary data
    genome = GenomeWorker(species, ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_CDS, SEQUENCE_LOAD.NOT_LOAD)

    OGs_by_type = {REVIEWED_OVERLAP_TYPE.CONVERGENT: {},
                   REVIEWED_OVERLAP_TYPE.DIFF_NESTED: {},
                   REVIEWED_OVERLAP_TYPE.DIVERGENT: {},
                   REVIEWED_OVERLAP_TYPE.PROMOTOR: {}}

    for chr_index in range(1, genome.chromosomes_count() + 1):
        genes_cnt = genome.genes_count_on_chr(chr_index)
        for i in range(0, genes_cnt):
            gene1 = genome.gene_by_indexes(chr_index, i)
            for j in range(i + 1, genes_cnt):
                gene2 = genome.gene_by_indexes(chr_index, j)
                ov_type = genome.get_features_REVIEWED_overlap_type(gene1, gene2)
                if ov_type != OVERLAP_TYPE.NONE:
                    OGs_by_type[ov_type][gene1.id] = True
                    OGs_by_type[ov_type][gene2.id] = True

    outfile = open(f"./generated_data/length(CDS)_vs_OGs_data_{species.short_name()}.txt", "w")

    length_and_gene_id = []
    for chr_index in range(1, genome.chromosomes_count() + 1):
        genes_cnt = genome.genes_count_on_chr(chr_index)
        for i in range(0, genes_cnt):
            gene = genome.gene_by_indexes(chr_index, i)
            transcript = genome.get_transcript_from_gene_by_criteria(gene.id, TRANSCRIPT_CRITERIA.LONGEST_CDS,
                                                                     TRANSCRIPT_CRITERIA.LONGEST)
            if transcript is None: continue
            cds_length = genome.get_transcript_CDS_length(transcript.id)
            length_and_gene_id.append((cds_length, gene.id))

            # length = gene.end - gene.start + 1
            # length_and_gene_id.append((length, gene.id))

    groups = []
    labels = []

    groups_count = 30
    length_and_gene_id.sort(key=lambda x: x[0])
    steps = len(length_and_gene_id) // groups_count
    for group_index in range(groups_count):
        l = group_index * steps
        r = (group_index + 1) * steps - 1 if group_index + 1 < groups_count else len(length_and_gene_id) - 1
        groups.append([])
        labels.append(f"{length_and_gene_id[l][0]}-{length_and_gene_id[r][0]}")
        for i in range(l, r + 1):
            groups[group_index].append(length_and_gene_id[i][1])

    for group_index in range(len(groups)):
        gene_ids = groups[group_index]
        group_label = labels[group_index]

        extra_data = ""
        for ov_type in ov_types:
            overlapped_ids = OGs_by_type[ov_type]
            ov_count = 0
            for gene_id in gene_ids:
                if overlapped_ids.__contains__(gene_id):
                    ov_count += 1
            extra_data += f"\t{ov_count}"
        outfile.write(f"{group_label}\t{len(gene_ids)}{extra_data}\n")
    outfile.close()
