import statistics
import numpy
from worker_genome import *
from GO_graph import *
from AgeData import*
import scipy.stats


specie_list = [SPECIES.Homo_sapiens,
               SPECIES.Mus_musculus,
               SPECIES.Drosophila_melanogaster,
               SPECIES.Arabidopsis_thaliana
               ]

for species in specie_list:
    # Load necessary data
    genome = GenomeWorker(species, ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_CDS, SEQUENCE_LOAD.NOT_LOAD)
    age_data = AgeData(species, genome)

    outfile = open(f"./generated_data/age_vs_length_data_{species.short_name()}.txt", "w")

    age_group_type = AgeData.AgeGroupsType.annotated_groups
    groups = age_data.calculate_groups(age_group_type)
    for group_index in range(len(groups)):
        group = groups[group_index]
        group_label = age_data.get_group_label(group, group_index, age_group_type)
        gene_ids = age_data.get_gene_ids_in_group(group)

        extra_data = ""

        length = []
        ov_count = 0
        for gene_id in gene_ids:
            gene = genome.feature_by_id(gene_id)
            #length.append(gene.end - gene.start + 1)

            transcript = genome.get_transcript_from_gene_by_criteria(gene.id, TRANSCRIPT_CRITERIA.LONGEST_CDS,
                                                                     TRANSCRIPT_CRITERIA.LONGEST)
            if transcript is None: continue
            frags = genome.get_fragments_from_transcript(transcript.id)

            cds_length = 0
            for frag in frags:
                if frag.featuretype == "CDS":
                    cds_length += frag.end - frag.start + 1
            length.append(cds_length)
        outfile.write(f"{group_label}\t{len(gene_ids)}\t{statistics.median(length)}\n")
    outfile.close()
