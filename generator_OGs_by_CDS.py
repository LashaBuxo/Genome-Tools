from overlaps_by_CDS_finder import *
from generator_writer import *

start_time = time.time()

assert len(sys.argv) == 2
species = SPECIES.from_string(sys.argv[1])

# Load necessary annotation data
genome = GenomeWorker(species, ANNOTATIONS.ENSEMBL, ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_FRAGMENTS,
                      SEQUENCE_LOAD.LOAD)

OGs_by_any_ISO_graph, OGs_records = overlaps_finder(genome, SEARCH_METHOD.BY_GENE_CDSs)

for interaction in OVERLAP_INTERACTION.get_overlap_interactions():
    addon = '(ATI)' if interaction is OVERLAP_INTERACTION.ATI else ''

    g_clusters = OGs_by_any_ISO_graph.get_connected_clusters(interaction)

    output_texts_path = f'generated_data/texts/OGs by CDS {addon} [{species.short_name()}].txt'
    with open(output_texts_path, 'w') as file:
        if interaction == OVERLAP_INTERACTION.ATI:
            write_summarised_data(genome, file, g_clusters)
            write_clusters_data(genome, file, g_clusters)
        else:
            write_summarised_data(genome, file, g_clusters, OGs_records=OGs_records, with_details=True)
            write_clusters_data(genome, file, g_clusters, OGs_records=OGs_records, with_records=True)

    output_plot_path = f'generated_data/plots/'

    output_plots_path = f'generated_data/plots/OGs by CDS {addon} [{species.short_name()}].png'
