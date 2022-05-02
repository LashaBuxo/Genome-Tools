# python ./generator_control_data.py Drosophila_melanogaster Ensembl


import plotly.graph_objects as go
from plotly.subplots import make_subplots
import sys
import math

from analyzer_worker import *
from genome_worker import *

species_list = [
    SPECIES.Homo_sapiens,
    SPECIES.Mus_musculus,
    SPECIES.Rattus_norvegicus,
    SPECIES.Drosophila_melanogaster,
    SPECIES.Danio_rerio,
]

annotation = ANNOTATIONS.NCBI if sys.argv[1] == 'NCBI' else ANNOTATIONS.ENSEMBL

# region GRAPHIC: variable description
species_color = ['#0CE83F', '#0CE83F', '#0CE83F', '#0CE83F', '#0CE83F']
figure_gc = make_subplots(rows=3, cols=2, subplot_titles=(str(species_list[0]), str(species_list[1]),
                                                          str(species_list[2]), str(species_list[3]),
                                                          str(species_list[4])))

figure_hd = make_subplots(rows=3, cols=2, subplot_titles=(str(species_list[0]), str(species_list[1]),
                                                          str(species_list[2]), str(species_list[3]),
                                                          str(species_list[4])))
# endregion

for specie_ind in range(len(species_list)):
    species = species_list[specie_ind]
    genome = GenomeWorker(species, annotation, ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_CDS, SEQUENCE_LOAD.LOAD)
    stats = AnalyzerData()

    gc_line_x, gc_line_y = [], []
    hd_line_x, hd_line_y = [], []

    for chr_id in range(1, genome.chromosomes_count() + 1):
        print(chr_id)
        genes_cnt = genome.genes_count_on_chr(chr_id)
        for gene_ind in range(0, genes_cnt):
            gene = genome.gene_by_ind(chr_id, gene_ind)
            transcripts = genome.get_gene_all_transcripts(gene)

            for transcript in transcripts:
                # make sure, only CDS features are loaded from the database at start
                CDSs_fragments = genome.get_transcript_all_fragments_sorted(transcript)
                # make sure, we have CDSs sorted, later we merge them with 'corresponding' frame
                CDSs_fragments.sort(key=lambda x: x.start)

                # generate/merge overall CDS from UTR5 (it must be reversed for '-' gene)
                overall_CDS_sequence = ""
                utr5_frame = CDSs_fragments[0].frame if gene.strand == '+' else CDSs_fragments[
                    len(CDSs_fragments) - 1].frame  # most cases it must be 0

                for i in range(0, len(CDSs_fragments)):
                    # from - strand, we need to merge from right
                    CDS = CDSs_fragments[i] if gene.strand == '+' else CDSs_fragments[len(CDSs_fragments) - 1 - i]
                    overall_CDS_sequence += genome.retrieve_interval_sequence(chr_id, CDS.start, CDS.end,
                                                                              CDS.strand)
                stats.analyze_sequence_stats(overall_CDS_sequence, frame=utr5_frame)

            gene_ind += 1

    output_data_path = f"./generated_data/control/CDS mean values [{species.short_name()}, {annotation.short_name()}].txt"
    with open(output_data_path, 'w') as f:
        f.write(stats.get_formatted_results(indent_tabs=0, use_percents=True))

    # region GRAPHIC: region HD and GC image generation
    gc_content = stats.get_regional_gc_contents()
    hd_content = stats.get_regional_hd_contents()
    gc_max_value, hd_max_value = 0, 0
    gc_min_value, hd_min_value = 100, 100
    for i in range(0, AnalyzerData.SUB_REGIONS):
        gc_line_x.append(i)
        gc_line_y.append(gc_content[i])
        hd_line_x.append(i)
        hd_line_y.append(hd_content[i])
        gc_max_value = max(gc_max_value, gc_content[i])
        hd_max_value = max(hd_max_value, hd_content[i])
        gc_min_value = min(gc_min_value, gc_content[i])
        hd_min_value = min(hd_min_value, hd_content[i])

    fig_row = 1 + specie_ind // 2
    fig_col = 1 + specie_ind % 2

    figure_gc.add_trace(go.Scatter(x=gc_line_x, y=gc_line_y, mode='lines', name="GC content",
                                   line=dict(color=species_color[specie_ind], width=1)), row=fig_row, col=fig_col)
    figure_hd.add_trace(go.Scatter(x=hd_line_x, y=hd_line_y, mode='lines', name="GC content",
                                   line=dict(color=species_color[specie_ind], width=1)), row=fig_row, col=fig_col)

    gc_mean_value = stats.get_mean_gc_content()
    hd_mean_value = stats.get_mean_hd_content()

    figure_gc.add_trace(go.Scatter(x=[0, AnalyzerData.SUB_REGIONS - 1], y=[gc_mean_value, gc_mean_value], mode='lines',
                                   name="mean value", line=dict(color='red', width=0.75, dash='dot')),
                        row=fig_row, col=fig_col)

    figure_hd.add_trace(go.Scatter(x=[0, AnalyzerData.SUB_REGIONS - 1], y=[hd_mean_value, hd_mean_value], mode='lines',
                                   name="mean value", line=dict(color='red', width=0.75, dash='dot')),
                        row=fig_row, col=fig_col)

    figure_gc.update_yaxes(
        tickmode='array',
        tickvals=[math.floor(gc_min_value), gc_mean_value, math.ceil(gc_max_value)],
        ticktext=["{:.1f}".format(math.floor(gc_min_value)) + "%", "<b>" + "{:.1f}".format(gc_mean_value) + "%</b>",
                  "{:.1f}".format(math.ceil(gc_max_value)) + "%"],
        row=fig_row, col=fig_col
    )
    figure_gc.add_hline(y=math.floor(gc_min_value), line_width=0.5, line_color='rgb(150,150,150)', row=fig_row,
                        col=fig_col)
    figure_gc.add_hline(y=math.ceil(gc_max_value), line_width=0.5, line_color='rgb(150,150,150)', row=fig_row,
                        col=fig_col)

    figure_hd.update_yaxes(
        tickmode='array',
        tickvals=[math.floor(hd_min_value), hd_mean_value,
                  math.ceil(hd_max_value)],
        ticktext=["{:.1f}".format(math.floor(hd_min_value)) + "%", "<b>" + "{:.1f}".format(hd_mean_value) + "%</b>",
                  "{:.1f}".format(math.ceil(hd_max_value)) + "%"],
        row=fig_row, col=fig_col
    )
    figure_hd.add_hline(y=math.floor(hd_min_value), line_width=0.5, line_color='rgb(150,150,150)', row=fig_row,
                        col=fig_col)
    figure_hd.add_hline(y=math.ceil(hd_max_value), line_width=0.5, line_color='rgb(150,150,150)', row=fig_row,
                        col=fig_col)
    # endregion

    stats.release_memory()
    genome.release_memory()

# region GRAPHIC: region HD and GC image generation

figure_gc.update_layout(title="Coding Sequences in species by GC-content", legend_title="lines", showlegend=False, )
figure_hd.update_layout(title="Coding Sequences in species by High-degeneracy amino acids", legend_title="lines",
                        showlegend=False)

figure_annotation = dict(xref='paper', yref='paper', x=0.5, y=-0.1, xanchor='center', yanchor='top',
                         text=f'annotation = {annotation.short_name()}, sub-regions = {AnalyzerData.SUB_REGIONS}',
                         font=dict(family='Arial', size=10, color='rgb(150,150,150)'),
                         showarrow=False)

figure_gc.add_annotation(figure_annotation)
figure_hd.add_annotation(figure_annotation)

figure_gc.update_annotations(font_size=10)
figure_hd.update_annotations(font_size=10)

figure_gc.update_xaxes(showticklabels=False)
figure_hd.update_xaxes(showticklabels=False)
figure_hd.update_yaxes(tickfont=dict(family='Rockwell', size=9), tickformat=',.1%')
figure_gc.update_yaxes(tickfont=dict(family='Rockwell', size=9), tickformat=',.1%')

gc_path = f'generated_data/control/GC content ({annotation.short_name()}, k={AnalyzerData.SUB_REGIONS}).png'
hd_path = f'generated_data/control/HD content ({annotation.short_name()}, k={AnalyzerData.SUB_REGIONS}).png'

figure_gc.write_image(gc_path, scale=4.0)
figure_hd.write_image(hd_path, scale=4.0)

# endregion
