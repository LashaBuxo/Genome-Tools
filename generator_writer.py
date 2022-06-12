from overlap_by_CDS_record import *
import plotly.express as px
import pandas as pd
import plotly.graph_objects as go


def match_score(s1: str, s2: str):
    s1 = s1.replace('.', '').replace('|', '')
    s2 = s2.replace('.', '').replace('|', '')
    dp = []
    for i in range(0, len(s1) + 1):
        arr = [0] * (len(s2) + 1)
        dp.append(arr)
    for i in range(1, len(s1) + 1):
        for j in range(1, len(s2) + 1):
            dp[i][j] = max(max(dp[i - 1][j], dp[i][j - 1]),
                           dp[i - 1][j - 1] + 1 if s1[i - 1] == s2[j - 1] else dp[i - 1][j - 1])

    return dp[len(s1)][len(s2)]


def write_excel_OGs_by_CDS_data(file, genome: GenomeWorker, mouse: GenomeWorker, OGs_records: list[OverlapRecord],
                                OGs_records_mouse: list[OverlapRecord], OGs_random_records: list[OverlapRecord],
                                control_groups_data,
                                control_homologue_species, observed_homologue_species, control_max_homologue_species,
                                observed_max_homologue_species
                                ):
    file.write('\n--------------------------EXCEL EXPORT DATA-----------------------------------\n\n')
    file.write(
        'gene1\tgene2\tt1\tt2\tt_ov\tex_ov\tov_len\tGC\tCONS1\tCONS2\tt_ov\tex_ov\tov_len\tGC\tCONS1\tCONS2\thas_mM?\tov_mM?\tov_len(mM)\tmatch_score\n')

    for index in range(0, len(OGs_records)):
        record = OGs_records[index]
        sym1, sym2 = record.gene1_symbol, record.gene2_symbol

        transcript1, transcript2 = record.transcript1, record.transcript2

        cons1, cons2 = record.get_transcript_conservation_score(
            transcript1.id), record.get_transcript_conservation_score(
            transcript2.id)

        presentedMM = "'+" if mouse.are_genes_presented(sym1, sym2) else "'-"
        overlappedMM = "'+" if mouse.are_genes_overlapped(sym1, sym2) else "'-"

        # complete1 = "'+" if genome.is_transcript_complete(transcript1.id) else "'-"
        # complete2 = "'+" if genome.is_transcript_complete(transcript2.id) else "'-"

        trans_ov_type, exon_ov_type = record.transcriptic_overlap_type.short_name(), record.exonic_overlap_type.short_name()
        ov_length = record.get_record_length()

        gc = record.get_record_GC_content()

        co_overlapped_len = 0
        ov_seq2 = ''
        for record_mouse in OGs_records_mouse:
            if record.gene1_symbol != record_mouse.gene1_symbol and record.gene1_symbol != record_mouse.gene2_symbol:
                continue
            if record.gene2_symbol != record_mouse.gene1_symbol and record.gene2_symbol != record_mouse.gene2_symbol:
                continue
            co_overlapped_len = record_mouse.get_record_length()
            ov_seq2 = record_mouse.get_overlapped_sequence()
        ov_seq1 = record.get_overlapped_sequence()
        score = match_score(ov_seq1, ov_seq2)

        random_record = OGs_random_records[index]
        random_length = random_record.get_record_length()
        random_GC = random_record.get_record_GC_content()
        random_cons1, random_cons2 = random_record.get_transcript_conservation_score(random_record.transcript1.id), \
                                     random_record.get_transcript_conservation_score(random_record.transcript2.id)
        random_trans_ov_type, random_exon_ov_type = random_record.transcriptic_overlap_type.short_name(), \
                                                    random_record.exonic_overlap_type.short_name()
        file.write(
            f"{sym1}\t{sym2}\t{transcript1.id.replace('transcript:', '')}\t{transcript2.id.replace('transcript:', '')}\t"
            f"{random_trans_ov_type}\t{random_exon_ov_type}\t{random_length}\t{random_GC}\t{random_cons1}\t{random_cons2}"
            f"\t{trans_ov_type}\t{exon_ov_type}\t{ov_length}\t{gc}\t{cons1}\t{cons2}\t{presentedMM}\t{overlappedMM}\t"
            f"{co_overlapped_len}\t{score}\n")

    file.write('\n\n')
    for data in control_groups_data:
        file.write(f"{data}\n")

    arr = list(control_homologue_species.items())
    arr.sort(key=lambda x: x[1], reverse=True)
    dict = {}
    for key, val in arr:
        dict[key] = (val, 0, 0, 0)
    arr = list(observed_homologue_species.items())
    for key, val in arr:
        if dict.__contains__(key):
            dict[key] = (dict[key][0], val, 0, 0)
        else:
            dict[key] = (0, val, 0, 0)
    arr = list(control_max_homologue_species.items())
    for key, val in arr:
        if dict.__contains__(key):
            dict[key] = (dict[key][0], dict[key][1], val, 0)
        else:
            dict[key] = (0, 0, val, 0)
    arr = list(observed_max_homologue_species.items())
    for key, val in arr:
        if dict.__contains__(key):
            dict[key] = (dict[key][0], dict[key][1], dict[key][2], val)
        else:
            dict[key] = (0, 0, 0, val)
    arr = list(dict.items())
    file.write(
        '\n\nspecies\thomologue count (control)\thomologue count (observed)\thomologue max count (control)\thomologue max count (observed)\n')
    for key, (val1, val2, val3, val4) in arr:
        file.write(f'{key}\t{val1}\t{val2}\t{val3}\t{val4}\n')

    file.write('\n------------------------------------------------------------------------------\n\n')


def write_excel_formatted_general_data(file, genome: GenomeWorker, gene_clusters):
    file.write('\n--------------------------EXCEL EXPORT DATA-----------------------------------\n\n')
    file.write('chr\tsize\tgenes\tOGs\tmito_genes\tmito_OGs\n')
    OGs, mito_genes, OGs_mito = {}, {}, {}

    mito_genes_list = []
    mito_OGs_list = []
    for cluster in gene_clusters:
        for gene_id in cluster:
            chr_id = genome.get_feature_chromosomal_position(gene_id)[0]
            OGs[chr_id] = OGs.get(chr_id, 0) + 1
            if genome.is_gene_MITO(gene_id):
                OGs_mito[chr_id] = OGs_mito.get(chr_id, 0) + 1
                gene_symbol = genome.get_gene_symbol(genome.feature_by_id(gene_id))
                mito_OGs_list.append(gene_symbol)

    chrs = genome.chromosomes_count()
    for chr_id in range(1, chrs + 1):
        chr_name = f'{chr_id}' if chr_id < chrs - 2 else 'X' if chr_id == chrs - 2 else 'Y' if chr_id == chrs - 1 else 'MT'
        chr_size = genome.chromosome_length(chr_id)
        genes_cnt = genome.genes_count_on_chr(chr_id)
        for ind in range(0, genes_cnt):
            gene = genome.gene_by_indexes(chr_id, ind)
            if genome.is_gene_MITO(gene.id):
                mito_genes[chr_id] = mito_genes.get(chr_id, 0) + 1
                mito_genes_list.append(genome.get_gene_symbol(gene))
        file.write(
            f'{chr_name}\t{chr_size}\t{genes_cnt}\t{OGs.get(chr_id, 0)}\t{mito_genes.get(chr_id, 0)}\t{OGs_mito.get(chr_id, 0)}\n')

    file.write('\nmito_genes list\n')
    for gene_sym in mito_genes_list:
        file.write(gene_sym + "\n")
    file.write('\nmito_OGs list\n')
    for gene_sym in mito_OGs_list:
        file.write(gene_sym + "\n")
    file.write('\n------------------------------------------------------------------------------\n\n')


LABELS_PARAMS = {
    'control (Divergent)': ("კონტროლ ჯგუფი (დივ.)", 'darkcyan'),
    'control (Convergent)': ("კონტროლ ჯგუფი (კონვ.)", 'darkcyan'),
    'control (Nested (same strand))': ("კონტროლ ჯგუფი (ჩადგმ.-)", 'darkcyan'),
    'control (Tandem)': ("კონტროლ ჯგუფი (ტანდ.)", 'darkcyan'),
    'control (Nested (diff strand))': ("კონტროლ ჯგუფი (ჩადგმ.+)", 'darkcyan'),

    'Convergent': ('კონვერგენციული', '#ff6f31')
    , 'Divergent': ('დივერგენციული', 'cornflowerblue'),
    'Tandem': ('ტანდემური', '#7e3794')
    , 'Nested (same strand)': ('ჩადგმული (იმავე ჯაჭვზე)', '#004561'),
    'Nested (diff strand)': ('ჩადგმული (სხვა ჯაჭვზე)', '#4cdc8b')}


def build_2d_plot(data, path):
    x, y, types, labels = [], [], [], []
    for item in data:
        x.append(item[0])
        y.append(item[1])
        types.append(LABELS_PARAMS[item[2]][0])
        labels.append(item[3])
    data3 = pd.DataFrame({'x_var': x, 'y_var': y, 'types': types, 'labels': labels})

    col_map = {}
    for key, (geo_label, color) in LABELS_PARAMS.items():
        col_map[geo_label] = color
    fig = px.scatter(data_frame=data3, x='x_var', y='y_var', color_discrete_map=col_map, color='types', text='labels')

    fig.update_traces(textposition='middle right', marker_size=10)
    # fig.update_yaxes(gridcolor='whi')
    # fig.update_xaxes(showline=True, linewidth=1, linecolor='grey')
    # fig.update_yaxes(showline=True, linewidth=1, linecolor='grey')

    fig.update_layout(
        # title={'text': "<b>CDS-ით გადაფარული გენები</b>", 'xanchor': 'center', 'y': 0.95, 'x': 0.5,  'yanchor': 'bottom'},
        yaxis_title='კონსერვაციის ხარისხი 1',
        xaxis_title='კონსერვაციის ხარისხი 2',
        legend_title='გადაფარვის ტიპი',
        font_family="Arial", font_color='black', font_size=24,
        paper_bgcolor="white", plot_bgcolor='#e5ecf6',
        legend=dict(yanchor="top", y=0.99, xanchor="right", x=0.99))

    fig.add_vline(x=3.5, line_color='white', line_width=0.75)
    fig.show()
    fig.write_image(path, scale=2.0, width=1500, height=750)


def build_histogram(data, mean_x, legend_title, annot_title, l, r, ax, ay, path):
    fig = go.Figure(data=[
        go.Histogram(x=data, name=legend_title, histnorm='probability', nbinsx=10,
                     marker=dict(color='lightblue'))])
    fig.update_layout(legend=dict(font=dict(size=25), yanchor="top", y=0.99, xanchor="left", x=0.01),
                      title_text=legend_title, showlegend=True)
    fig.update_yaxes(showticklabels=False, )
    fig.update_xaxes(tickfont=dict(size=20), range=[l, r])
    fig.add_vline(x=mean_x, line_color='red', line_width=1)

    fig.add_annotation(
        x=mean_x, y=0.18, xref="x", yref="y", text=annot_title, showarrow=True,
        font=dict(family="Courier New, monospace", size=25, color="black"),
        align="center", arrowhead=2, arrowsize=1, arrowwidth=2, arrowcolor="black", ax=ax,
        ay=ay, bordercolor="black", borderwidth=2, borderpad=4, opacity=0.8
    )

    fig.write_image(path, scale=2.0, width=1500, height=750)


# region clusters data writers


def write_clusters_data(genome: GenomeWorker, file, clusters, OGs_records=None, with_records=False):
    file.write(f'--------------------------clusters and their parameters-----------------------------------\n\n')
    for i in range(len(clusters)):
        cluster = clusters[i]
        file.write(
            f'cluster #{i + 1}  |  size: {len(cluster)}  | chr: {genome.get_feature_chromosomal_position(cluster[0])[1]}\n')
        for gene_index in range(len(cluster)):
            g_id = cluster[gene_index]
            g = genome.feature_by_id(g_id)
            isos = len(genome.get_transcripts_from_gene(g.id))

            g_name = genome.get_gene_symbol(g)
            g_desc = genome.get_gene_description(g)
            mito_status = 'Yes' if genome.is_gene_MITO(g_id) else 'No'
            file.write(
                f'\t{g_name} = strand {g.strand} | isMitoGene: {mito_status} | isoforms: {isos} | desc: {g_desc}\n')

        if with_records: __write_cluster_record(file, cluster, OGs_records, 1)
        file.write('\n')


def __write_cluster_record(file, cluster, OGs_records, indent_tabs=0):
    for record in OGs_records:
        belongs_to_cluster = False
        for g_id in cluster:
            if record.is_record_about_gene(g_id): belongs_to_cluster = True
        if not belongs_to_cluster: continue
        file.write(record.get_formatted_stats(indent_tabs))


# endregion


# region summarised data writers

def write_summarised_data(genome, file, gene_clusters, OGs_records=None, with_details=False):
    OGs, OGs_count, OGs_count_formatted, og_nuclear_mito_proteome, og_mtDNA_mito_proteome = __get_OGs_count_by_cluster(
        gene_clusters, genome)

    file.write(f'--------------------------summarised values---------------------------------\n')
    file.write(f'total overlapped (by CDS) genes: {OGs_count} = {OGs_count_formatted}\n\n')
    file.write(
        f"\tmito proteome coding genes: {og_nuclear_mito_proteome + og_mtDNA_mito_proteome} = {og_nuclear_mito_proteome} (nuclear DNA) + {og_mtDNA_mito_proteome} (mtDNA)\n\n")

    if not with_details:  # if we don't require other information
        return

    diff_strand_OGs, same_strand_OGs = __get_OGs_count_by_strands(OGs_records, OGs)
    OGs_len_est = __get_OGs_gene_length_estimation(genome, OGs)
    OGs_iso_ests = __get_OGs_isoforms_estimations(genome, OGs_records, OGs)

    file.write(f'diff-strand overlapped (by CDS) genes: {diff_strand_OGs}\n')
    file.write(f'same-strand overlapped (by CDS) genes with different frame: {same_strand_OGs}\n\n')

    file.write(f'overlapped gene length estimation (mean, std, median, interquartile): {OGs_len_est}\n\n')

    file.write(f'isoforms per overlapped gene (mean, std, median, interquartile): {OGs_iso_ests[0]}\n')
    file.write(f'overlapped isoforms per overlapped gene (mean, std, median, interquartile): {OGs_iso_ests[1]}\n')
    file.write(f'non-Overlapped isoforms per overlapped gene (mean, std, median, interquartile): {OGs_iso_ests[2]}\n\n')


def __get_OGs_count_by_cluster(clusters, genome):
    OGs_by_cluster_size = {}
    OGs_count = 0
    OGs = []
    og_nuclear_mito_proteome = 0
    og_mtDNA_mito_proteome = 0

    for cluster in clusters:
        OGs_count += len(cluster)
        for g_id in cluster:
            OGs.append(g_id)
            is_mito = genome.get_feature_chromosomal_position(g_id)[0] == genome.chromosomes_count()
            if genome.is_gene_MITO(g_id) and is_mito: og_mtDNA_mito_proteome += 1
            if genome.is_gene_MITO(g_id) and not is_mito: og_nuclear_mito_proteome += 1
        OGs_by_cluster_size[len(cluster)] = OGs_by_cluster_size.get(len(cluster), 0) + 1
    OGs_count_formatted = ''
    for size in sorted(list(OGs_by_cluster_size.keys())):
        OGs_count_formatted += f'{" +" if len(OGs_count_formatted) > 0 else ""} {OGs_by_cluster_size[size]}x[{size}]'
    return OGs, OGs_count, OGs_count_formatted, og_nuclear_mito_proteome, og_mtDNA_mito_proteome


def __get_OGs_isoforms_estimations(genome, OGs_records, OGs_ids):
    genes_isoforms = []
    genes_overlapped_isoforms = []
    genes_non_overlapped_isoforms = []
    for g_id in OGs_ids:
        gene = genome.feature_by_id(g_id)
        transcripts = genome.get_transcripts_from_gene(gene.id)
        genes_isoforms.append(len(transcripts))
        cnt = 0
        for record in OGs_records:
            if record.is_record_about_gene(g_id): cnt += 1
        genes_overlapped_isoforms.append(cnt)
        genes_non_overlapped_isoforms.append(len(transcripts) - cnt)

    est1 = get_value_estimation(genes_isoforms)
    est2 = get_value_estimation(genes_overlapped_isoforms)
    est3 = get_value_estimation(genes_non_overlapped_isoforms)
    return f'{est1[0]}, {est1[1]}, {est1[2]}, {est1[3]}', f'{est2[0]}, {est2[1]}, {est2[2]}, {est2[3]}', f'{est3[0]}, {est3[1]}, {est3[2]}, {est3[3]}'


def __get_OGs_gene_length_estimation(genome, OGs_ids):
    genes_length = []
    for g_id in OGs_ids:
        gene = genome.feature_by_id(g_id)
        genes_length.append(gene.end - gene.start + 1)
    est = get_value_estimation(genes_length)
    return f'{est[0]}, {est[1]}, {est[2]}, {est[3]}'


def __get_OGs_count_by_strands(OGs_records, OGs_ids):
    genes_having_diff_strand_overlap = 0
    genes_having_same_strand_overlap = 0
    for g_id in OGs_ids:
        has_diff_strand_overlap = False
        for record in OGs_records:
            if record.is_record_about_gene(g_id):
                if record.is_record_diff_strand(): has_diff_strand_overlap = True
        genes_having_diff_strand_overlap += 1 if has_diff_strand_overlap else 0
        genes_having_same_strand_overlap += 1 if not has_diff_strand_overlap else 0
    return genes_having_diff_strand_overlap, genes_having_same_strand_overlap

# endregion
