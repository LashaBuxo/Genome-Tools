from os.path import exists
import plotly.graph_objects as go

control_data_prefix = ".\generated_data\control\CDS mean values"
observed_data_prefix = ".\generated_data\observed\OGs by CDS"
summarized_data_prefix = ".\generated_data\summarized\summary"

overlaps_by_strands = ['diff_stranded', 'same_stranded', 'both_stranded']
overlaps_by_ORF = ['diff_ORF', 'same_ORF', 'both_ORF']
annotations = [
    'NCBI',
    'Ensembl'
]
species_list = ['Human', 'Mouse', 'Rat', 'Drosophila', 'ZebraFish', ]


# region hard-coded methods
def hard_coded_nucleotide_data(str):
    #		C - 28.55%	G - 28.55%	A - 21.45%	T - 21.45%
    cg_content = 0
    parts = str.split('\t')
    for k in range(0, len(parts)):
        if len(parts[k]) > 0 and (parts[k][0] == 'C' or parts[k][0] == 'G'):
            num = parts[k].split(' ')[2].replace('\n', '')
            num = float(num[0:len(num) - 1])
            cg_content += num
    cg_content = float("{:.2f}".format(cg_content))
    return cg_content


def hard_coded_degeneracy_data(str):
    #		medium - 36.32%	low - 34.58%	high - 29.1%
    hd_content = 0
    parts = str.split('\t')
    for k in range(0, len(parts)):
        if len(parts[k]) > 0 and parts[k][0] == 'h':
            num = parts[k].split(' ')[2].replace('\n', '')
            num = float(num[0:len(num) - 1])
            hd_content = num
    return hd_content


def hard_coded_genes_count(s: str):
    num = s.split(':')[1].replace(' ', '')
    return num


def hard_coded_sequence_length(s: str):
    num = s.split(':')[1].replace(' ', '').replace('\n', '')
    return int(num)


def color(color, text):
    return f"<span style='color:{str(color)}'> {str(text)} </span>"


def generate_cell_text_from_observed_data(observed_data, control_data):
    mean_gc = control_data[0]
    mean_hd = control_data[1]
    gc_content_gross = "{:.2f}".format(observed_data[2] - mean_gc) + '%'
    hd_content_gross = "{:.2f}".format(observed_data[3] - mean_hd) + '%'
    gc_content_gross = color('green', '+' + gc_content_gross) if observed_data[2] - mean_gc > 0 else color('red',
                                                                                                           gc_content_gross)

    hd_content_gross = color('green', '+' + hd_content_gross) if observed_data[3] - mean_hd > 0 else color('red',
                                                                                                           hd_content_gross)
    s = f'Overlapped genes: {observed_data[0]}<br>Overlapped length: {observed_data[1]}<br>' \
        f'GC-content: <b>{observed_data[2]}%</b> ({gc_content_gross})<br>HD-content: <b>{observed_data[3]}%</b> ({hd_content_gross})'
    return s


def read_necessary_observed_data(file_path):
    if not exists(file_path): return (-1, -1, -1, -1)
    file = open(file_path, 'r')
    lines = file.readlines()

    overlapped_genes = hard_coded_genes_count(lines[0])
    overlapped_seq_length = hard_coded_sequence_length(lines[1])

    cg_content = 0
    hd_content = 0
    for i in range(len(lines)):
        if lines[i].__contains__('nucleotides:') and not lines[i].__contains__('di-nucleotides:'):
            cg_content = hard_coded_nucleotide_data(lines[i + 1])
        if lines[i].__contains__('amino-acids by codon degeneracy:'):
            hd_content = hard_coded_degeneracy_data(lines[i + 1])
    return (overlapped_genes, overlapped_seq_length, cg_content, hd_content)


def read_necessary_control_data(file_path):
    if not exists(file_path): return (-1, -1)
    file = open(file_path, 'r')
    lines = file.readlines()
    gc_content = 0
    hd_content = 0
    for i in range(len(lines)):
        if lines[i].__contains__('nucleotides:') and not lines[i].__contains__('di-nucleotides:'):
            gc_content = hard_coded_nucleotide_data(lines[i + 1])
        if lines[i].__contains__('amino-acids by codon degeneracy:'):
            hd_content = hard_coded_degeneracy_data(lines[i + 1])
    return (gc_content, hd_content)


# endregion

for annotation in annotations:
    for species in species_list:

        values = [['<b>out-of-phase</b>', '<b>in-phase</b>', '<b>any phase</b>'], ['', '', ''], ['', '', ''],
                  ['', '', '']]

        control_file_path = f'{control_data_prefix} [{species}, {annotation}].txt'
        control_data = read_necessary_control_data(control_file_path)
        # if control_data[0] == -1 or control_data[1] == -1: continue

        for i in range(0, len(overlaps_by_strands)):
            strands_type = overlaps_by_strands[i]
            for j in range(0, len(overlaps_by_ORF)):
                ORF_type = overlaps_by_ORF[j]
                observed_file_path = f'{observed_data_prefix} [{species}, {annotation}, {strands_type}, {ORF_type}].txt'

                observed_data = read_necessary_observed_data(observed_file_path)
                if observed_data[0] == -1 or observed_data[1] == -1 or observed_data[2] == -1 or observed_data[3] == -1:
                    continue

                values[i + 1][j] = generate_cell_text_from_observed_data(observed_data, control_data)

        fig = go.Figure(data=[go.Table(
            # columnorder=[1, 2, 3, 4],
            columnwidth=[40, 100, 100, 100],
            header=dict(
                values=[[f'<b>OGs by CDS<br>{species}</b> <b>({annotation})</b>'],
                        ['<b>diff stranded</b>'],
                        ['<b>same strands</b>'],
                        ['<b>any stranded</b>']],
                line_color='darkslategray',
                fill=dict(color=['royalblue', 'paleturquoise', 'paleturquoise', 'paleturquoise']),
                align=['center', 'center', 'center', 'center'],
                font=dict(color='black', size=12),
                # height=40
            ),
            cells=dict(
                values=values,
                line_color='darkslategray',
                fill=dict(color=['paleturquoise', 'white', 'white', 'white']),
                align=['center', 'left', 'left', 'left'],
                font=dict(size=[12, 10, 10, 10], color=['black', 'black', 'black', 'black']),
                height=30
            )
        )
        ])
        summarized_data_path = f'{summarized_data_prefix} [{species}, {annotation}].png'
        fig.write_image(summarized_data_path,
                        width=1400, height=500,
                        scale=2.0)
