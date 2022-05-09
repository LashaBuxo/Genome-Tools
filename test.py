from os.path import exists
import plotly.graph_objects as go
import glob
import pathlib

PWMs_common_path = "./motif_data/Homo_sapiens_2022_05_05_12_40_pm/pwms_all_motifs/"
sequences_path = "./overlapped_sequences.txt"

my_sequences_data = []

motif_PWMs = {}
motif_align_scores = {}


def read_all_sequences(path):
    file = open(path, 'r')
    lines = file.readlines()
    for line in lines:
        parts = line.replace('\n', '').split(' ')
        seq = []
        for i in range(0, len(parts[1])):
            seq.append(0 if parts[1][i] == 'A' else 1 if parts[1][i] == 'C' else 2 if parts[1][i] == 'G' else 3)
        my_sequences_data.append((parts[0], seq))


def get_hardcoded_PWM_from_file(path):
    file = open(path, 'r')
    lines = file.readlines()
    PWM = []
    for i in range(2, len(lines)):
        line = lines[i].replace('\n', '').split('\t')
        PWM.append((float(line[1]), float(line[2]), float(line[3]), float(line[4])))
    return PWM


def get_motif_three_best_alignment(PWM, id):
    motif_length = len(PWM)
    align_scores = []

    for sequence_data in my_sequences_data:
        gene_id = sequence_data[0]
        seq = sequence_data[1]
        l = 0

        while l + motif_length <= len(seq):
            align_score = 0
            for i in range(0, len(PWM)):
                align_score += PWM[i][seq[i + l]]
            align_score = align_score / motif_length if motif_length != 0 else 0
            align_scores.append((align_score, gene_id, l))
            l += 1
    align_scores.sort(key=lambda align: align[0], reverse=False)
    return [align_scores[0], align_scores[1], align_scores[2]]


read_all_sequences(sequences_path)

# All files and directories ending with .txt and that don't begin with a dot:
PWMs_path_list = glob.glob(f'{PWMs_common_path}*.txt')

for PWM_path in PWMs_path_list:
    motif_ID = pathlib.PurePath(PWM_path).name.replace('.txt', '')
    motif_PWMs[motif_ID] = get_hardcoded_PWM_from_file(PWM_path)
print(f'motifs loaded: {len(PWMs_path_list)}')

ind = 0
for (motif_id, motif_PWM) in motif_PWMs.items():
    ind += 1
    if ind % 1000 == 0: print(f'{ind}/{len(PWMs_path_list)}')
    motif_align_scores[motif_id] = get_motif_three_best_alignment(motif_PWM, motif_id)

sorted_motifs_by_score = sorted(motif_align_scores.items(),
                                key=lambda score_item: score_item[1][0][0] + score_item[1][1][0] + score_item[1][2][0],
                                reverse=True)

for i in range(0, 10):
    motif_id = sorted_motifs_by_score[i][0]
    align_data = sorted_motifs_by_score[i][1]

    score = align_data[0][0]
    print(f'{i}) {motif_id}:  score - {score}   genes - [{align_data[0][1]}/{align_data[0][2]},'
          f'{align_data[1][1]}/{align_data[1][2]},{align_data[2][1]}/{align_data[2][2]}]\n')
