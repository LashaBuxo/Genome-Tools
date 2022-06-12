import scipy
from mlxtend.evaluate import permutation_test
from overlaps_by_CDS_finder import *
from generator_writer import *
import numpy as np

start_time = time.time()

# Load necessary annotation data
genome = GenomeWorker(SPECIES.Homo_sapiens, ANNOTATIONS.ENSEMBL, ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_FRAGMENTS,
                      SEQUENCE_LOAD.LOAD)

mouse_genome = GenomeWorker(SPECIES.Mus_musculus, ANNOTATIONS.ENSEMBL,
                            ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_FRAGMENTS, SEQUENCE_LOAD.LOAD)


def random_CDS_overlap(record_A: OverlapRecord, record_B: OverlapRecord):
    assert record_A is not None or record_B is not None
    if record_A is None: return record_B
    if record_B is None: return record_A
    return record_A if random.randint(0, 1000000) < random.randint(0, 1000000) else record_B


def most_conserved_CDS_overlap(record_A: OverlapRecord, record_B: OverlapRecord):
    assert record_A is not None or record_B is not None
    if record_A is None: return record_B
    if record_B is None: return record_A

    consA_1 = record_A.get_transcript_conservation_score(record_A.transcript1.id)
    consA_2 = record_A.get_transcript_conservation_score(record_A.transcript2.id)
    consB_1 = record_B.get_transcript_conservation_score(record_B.transcript1.id)
    consB_2 = record_B.get_transcript_conservation_score(record_B.transcript2.id)

    if min(consA_1, consA_2) == min(consB_1, consB_2):
        if max(consA_1, consA_2) == max(consB_1, consB_2):
            return record_A if record_A.get_record_length() > record_B.get_record_length() else record_B
        else:
            return record_A if max(consA_1, consA_2) > max(consB_1, consB_2) else record_B
    else:
        return record_A if min(consA_1, consA_2) > min(consB_1, consB_2) else record_B


mouse_graph, OGs_records_mouse = overlaps_finder(mouse_genome, SEARCH_METHOD.BY_GENE_CDSs_EXCEPT_ATI, min_length=1,
                                                 record_compare_function=most_conserved_CDS_overlap)
mouse_clusters = mouse_graph.get_connected_clusters()
cnt = 0
for c in mouse_clusters:
    cnt += len(c)
print(f'Mouse CDS overlapped genes: {cnt} clusters: {len(mouse_clusters)}')

OGs_by_any_ISO_graph, OGs_records = overlaps_finder(genome, SEARCH_METHOD.BY_GENE_CDSs_EXCEPT_ATI, min_length=1,
                                                    record_compare_function=most_conserved_CDS_overlap)
_, OGs_random_records = overlaps_finder(genome, SEARCH_METHOD.BY_GENE_CDSs_EXCEPT_ATI, min_length=1,
                                        record_compare_function=random_CDS_overlap)

g_clusters = OGs_by_any_ISO_graph.get_connected_clusters()

OGs_random_records.sort(key=lambda x: x.get_record_length(), reverse=True)
OGs_records.sort(key=lambda x: x.get_record_length(), reverse=True)

output_texts_path = f'generated_data/texts/OGs by CDS (filtered).txt'
file = open(output_texts_path, 'w')
write_summarised_data(genome, file, g_clusters, OGs_records=OGs_records, with_details=True)
write_clusters_data(genome, file, g_clusters, OGs_records=OGs_records, with_records=True)

output_plot_path = f'generated_data/plots/'

overlapped = {}
for c in g_clusters:
    for g in c:
        overlapped[g] = True

transcripts_by_cds_length = []
for chr_id in range(1, genome.chromosomes_count() + 1):
    genes_cnt = genome.genes_count_on_chr(chr_id)
    for i in range(0, genes_cnt):
        gene = genome.gene_by_indexes(chr_id, i)
        # if overlapped.__contains__(gene.id): continue
        transcripts = genome.get_transcripts_from_gene(gene.id)
        for transcript in transcripts:
            CDS_length = genome.get_transcript_CDS_length(transcript.id)
            transcripts_by_cds_length.append((transcript.id, CDS_length))
transcripts_by_cds_length.sort(key=lambda x: x[1])


def overlap_starting_pos(transcript_id, ov_segments):
    transcript = genome.feature_by_id(transcript_id)
    start_genomic_pos = ov_segments[0][0][0] if transcript.strand == '+' else ov_segments[len(ov_segments) - 1][0][1]

    frags = genome.get_fragments_from_transcript(transcript_id)

    cur_cds_index = 0
    for index in range(0, len(frags)):
        frag = frags[index] if transcript.strand == '+' else frags[len(frags) - index - 1]
        if frag.featuretype != 'CDS': continue
        if not genome.are_segments_overlapped((frag.start, frag.end), (start_genomic_pos, start_genomic_pos)):
            cur_cds_index += frag.end - frag.start + 1
        else:
            cur_cds_index += start_genomic_pos - frag.start + 1 if transcript.strand == '+' else frag.end - start_genomic_pos + 1
            break
    return cur_cds_index


def get_transcript_similar_CDS(control_id, ov_segments, ov_id, ov_type):
    control_transcript = genome.feature_by_id(control_id)
    chr_id = genome.get_feature_chromosomal_position(control_transcript.id)[0]
    frags = genome.get_fragments_from_transcript(control_transcript.id)

    cds = ''
    for ind in range(0, len(frags)):
        frag = frags[ind] if control_transcript.strand == '+' else frags[len(frags) - ind - 1]
        if frag.featuretype != 'CDS': continue
        cds += genome.retrieve_interval_sequence(chr_id, frag.start, frag.end, control_transcript.strand)

    # if ov_type == OVERLAP_TYPE.DIVERGENT:
    #     start_cds_index = 1
    # else:
    #
    start_cds_index = overlap_starting_pos(ov_id, ov_segments)

    ov_length = 0
    for (l, r), _ in ov_segments:
        ov_length += r - l + 1

    return cds[start_cds_index - 1:start_cds_index + ov_length - 1]


def collect_control_data(rec: OverlapRecord, transcript_id, ov_type):
    observed_length = genome.get_transcript_CDS_length(transcript_id)

    desired_lower_limit, desired_upper_limit = observed_length, observed_length * 1.1
    l = 0
    r = len(transcripts_by_cds_length) - 1
    while l < r:
        mid = (l + r) // 2
        if transcripts_by_cds_length[mid][1] < desired_lower_limit:
            l = mid + 1
        else:
            r = mid
    index_l = l
    l = 0
    r = len(transcripts_by_cds_length) - 1
    while l < r:
        mid = (l + r) // 2
        if transcripts_by_cds_length[mid][1] < desired_upper_limit:
            l = mid + 1
        else:
            r = mid
    index_r = l

    control_transcript_id = transcripts_by_cds_length[random.randint(index_l, index_r)][0]
    cds = get_transcript_similar_CDS(control_transcript_id, rec.overlapped_segments, transcript_id, ov_type)
    assert len(cds) == rec.get_record_length()
    stats = AnalyzerData()
    stats.analyze_sequence_stats(cds, 0)
    gc = stats.get_gc_content()

    random_transcript = transcripts_by_cds_length[random.randint(0, len(transcripts_by_cds_length) - 1)][0]
    cons = genome.get_transcript_conservation_info(random_transcript)
    hom_species = genome.get_transcript_homologue_species(random_transcript)
    parent_id = genome.feature_by_id(random_transcript).attributes['Parent'][0]
    hom_max_species = genome.get_gene_max_conserved_homologue_species(parent_id)
    return gc, cons, hom_species, hom_max_species


RECORD_MIN_LENGTH = 60
N = len(OGs_records)
RUNS = 1000

control_groups_data = []

observed_homologue_species = {}
observed_max_homologue_species = {}
control_homologue_species_runs = []
control_max_homologue_species_runs = []

GC_means, CONS_medians = [], []
GC_observed_X = 0
GC_observed_count = 0
CONS_observed = []


def add_homologs(hom_species, data):
    for x in hom_species:
        if not data.__contains__(x): data[x] = 0
        data[x] += 1
    data['count'] = data.get('count', 0) + 1


def normalize_runs_and_counts(run_results):
    total_specie_count = {}
    for r in run_results:
        for specie, count in r.items():
            if specie == 'count': continue
            if not total_specie_count.__contains__(specie): total_specie_count[specie] = 0
            total_specie_count[specie] += count / r['count']  # normalize counts, by transcripts count
    for specie, count in total_specie_count.items():
        total_specie_count[specie] = count / len(run_results)  # normalize runs
    return total_specie_count


for run in range(0, RUNS):
    gc = {'All': [], 'Convergent': [], 'Divergent': [], 'Nested (diff strand)': [], 'Tandem': [],
          'Nested (same strand)': []}
    cn = {'All': [], 'Convergent': [], 'Divergent': [], 'Nested (diff strand)': [], 'Tandem': [],
          'Nested (same strand)': []}
    cn_zeroes = {'All': 0, 'Convergent': 0, 'Divergent': 0, 'Nested (diff strand)': 0, 'Tandem': 0,
                 'Nested (same strand)': 0}

    control_homologue_species_run = {}
    control_max_homologue_species_run = {}

    for record in OGs_random_records:
        ov_type = record.transcriptic_overlap_type.short_name()
        gc1, cons1, homs1, homs_max1, = collect_control_data(record, record.transcript1.id,
                                                             record.transcriptic_overlap_type)
        gc2, cons2, homs2, homs_max2 = collect_control_data(record, record.transcript2.id,
                                                            record.transcriptic_overlap_type)
        if record.get_record_length() >= RECORD_MIN_LENGTH:
            gc['All'].append(gc1)
            gc['All'].append(gc2)
            gc[ov_type].append(gc1)
            gc[ov_type].append(gc2)

        if cons1 != 0:
            cn['All'].append(cons1)
            cn[ov_type].append(cons1)
            add_homologs(homs1, control_homologue_species_run)
        else:
            cn_zeroes['All'] += 1
            cn_zeroes[ov_type] += 1
        if cons2 != 0:
            cn['All'].append(cons2)
            cn[ov_type].append(cons2)
            add_homologs(homs2, control_homologue_species_run)
        else:
            cn_zeroes['All'] += 1
            cn_zeroes[ov_type] += 1

        add_homologs(homs_max1, control_max_homologue_species_run)
        add_homologs(homs_max2, control_max_homologue_species_run)

    control_homologue_species_runs.append(control_homologue_species_run)
    control_max_homologue_species_runs.append(control_max_homologue_species_run)

    if len(control_groups_data) == 0: control_groups_data.append(
        f"Random group\tGC (All-{len(gc['All'])})\tGC (Divergent-{len(gc['Divergent'])})\tGC (Convergent-{len(gc['Convergent'])})\tGC (Nested d.-{len(gc['Nested (diff strand)'])})\tGC (Tandem-{len(gc['Tandem'])})\tGC (Nested s.-{len(gc['Nested (same strand)'])})\tCons-0 (All)\tCons-0 (Divergent)\tCons-0 (Convergent)\tCons-0 (Nested d.)\tCons-0 (Tandem)\tCons-0 (Nested s.)\tCons-!0 (All)\tCons-!0 (Divergent)\tCons-!0 (Convergent)\tCons-!0 (Nested d.)\tCons-!0 (Tandem)\tCons-!0 (Nested s.)")
    control_groups_data.append(
        f"{run + 1}\t{np.mean(gc['All'])}\t{np.mean(gc['Divergent'])}\t{np.mean(gc['Convergent'])}\t{np.mean(gc['Nested (diff strand)'])}\t{np.mean(gc['Tandem'])}\t{np.mean(gc['Nested (same strand)'])}"
        f"\t{cn_zeroes['All']}\t{cn_zeroes['Divergent']}\t{cn_zeroes['Convergent']}\t{cn_zeroes['Nested (diff strand)']}\t{cn_zeroes['Tandem']}\t{cn_zeroes['Nested (same strand)']}"
        f"\t{np.median(cn['All']) if len(cn['All']) != 0 else -1}\t{np.median(cn['Divergent'])}\t{np.median(cn['Convergent']) if len(cn['Convergent']) != 0 else -1}\t{np.median(cn['Nested (diff strand)']) if len(cn['Nested (diff strand)']) != 0 else -1}\t{np.median(cn['Tandem']) if len(cn['Tandem']) != 0 else -1}\t{np.median(cn['Nested (same strand)']) if len(cn['Nested (same strand)']) != 0 else -1}")
    GC_means.append(np.mean(gc['All']))
    CONS_medians.append(np.median(cn['All']))
# observed_groups_overall = []
conservation_2d_data = []
for record in OGs_records:
    ov_type = record.transcriptic_overlap_type.short_name()
    cons1, cons2 = record.get_transcript_conservation_score(
        record.transcript1.id), record.get_transcript_conservation_score(record.transcript2.id)

    g1_label, g2_label = record.gene1_symbol, record.gene2_symbol
    item_label = f'{g1_label} - {g2_label}' if record.get_record_conservation_score() > 0 else ''
    if record.get_record_length() >= RECORD_MIN_LENGTH and cons1 > 0 and cons2 > 0:
        cons1, cons2 = (cons2, cons1) if cons1 < cons2 else (cons1, cons2)
        conservation_2d_data.append((cons2, cons1, ov_type, item_label))

for record in OGs_random_records:
    if record.get_record_length() >= RECORD_MIN_LENGTH:
        GC_observed_X += record.get_record_GC_content()
        GC_observed_count += 1
    cons1, cons2 = record.get_transcript_conservation_score(
        record.transcript1.id), record.get_transcript_conservation_score(record.transcript2.id)
    if cons1 > 0:
        CONS_observed.append(cons1)
        add_homologs(genome.get_transcript_homologue_species(record.transcript1.id), observed_homologue_species)
    if cons2 > 0:
        CONS_observed.append(cons2)
        add_homologs(genome.get_transcript_homologue_species(record.transcript2.id), observed_homologue_species)

    add_homologs(genome.get_gene_max_conserved_homologue_species(record.gene1.id), observed_max_homologue_species)
    add_homologs(genome.get_gene_max_conserved_homologue_species(record.gene2.id), observed_max_homologue_species)

write_excel_OGs_by_CDS_data(file, genome, mouse_genome, OGs_records, OGs_records_mouse, OGs_random_records,
                            control_groups_data,
                            normalize_runs_and_counts(control_homologue_species_runs),
                            normalize_runs_and_counts([observed_homologue_species]),
                            normalize_runs_and_counts(control_max_homologue_species_runs),
                            normalize_runs_and_counts([observed_max_homologue_species])
                            )

build_2d_plot(conservation_2d_data,
              path=output_plot_path + "OGs by CDS (filtered) [2d conservation].png")

build_histogram(GC_means, GC_observed_X / GC_observed_count,
                "კონტროლ ჯგუფების საშუალო <br>GC შემცველობები", "გადაფარული გენების <br> საშუალო GC შემცველობა",
                0.48, 0.6, -200, -30, output_plot_path + "GC Means.png")
build_histogram(CONS_medians, np.median(CONS_observed),
                "კონტროლ ჯგუფების კონსერვაციის მედიანები", "გადაფარული გენების<br>კონსერვაციის მედიანა",
                2, 7, 200, -30, output_plot_path + "Conservation Medians.png")
