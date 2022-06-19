from worker_genome import *
from worker_genome_enums import *
from overlap_by_CDS_record import OverlapRecord


class SEARCH_METHOD(IntEnum):
    BY_GENE_BOUNDARIES = 0
    BY_GENE_CDSs = 1
    BY_GENE_CDSs_EXCEPT_ATI = 2


# min length and record_compare_function can be used when CDS overlaps are collected
def overlaps_finder(genome: GenomeWorker, method: SEARCH_METHOD, min_length=1, record_compare_function=None):
    OGs_by_any_ISO_graph = AnalyzerGraph()
    OGs_records = []
    # region observed data ( OGs by CDS) collector
    for chr_id in range(1, genome.chromosomes_count() + 1):
        genes_cnt = genome.genes_count_on_chr(chr_id)

        for i in range(0, genes_cnt):
            for j in range(i + 1, genes_cnt):
                gene_A, gene_B = genome.gene_by_indexes(chr_id, i), genome.gene_by_indexes(chr_id, j)

                overlap_type = genome.get_overlap_type(gene_A, gene_B)

                # if genes don't overlap, there is no point in searching overlaps by fragments (CDS in this case)
                if overlap_type is OVERLAP_TYPE.NONE: continue

                # if we are searching genes overlap by boundaries we just make edge and continue searching
                if method == SEARCH_METHOD.BY_GENE_BOUNDARIES:
                    edge = AnalyzerGraph.GraphEdge(gene_A.id, gene_B.id)
                    OGs_by_any_ISO_graph.add_edge(edge)
                    continue

                genes_overlapped_except_ATI = False
                occurred_in_phase_overlap = False

                transcripts_A = genome.get_transcripts_from_gene(gene_A.id)
                transcripts_B = genome.get_transcripts_from_gene(gene_B.id)

                best_record = None
                for transcript1 in transcripts_A:
                    for transcript2 in transcripts_B:
                        if transcript1.strand == '+':
                            transcript_a = transcript1
                            transcript_b = transcript2
                        else:
                            transcript_a = transcript2
                            transcript_b = transcript1
                        non_ati_segments, ov_l, ATI = genome.get_overlaps_between_transcripts(transcript_a.id,
                                                                                              transcript_b.id)
                        # if overlaps contains in_phase overlaps
                        if ATI: occurred_in_phase_overlap = True

                        # isoforms overlapped by at least one CDS without same-strand in-phase manner
                        if len(non_ati_segments) > 0:
                            record = OverlapRecord(transcript_a, transcript_b, chr_id, genome, non_ati_segments)

                            # if overlap is below our desired number, then we should not connect genes
                            if record.get_record_length() < min_length:
                                continue
                            genes_overlapped_except_ATI = True

                            # if there is no comparison function, then we must include all records
                            if record_compare_function is None:
                                OGs_records.append(record)
                            else:
                                best_record = record_compare_function(best_record, record)

                # isoforms overlapped by at least one CDS without same-strand in-phase manner
                if genes_overlapped_except_ATI:
                    if record_compare_function is not None:  # if there is comparison function, then best is added
                        OGs_records.append(best_record)
                    edge = AnalyzerGraph.GraphEdge(gene_A.id, gene_B.id, OVERLAP_INTERACTION.ANY_EXCEPT_ATI)
                    OGs_by_any_ISO_graph.add_edge(edge)
                else:
                    if method == SEARCH_METHOD.BY_GENE_CDSs_EXCEPT_ATI: continue
                    if occurred_in_phase_overlap:  # all overlaps were same strand in-phase
                        edge = AnalyzerGraph.GraphEdge(gene_A.id, gene_B.id, OVERLAP_INTERACTION.ATI)
                        OGs_by_any_ISO_graph.add_edge(edge)
    return OGs_by_any_ISO_graph, OGs_records
