import statistics

import numpy

from worker_genome import *
from GO_graph import *
import scipy.stats

from AgeData import *


def ShinyGOLine2Data(line: str):
    arr = line.replace('\n', '').split(' "," ')
    if len(arr) == 1: return None
    assert len(arr) == 2
    info_data, site_and_genes = arr[0], arr[1]

    detailed_arr = info_data.split(',')
    FDR, nGenes, pathway_genes, fold = detailed_arr[0], detailed_arr[1], detailed_arr[2], detailed_arr[3]

    cnt = 0
    ind = -1
    for i in range(len(info_data)):
        if info_data[i] == ',':
            cnt += 1
            if cnt == 4:
                ind = i
                break
    assert ind != -1
    go_name = info_data[ind + 2:len(info_data)]
    ########################
    site_str, genes_str = site_and_genes.split(',')[0], site_and_genes.split(',')[1]
    go_id = site_str[len(site_str) - 11:len(site_str) - 1]

    genes_str = genes_str.replace('"', '').split(' ')
    enriched_gene_symbols = []
    for gene_str in genes_str:
        if len(gene_str) > 0:
            enriched_gene_symbols.append(gene_str)

    return float(FDR), int(nGenes), int(pathway_genes), float(fold), go_name, go_id, enriched_gene_symbols


def filter_out_children_terms(terms, graph: GOGraph):
    filtered_list = []

    not_found_terms = {}
    for i in range(0, len(terms)):
        flag = False

        for j in range(0, len(terms)):
            if j == i: continue
            if not graph.is_term_found(terms[i]):
                not_found_terms[terms[i]] = True
            elif not graph.is_term_found(terms[j]):
                not_found_terms[terms[j]] = True
            elif graph.is_first_term_child_of_second(terms[i], terms[j]):
                flag = True
                break
        if not flag:
            filtered_list.append(terms[i])
    print(f'From {len(terms)} terms {len(not_found_terms)} not found in ontology graph!')
    return filtered_list


def sort_value(term, pathways_enriched, graph: GOGraph):
    enriched_species_dict = pathways_enriched[term]

    namespace = graph.get_term_namespace(term) if graph.is_term_found(term) else "biological_process"
    namespace_value = 0 if namespace == "biological_process" else (
        1000 if namespace == "cellular_component" else 1000000)

    species_cnt_occurred = 0
    p_values = []
    fold_values = []
    nGenes_values = []
    for species in specie_list:
        if enriched_species_dict.__contains__(species):
            species_cnt_occurred += 1
            (FDR, nGenes, Pathway_Genes, Fold, _) = enriched_species_dict[species]
            p_values.append(FDR)
            fold_values.append(Fold)
            nGenes_values.append(nGenes)
    mean_p_value = statistics.mean(p_values)
    mean_fold_value = statistics.mean(fold_values)
    mean_nGenes_value = statistics.mean(nGenes_values)
    # namespace is greater, cnt is greater, mean_p_smaller

    return -namespace_value + species_cnt_occurred - mean_p_value


specie_list = [SPECIES.Homo_sapiens,
               SPECIES.Mus_musculus,
               SPECIES.Drosophila_melanogaster,
               SPECIES.Arabidopsis_thaliana
               ]

cut_off = 0.01
min_number_of_organisms = 3

GO_graph = GOGraph()

gene_set_types = ["convergent", "nested (diff)", "divergent",
                   "OGs"]

ontology_types = ["enrichment_all.csv", "enrichment_all (1).csv", "enrichment_all (2).csv"]

pathway_id2_pathway = {}

for gene_set_type in gene_set_types:
    out_file = open(f"generated_data/Enriched_GOs ({gene_set_type}).txt", "w")
    pathways_enriched = {}
    for species in specie_list:
        for ontology_type in ontology_types:
            file_path = f"./used_data/ShinyGO/{species.short_name()}/{gene_set_type}/{ontology_type}"
            file = open(file_path, "r")
            lines = file.readlines()
            for index in range(1, len(lines)):
                data = ShinyGOLine2Data(lines[index])
                if data is None: continue
                FDR, nGenes, pathway_Genes, fold, pathway, pathway_id, gene_set, = data

                if FDR > cut_off: continue
                pathway_id2_pathway[pathway_id] = pathway

                if not pathways_enriched.__contains__(pathway_id):
                    pathways_enriched[pathway_id] = {}
                pathways_enriched[pathway_id][species] = (FDR, nGenes, pathway_Genes, fold, gene_set)

            file.close()

    valid_pathways = []
    for pathway_id, enriched_species_dict in pathways_enriched.items():
        if len(enriched_species_dict) >= min_number_of_organisms:
            valid_pathways.append(pathway_id)

    valid_pathways = filter_out_children_terms(valid_pathways, graph=GO_graph)

    valid_pathways.sort(key=lambda x: sort_value(x, pathways_enriched, graph=GO_graph), reverse=True)

    for pathway_id in valid_pathways:
        enriched_species_dict = pathways_enriched[pathway_id]
        assert len(enriched_species_dict) >= min_number_of_organisms
        pathway = pathway_id2_pathway[pathway_id]
        out_file.write(f"{pathway}\t")
        for species in specie_list:
            if not enriched_species_dict.__contains__(species):
                out_file.write("-\t-\t-\t-\t")
            else:
                (FDR, nGenes, Pathway_Genes, Fold, _) = enriched_species_dict[species]
                out_file.write(f"{FDR}\t{nGenes}\t{Pathway_Genes}\t{Fold}\t")
        out_file.write(f"\n")

    out_file.close()

    if gene_set_type == "OGs":
        out_file = open("generated_data/human_pathway_details.txt", "w")

        genome = GenomeWorker(SPECIES.Homo_sapiens, ANNOTATION_LOAD.GENES, SEQUENCE_LOAD.NOT_LOAD)
        age_data = AgeData(SPECIES.Homo_sapiens, genome)

        for pathway_id in valid_pathways:
            enriched_species_dict = pathways_enriched[pathway_id]
            if not enriched_species_dict.__contains__(SPECIES.Homo_sapiens): continue
            (FDR, nGenes, Pathway_Genes, Fold, gene_set) = enriched_species_dict[SPECIES.Homo_sapiens]

            pathway = pathway_id2_pathway[pathway_id]
            out_file.write(f"{pathway}\t{pathway_id}\t{nGenes}\t{Pathway_Genes}\t{Fold}\t")

            not_found_genes = 0
            length = []
            ages = []
            for gene_sym in gene_set:
                gene = genome.gene_by_symbol(gene_sym)
                if gene is None:
                    not_found_genes += 1
                    continue
                length.append(gene.end - gene.start + 1)
                age = age_data.get_gene_age(gene.id)
                if age is not None:
                    ages.append(age)
            out_file.write(f"{not_found_genes}\t{statistics.median(length)}\t{statistics.median(ages)}\n")
        out_file.close()
