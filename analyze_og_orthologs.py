# Author: "Lasha Bukhnikashvili"
#
# Description:
#   ***
#
# Usage:
#   analyze_og.py <species1> <species2> <annotation>
#
# Params (possible) to run:
#   species1/species2:  Homo sapiens / Rattus norvegicus / Mus musculus / Danio rerio / Drosophila melanogaster,
#               Caenorhabditis elegans
#   annotation: NCBI / Ensembl
#
# Example:
#   python analyze_og_orthologs.py 'Homo sapiens' 'Mus musculus' NCBI
#
# Output:
#   prints stats in console

import sys
import time

from genome_worker import *

id_prefixes_ortho_to_ncbi = {'HGNC': 'HGNC',
                             'FB': 'FLYBASE',
                             'ZFIN': 'ZFIN',
                             'RGD': 'RGD',
                             'WB': 'WormBase',
                             'MGI': 'MGI'}

orthologs_path = "./genome_data/ORTHOLOGY-ALLIANCE_COMBINED.tsv"

symbol_storage = {}
ortholog_for_species1 = {}
ortholog_for_species2 = {}


def store_symbol(gene_id, gene_symbol):
    if symbol_storage.__contains__(gene_id):
        assert symbol_storage[gene_id] == gene_symbol
        return
    symbol_storage[gene_id] = gene_symbol


def store_ortholog_feature(target_species1, target_species2,
                           gene1_id, gene1_symbol, gene1_species_taxon_id, gene1_species_name,
                           gene2_id, gene2_symbol, gene2_species_taxon_id, gene2_species_name):
    if target_species1.to_taxonID() != gene1_species_taxon_id: return
    if target_species2.to_taxonID() != gene2_species_taxon_id: return

    assert id_prefixes_ortho_to_ncbi.keys().__contains__(gene1_id.split(':')[0])
    assert id_prefixes_ortho_to_ncbi.keys().__contains__(gene2_id.split(':')[0])

    # WB:WBGene00011502 -> WormBase:WBGene00011502
    gene1_id = id_prefixes_ortho_to_ncbi[gene1_id.split(':')[0]] + ":" + gene1_id.split(':')[1]
    gene2_id = id_prefixes_ortho_to_ncbi[gene2_id.split(':')[0]] + ":" + gene2_id.split(':')[1]

    store_symbol(gene1_id, gene1_symbol)
    store_symbol(gene2_id, gene2_symbol)

    orthologs_for_gene1 = [] if not ortholog_for_species1.__contains__(gene1_id) else ortholog_for_species1[gene1_id]
    orthologs_for_gene2 = [] if not ortholog_for_species2.__contains__(gene2_id) else ortholog_for_species2[gene2_id]

    assert not orthologs_for_gene1.__contains__(gene2_id)
    assert not orthologs_for_gene2.__contains__(gene1_id)

    orthologs_for_gene1.append(gene2_id)
    orthologs_for_gene2.append(gene1_id)

    ortholog_for_species1[gene1_id] = orthologs_for_gene1
    ortholog_for_species2[gene2_id] = orthologs_for_gene2


def load_orthologs(target_species1, target_species2, ):
    lines = open(orthologs_path, 'r').readlines()

    past_description = False
    for line in lines:
        if past_description:
            values = line.split("\t")
            store_ortholog_feature(target_species1, target_species2,
                                   values[0], values[1], values[2], values[3],
                                   values[4], values[5], values[6], values[7])
        if line.startswith("Gene1ID"): past_description = True


start_time = time.time()

assert len(sys.argv) == 4
assert sys.argv[1] != sys.argv[2]
assert sys.argv[3] == 'NCBI' or sys.argv[3] == 'Ensembl'

species1 = SPECIES.from_string(sys.argv[1])
species2 = SPECIES.from_string(sys.argv[2])
annotation_source = ANNOTATIONS.NCBI if sys.argv[3] == 'NCBI' else ANNOTATIONS.ENSEMBL

load_orthologs(species1, species2)

genome1 = GenomeWorker(species1, annotation_source, ANNOTATION_LOAD.GENES, SEQUENCE_LOAD.NOT_LOAD)
genome2 = GenomeWorker(species2, annotation_source, ANNOTATION_LOAD.GENES, SEQUENCE_LOAD.NOT_LOAD)

genes_with_correct_prefixes = 0
genes_with_correct_prefixes_and_homologs = 0
possible_prefixes = list(id_prefixes_ortho_to_ncbi.values())

####################################################
human_OGs_by_CDS_path = "./results/OGs_Homo_sapiens (+descriptions) - NCBI [diff_stranded, both_ORF].txt"
lines = open(human_OGs_by_CDS_path, 'r').readlines()
OGPs = []
OG_corrected_id = {}
for line in lines:
    if not line.__contains__('genes: ['): continue
    part = line.split('genes: [')[1].replace(']', '')
    gene1, gene2 = part.split(';')[0].rstrip('\n'), part.split(';')[1].rstrip('\n')

    OGPs.append((gene1, gene2))
    OG_corrected_id[gene1] = 'xxx'
    OG_corrected_id[gene2] = 'xxx'
    # print(gene1 + " " + gene2)
print(OG_corrected_id.keys())
####################################################

# excluding mitochondria
for chr_id in range(1, genome1.chromosomes_count()):
    genes_cnt = genome1.genes_count_on_chr(chr_id)
    for i in range(0, genes_cnt):
        gene = genome1.gene_by_ind(chr_id, i)
        assert gene.attributes.__contains__('Dbxref')
        ids = gene.attributes['Dbxref']

        contains_valid_id = False
        corrected_id = ""

        for id in ids:
            parts = id.split(':')
            if possible_prefixes.__contains__(parts[0]):
                corrected_id = parts[0] + ":" + parts[len(parts) - 1]

        if corrected_id != "":
            genes_with_correct_prefixes += 1
            if ortholog_for_species1.__contains__(corrected_id):
                genes_with_correct_prefixes_and_homologs += 1
            if OG_corrected_id.__contains__(gene.id):
                OG_corrected_id[gene.id] = corrected_id

print(str(genes_with_correct_prefixes_and_homologs) + "/" + str(genes_with_correct_prefixes))

####################################################

OGs = 0
OGs_with_homolog = 0
for OGP in OGPs:
    # print(str(OGP[0]) + " " + str(OGP[1]))
    corr_id1 = '' if not OG_corrected_id.__contains__(OGP[0]) else OG_corrected_id[OGP[0]]
    corr_id2 = '' if not OG_corrected_id.__contains__(OGP[1]) else OG_corrected_id[OGP[1]]
    if corr_id1 == 'xxx' or corr_id2 == 'xxx': continue
    OGs += 2
    if ortholog_for_species1.__contains__(corr_id1):
        print(ortholog_for_species1[corr_id1])
        OGs_with_homolog += 1
    if ortholog_for_species1.__contains__(corr_id2):
        print(ortholog_for_species1[corr_id2])
        OGs_with_homolog += 1
    print(" ")
print(str(OGs_with_homolog) + "/" + str(OGs))

####################################################

print("--- %s seconds ---" % (time.time() - start_time))
