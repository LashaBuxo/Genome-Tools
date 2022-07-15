import requests, sys
from worker_genome import *
# now we have two json STRINGS
import json

s="'GO:0003677', 'GO:0003950', 'GO:0005634', 'GO:0006281', 'GO:0006366', 'GO:0006471', 'GO:0003950', 'GO:0005654', 'GO:0005654', 'GO:0005654', 'GO:0005654', 'GO:0005654', 'GO:0005654', 'GO:0005654', 'GO:0005654', 'GO:0005654', 'GO:0005654', 'GO:0005654', 'GO:0005654', 'GO:0005654', 'GO:0005654', 'GO:0005654', 'GO:0005654', 'GO:0005654', 'GO:0005654', 'GO:0005654', 'GO:0005654', 'GO:0005654', 'GO:0005654', 'GO:0005654', 'GO:0005654', 'GO:0005654', 'GO:0000122', 'GO:0000723', 'GO:0030225', 'GO:0006281'"
arr=s.split(', ')
z=""
for x in arr:
    s=x.replace("'","")
    if z=="": z=s
    else:
        z+=f',{s}'
print(z)
assert False
class AnnotationRecordKey:
    def __init__(self, gene_symbol, taxon_ID, taxon_Name):
        self.gene_symbol = gene_symbol
        self.taxon_ID = taxon_ID
        self.taxon_Name = taxon_Name

    def __eq__(self, another):
        return self.gene_symbol + self.taxon_Name == another.gene_symbol + another.taxon_Name

    def __hash__(self):
        return hash(self.gene_symbol + self.taxon_Name)


class AnnotationRecordValue:
    def __init__(self):
        self.go_IDs = []
        self.go_Names = []
        self.go_Qualifiers = []
        self.go_EvidenceCodes = []

    def add_value(self, id, name, qualifier, evidence_code):
        self.go_IDs.append(id)
        self.go_Names.append(name)
        self.go_Qualifiers.append(qualifier)
        self.go_EvidenceCodes.append(evidence_code)


def DownloadAnnotations(taxonID, entries, page):
    # ECO:0000304 author statement supported by traceable reference used in manual assertion
    # "reference=PMID:21873635"

    slim_generic_ids = "GO:1901135,GO:0140053,GO:0140014,GO:0140013,GO:0098754,GO:0098542,GO:0072659,GO:0071941,GO:0071554,GO:0065003,GO:0061024,GO:0061007,GO:0055086,GO:0055085,GO:0055065,GO:0051604,GO:0050886,GO:0050877,GO:0048870,GO:0048856,GO:0044782,GO:0042254,GO:0042060,GO:0036211,GO:0034330,GO:0032200,GO:0031047,GO:0030198,GO:0030163,GO:0030154,GO:0023052,GO:0022600,GO:0022414,GO:0016192,GO:0016073,GO:0016071,GO:0015979,GO:0012501,GO:0007568,GO:0007163,GO:0007155,GO:0007059,GO:0007040,GO:0007031,GO:0007018,GO:0007010,GO:0007005,GO:0006954,GO:0006914,GO:0006913,GO:0006886,GO:0006790,GO:0006766,GO:0006629,GO:0006575,GO:0006520,GO:0006486,GO:0006457,GO:0006399,GO:0006355,GO:0006351,GO:0006325,GO:0006310,GO:0006260,GO:0006091,GO:0005975,GO:0003016,GO:0003014,GO:0003013,GO:0003012,GO:0002376,GO:0002181,GO:0000910,GO:0000278,GO:0140691,GO:0140657,GO:0140313,GO:0140299,GO:0140223,GO:0140110,GO:0140104,GO:0140098,GO:0140097,GO:0140096,GO:0120274,GO:0098772,GO:0098631,GO:0090729,GO:0060090,GO:0060089,GO:0048018,GO:0045735,GO:0045182,GO:0044183,GO:0042393,GO:0038024,GO:0031386,GO:0016874,GO:0016853,GO:0016829,GO:0016787,GO:0016740,GO:0016491,GO:0016209,GO:0009975,GO:0008289,GO:0008092,GO:0005215,GO:0005198,GO:0003924,GO:0003824,GO:0003774,GO:0003723,GO:0003677,GO:0001618,GO:0043226,GO:0031410,GO:0031012,GO:0030312,GO:0009579,GO:0009536,GO:0005929,GO:0005886,GO:0005856,GO:0005840,GO:0005829,GO:0005815,GO:0005811,GO:0005794,GO:0005783,GO:0005777,GO:0005773,GO:0005768,GO:0005764,GO:0005730,GO:0005694,GO:0005654,GO:0005635,GO:0005634,GO:0005618,GO:0005615,GO:0005576,GO:0000228,GO:0005739,GO:0006281"
    mito_DNA_ids = "GO:0006281,GO:0005739"
    # f"goId={mito_DNA_ids}&goUsageRelationships=is_a%2Cpart_of%2Coccurs_in&" \
    requestURL = "https://www.ebi.ac.uk/QuickGO/services/annotation/search?" \
                 "includeFields=goName&includeFields=taxonName&" \
                 f"taxonId={taxonID}&" \
                 "evidenceCode=ECO%3A0000304&" \
                 f"limit={entries}&page={page}"

    r = requests.get(requestURL, headers={"Accept": "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    responseBody = r.text

    dictData = json.loads(responseBody)
    return dictData


#
# file = open("used_data/GOAnnotations_generic.txt", 'w')
# init_data = DownloadAnnotations(SPECIES.Homo_sapiens.taxon_ID(), 200, 11)
# file.write(json.dumps(init_data,indent=2))

species_list = [SPECIES.Homo_sapiens
                # , SPECIES.Pan_troglodytes, SPECIES.Pongo_abelii,
                #             SPECIES.Papio_anubis, SPECIES.Nomascus_leucogenys, SPECIES.Callithrix_jacchus,
                #             SPECIES.Rattus_norvegicus, SPECIES.Mus_musculus
                ]

GO_annotations = {}

page_max_size = 200

for species in species_list:
    init_data = DownloadAnnotations(species.taxon_ID(), page_max_size, 1)
    page_count = init_data["pageInfo"]['total']

    if page_count == 0:
        print(f"no data for {species.name.replace('_', ' ')}")
        continue
    else:
        print(f"downloading GO annotations for {species.name.replace('_', ' ')}")

    for page_index in range(1, page_count + 1):
        print(f'\tpages downloaded {page_index}/{page_count}')
        retrieved_data = DownloadAnnotations(species.taxon_ID(), page_max_size, page_index)
        for result in retrieved_data["results"]:
            record_key = AnnotationRecordKey(result['symbol'], result['taxonId'], result['taxonName'])

            if not GO_annotations.__contains__(record_key):
                GO_annotations[record_key] = AnnotationRecordValue()
            GO_annotations[record_key].add_value(id=result['goId'],
                                                 name=result['goName'],
                                                 qualifier=result['qualifier'],
                                                 evidence_code=result['evidenceCode'])

file = open("used_data/GO_data/human_GO_annotations/GOAnnotations_all.txt", 'w')

for key, value in GO_annotations.items():
    file.write(f'{key.gene_symbol}\t{key.taxon_ID}\t{key.taxon_Name}\t'
               f'{value.go_IDs}\t{value.go_Names}\t{value.go_Qualifiers}\t{value.go_EvidenceCodes}\n')
