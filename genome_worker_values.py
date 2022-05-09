# region Paths
WORKING_DATABASES_DIRECTORY = './databases/'

EXPRESSION_DATA_PATH = './expression_data/transcripts_all_ens75.tsv'


def EXPRESSION_DOWNLOAD_URL(gene_id: str):
    return f'https://www.proteinatlas.org/{gene_id}.json'


ENSEMBL_ANNOTATIONS = ['./genome_data/Ensembl/Homo sapiens/106/Homo_sapiens.GRCh38.106.chr.gff3',
                       './genome_data/Ensembl/Mus musculus/106/Mus_musculus.GRCm39.106.chr.gff3',
                       './genome_data/Ensembl/Rattus norvegicus/106/Rattus_norvegicus.mRatBN7.2.106.gff3',
                       './genome_data/Ensembl/Danio rerio/106/Danio_rerio.GRCz11.106.chr.gff3',
                       './genome_data/Ensembl/Drosophila melanogaster/106/Drosophila_melanogaster.BDGP6.32.106.gff3']

ENSEMBL_SEQUENCES = ['./genome_data/Ensembl/Homo sapiens/106/Homo_sapiens.GRCh38.dna.primary_assembly.fa',
                     './genome_data/Ensembl/Mus musculus/106/Mus_musculus.GRCm39.dna.primary_assembly.fa',
                     './genome_data/Ensembl/Rattus norvegicus/106/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa',
                     './genome_data/Ensembl/Danio rerio/106/Danio_rerio.GRCz11.dna.primary_assembly.fa',
                     './genome_data/Ensembl/Drosophila melanogaster/106/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa',
                     ]

NCBI_ANNOTATIONS = ['./genome_data/NCBI/Homo sapiens/110/(+utrs)GCF_000001405.40_GRCh38.p14_genomic.gff',
                    './genome_data/NCBI/Mus musculus/(+utrs)GCF_000001635.27_GRCm39_genomic.gff',
                    './genome_data/NCBI/Rattus norvegicus/(+utrs)GCF_015227675.2_mRatBN7.2_genomic.gff',
                    './genome_data/NCBI/Danio rerio/(+utrs)GCF_000002035.6_GRCz11_genomic.gff',
                    './genome_data/NCBI/Drosophila melanogaster/(+utrs)GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff']

NCBI_SEQUENCES = ['./genome_data/NCBI/Homo sapiens/110/GCF_000001405.40_GRCh38.p14_genomic.fna',
                  './genome_data/NCBI/Mus musculus/GCF_000001635.27_GRCm39_genomic.fna',
                  './genome_data/NCBI/Rattus norvegicus/GCF_015227675.2_mRatBN7.2_genomic.fna',
                  './genome_data/NCBI/Danio rerio/GCF_000002035.6_GRCz11_genomic.fna',
                  './genome_data/NCBI/Drosophila melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna']
# endregion

# region Chromosome numbers

NUMBER_OF_CHROMOSOMES = [24, 21, 22, 25, 6, 6]
#   Homo sapiens - 25: 22 autosome + 2 sex + (1 mitochondrion chromosome IGNORED)
#   Mus musculus - 22: 19 autosome + 2 sex + (1 mitochondrion chromosome IGNORED)
#   Rattus norvegicus - 23: 20 autosome + 2 sex + (1 mitochondrion chromosome IGNORED)
#   Danio rerio - 26: 25 autosome + (1 mitochondrion chromosome
#   Drosophila melanogaster - 7: 5 autosome + 2 sex + (1 mitochondrion chromosome IGNORED) - 1 Y
#                     NOTE: excluding 1 Y chromosome, coz NCBI does not have one
#   Caenorhabditis elegans - 7: 5 autosome + 1 sex + (1 mitochondrion chromosome IGNORED)
ENSEMBL_CHR_MAP_FOR_DROSOPHILA = {"2L": 1, "2R": 2, "3L": 3, "3R": 4, "4": 5, "X": 6,
                                  "mitochondrion_genome": -1}

ENSEMBL_CHR_MAP_FOR_Caenorhabditis = {"I": 1, "II": 2, "III": 3, "IV": 4, "V": 5, "X": 6, "MtDNA": -1}

# endregion
