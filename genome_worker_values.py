# region Paths
WORKING_DATABASES_DIRECTORY = './databases/'

ENSEMBL_ANNOTATIONS = ['./genome_data/Ensembl/Homo sapiens/106/Homo_sapiens.GRCh38.106.chr.gff3',
                       './genome_data/Ensembl/Mus musculus/106/Mus_musculus.GRCm39.106.chr.gff3',
                       './genome_data/Ensembl/Rattus norvegicus/106/Rattus_norvegicus.mRatBN7.2.106.gff3',
                       './genome_data/Ensembl/Danio rerio/106/Danio_rerio.GRCz11.106.chr.gff3',
                       './genome_data/Ensembl/Drosophila melanogaster/106/Drosophila_melanogaster.BDGP6.32.106.gff3',
                       './genome_data/Ensembl/Caenorhabditis elegans/106/Caenorhabditis_elegans.WBcel235.106.gff3',
                       './genome_data/Ensembl/Saccharomyces cerevisiae/106/Saccharomyces_cerevisiae.R64-1-1.106.gff3']

ENSEMBL_SEQUENCES = ['./genome_data/Ensembl/Homo sapiens/106/Homo_sapiens.GRCh38.dna.primary_assembly.fa',
                     './genome_data/Ensembl/Mus musculus/106/Mus_musculus.GRCm39.dna.primary_assembly.fa',
                     './genome_data/Ensembl/Rattus norvegicus/106/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa',
                     './genome_data/Ensembl/Danio rerio/106/Danio_rerio.GRCz11.dna.primary_assembly.fa',
                     './genome_data/Ensembl/Drosophila melanogaster/106/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa',
                     './genome_data/Ensembl/Caenorhabditis elegans/106/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa',
                     './genome_data/Ensembl/Saccharomyces cerevisiae/106/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa']

NCBI_ANNOTATIONS = ['./genome_data/NCBI/Homo sapiens/110/(+utrs)GCF_000001405.40_GRCh38.p14_genomic.gff', ]
NCBI_SEQUENCES = ['./genome_data/NCBI/Homo sapiens/110/GCF_000001405.40_GRCh38.p14_genomic.fna']
# endregion

# region Chromosome numbers

NUMBER_OF_CHROMOSOMES = [25, 23, 24, 26, 8, 7, 17]
#   Homo sapiens - 25: 22 autosome + 2 sex + 1 mitochondrion chromosome
#   Mus musculus - 23: 20 autosome + 2 sex + 1 mitochondrion chromosome
#   Rattus norvegicus - 24: 21 autosome + 2 sex + 1 mitochondrion chromosome
#   Danio rerio - 26: 25 autosome + 1 mitochondrion chromosome
#   Drosophila melanogaster - 8: 5 autosome + 2 sex + 1 mitochondrion chromosome
#   Caenorhabditis elegans - 7: 5 autosome + 1 sex + 1 mitochondrion chromosome
#   Saccharomyces cerevisiae - 17: 16 autosome + 1 mitochondrion chromosome
ENSEMBL_CHR_MAP_FOR_DROSOPHILA = {"2L": 1, "2R": 2, "3L": 3, "3R": 4, "4": 5, "X": 6, "Y": 7,
                                  "mitochondrion_genome": 8}

ENSEMBL_CHR_MAP_FOR_Caenorhabditis = {"I": 1, "II": 2, "III": 3, "IV": 4, "V": 5, "X": 6, "MtDNA": 7}

ENSEMBL_CHR_MAP_FOR_Saccharomyces = {"I": 1, "II": 2, "III": 3, "IV": 4, "V": 5, "VI": 6, "VII": 7, "VIII": 8, "IX": 9,
                                     "X": 10, "XI": 11, "XII": 12, "XIII": 13, "XIV": 14, "XV": 15, "XVI": 16,
                                     "Mito": 17}

# endregion
