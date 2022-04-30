# region ENUM classes
import enum


class TRANSCRIPT_CRITERIA(enum.Enum):
    NONE = 0
    LONGEST = 1
    LONGEST_CDS = 2
    LONGEST_CDS_AND_UTRs = 3  # similar to exons
    RANDOM = 4


class ANNOTATIONS(enum.Enum):
    NCBI = 1
    ENSEMBL = 2

    def short_name(self):
        if self == ANNOTATIONS.ENSEMBL:
            return 'Ensembl'
        else:
            return 'NCBI'

    def full_name(self):
        if self == ANNOTATIONS.ENSEMBL:
            return 'Ensembl'
        else:
            return 'RefSeq (NCBI)'


class ANNOTATION_LOAD(enum.Enum):
    GENES = 1
    GENES_AND_TRANSCRIPTS = 2
    GENES_AND_TRANSCRIPTS_AND_CDS = 3
    GENES_AND_TRANSCRIPTS_AND_FRAGMENTS = 4


class SEQUENCE_LOAD(enum.Enum):
    LOAD = 1
    NOT_LOAD = 2


class SPECIES(enum.Enum):
    Homo_sapiens = 0
    Mus_musculus = 1
    Rattus_norvegicus = 2
    Danio_rerio = 3
    Drosophila_melanogaster = 4

    @staticmethod
    def from_string(str):
        if str == "Homo sapiens":
            return SPECIES.Homo_sapiens
        elif str == "Mus musculus":
            return SPECIES.Mus_musculus
        elif str == "Rattus norvegicus":
            return SPECIES.Rattus_norvegicus
        elif str == "Danio rerio":
            return SPECIES.Danio_rerio
        elif str == "Drosophila melanogaster":
            return SPECIES.Drosophila_melanogaster
        else:
            if not str.__contains__('_'): assert False
            return SPECIES.from_string(str.replace('_', ' '))

    # @staticmethod
    # def from_taxonID(taxon_ID):
    #     if taxon_ID == "NCBITaxon:9606":
    #         return SPECIES.Homo_sapiens
    #     elif taxon_ID == "NCBITaxon:10090":
    #         return SPECIES.Mus_musculus
    #     elif taxon_ID == "NCBITaxon:10116":
    #         return SPECIES.Rattus_norvegicus
    #     elif taxon_ID == "NCBITaxon:7955":
    #         return SPECIES.Danio_rerio
    #     elif taxon_ID == "NCBITaxon:7227":
    #         return SPECIES.Drosophila_melanogaster
    #     else:
    #         assert False

    def short_name(self):
        if self == SPECIES.Homo_sapiens:
            return "Human"
        elif self == SPECIES.Mus_musculus:
            return "Mouse"
        elif self == SPECIES.Rattus_norvegicus:
            return "Rat"
        elif self == SPECIES.Danio_rerio:
            return "ZebraFish"
        elif self == SPECIES.Drosophila_melanogaster:
            return "Drosophila"
        else:
            assert False

    #
    # def to_taxonID(self):
    #     if self == SPECIES.Homo_sapiens:
    #         return "NCBITaxon:9606"
    #     elif self == SPECIES.Mus_musculus:
    #         return "NCBITaxon:10090"
    #     elif self == SPECIES.Rattus_norvegicus:
    #         return "NCBITaxon:10116"
    #     elif self == SPECIES.Danio_rerio:
    #         return "NCBITaxon:7955"
    #     elif self == SPECIES.Drosophila_melanogaster:
    #         return "NCBITaxon:7227"
    #     else:
    #         assert False

    def __str__(self):
        return str(self.name).replace('_', ' ')
