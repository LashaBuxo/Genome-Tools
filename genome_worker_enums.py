# region ENUM classes
import enum


class SPECIES(enum.Enum):
    Homo_sapiens = 0
    Mus_musculus = 1
    Rattus_norvegicus = 2
    Danio_rerio = 3
    Drosophila_melanogaster = 4
    Caenorhabditis_elegans = 5
    Saccharomyces_cerevisiae = 6

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
        elif str == "Caenorhabditis elegans":
            return SPECIES.Caenorhabditis_elegans
        elif str == "Saccharomyces cerevisiae":
            return SPECIES.Saccharomyces_cerevisiae
        else:
            assert False

    @staticmethod
    def from_taxonID(taxon_ID):
        if taxon_ID == "NCBITaxon:9606":
            return SPECIES.Homo_sapiens
        elif taxon_ID == "NCBITaxon:10090":
            return SPECIES.Mus_musculus
        elif taxon_ID == "NCBITaxon:10116":
            return SPECIES.Rattus_norvegicus
        elif taxon_ID == "NCBITaxon:7955":
            return SPECIES.Danio_rerio
        elif taxon_ID == "NCBITaxon:7227":
            return SPECIES.Drosophila_melanogaster
        elif taxon_ID == "NCBITaxon:6239":
            return SPECIES.Caenorhabditis_elegans
        elif taxon_ID == "NCBITaxon:559292":
            return SPECIES.Saccharomyces_cerevisiae
        else:
            assert False

    def to_taxonID(self):
        if self == SPECIES.Homo_sapiens:
            return "NCBITaxon:9606"
        elif self == SPECIES.Mus_musculus:
            return "NCBITaxon:10090"
        elif self == SPECIES.Rattus_norvegicus:
            return "NCBITaxon:10116"
        elif self == SPECIES.Danio_rerio:
            return "NCBITaxon:7955"
        elif self == SPECIES.Drosophila_melanogaster:
            return "NCBITaxon:7227"
        elif self == SPECIES.Caenorhabditis_elegans:
            return "NCBITaxon:6239"
        elif self == SPECIES.Saccharomyces_cerevisiae:
            return "NCBITaxon:559292"
        else:
            assert False

    def __str__(self):
        return str(self.name).replace('_', ' ')


class ANNOTATIONS(enum.Enum):
    NCBI = 1
    ENSEMBL = 2


class ANNOTATION_LOAD(enum.Enum):
    GENES = 1
    GENES_AND_TRANSCRIPTS = 2
    GENES_AND_TRANSCRIPTS_AND_CDS = 3
    GENES_AND_TRANSCRIPTS_AND_FRAGMENTS = 4


class SEQUENCE_LOAD(enum.Enum):
    LOAD = 1
    NOT_LOAD = 2
