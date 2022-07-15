# region ENUM classes
import enum
from enum import IntEnum


class OVERLAP_INTERACTION(enum.Enum):
    # General classification of overlaps by CDS
    ANY_EXCEPT_ATI = 0  # any type of overlap, except same strand in-phase overlap (ATI)
    ATI = 1  # only same strand in-phase overlap

    @staticmethod
    def get_overlap_interactions():
        return [OVERLAP_INTERACTION.ANY_EXCEPT_ATI, OVERLAP_INTERACTION.ATI]


class OVERLAP_TYPE(IntEnum):
    NONE = 0
    # Specific classification of overlaps
    CONVERGENT = 1  # different strand + convergent
    DIFF_NESTED = 2  # different strand + nested
    DIVERGENT = 3  # different strand + divergent
    SAME_NESTED = 4  # one of the gene entirely located into others boundaries on same strand
    TANDEM = 5  # genes located on same strand and neither is nested
    MULTI = 6

    def short_name(self):
        if self == OVERLAP_TYPE.CONVERGENT:
            return "Convergent"
        elif self == OVERLAP_TYPE.DIFF_NESTED:
            return "Nested (diff strand)"
        elif self == OVERLAP_TYPE.DIVERGENT:
            return "Divergent"
        elif self == OVERLAP_TYPE.SAME_NESTED:
            return "Nested (same strand)"
        elif self == OVERLAP_TYPE.TANDEM:
            return "Tandem"
        elif self == OVERLAP_TYPE.MULTI:
            return "MULTI TYPE"
        else:
            assert False

    @staticmethod
    def get_overlap_types():
        return [OVERLAP_TYPE.CONVERGENT, OVERLAP_TYPE.DIFF_NESTED, OVERLAP_TYPE.DIVERGENT,
                OVERLAP_TYPE.SAME_NESTED, OVERLAP_TYPE.TANDEM, OVERLAP_TYPE.MULTI]


class ANNOTATION_LOAD(enum.Enum):
    GENES = 1
    GENES_AND_TRANSCRIPTS = 2
    GENES_AND_TRANSCRIPTS_AND_CDS = 3
    GENES_AND_TRANSCRIPTS_AND_FRAGMENTS = 4


class SEQUENCE_LOAD(enum.Enum):
    LOAD = 1
    NOT_LOAD = 2


class EXPRESSIONS_LOAD(enum.Enum):
    LOAD = 1
    NOT_LOAD = 2


class SPECIES(enum.Enum):
    Homo_sapiens = 0
    Pan_troglodytes = 1
    Pongo_abelii = 2
    Papio_anubis = 3
    Nomascus_leucogenys = 4
    Callithrix_jacchus = 5
    Rattus_norvegicus = 6
    Mus_musculus = 7
    Sus_scrofa = 8
    Canis_lupus_familiaris = 9
    Gallus_gallus = 10
    Danio_rerio = 11
    Drosophila_melanogaster = 12
    E_coli = 13
    Arabidopsis_thaliana = 14
    Saccharomyces_cerevisiae = 15

    @staticmethod
    def from_string(str):
        if str == "Homo sapiens":
            return SPECIES.Homo_sapiens
        elif str == "Pan troglodytes":
            return SPECIES.Pan_troglodytes
        elif str == "Pongo abelii":
            return SPECIES.Pongo_abelii
        elif str == "Papio anubis":
            return SPECIES.Papio_anubis
        elif str == "Nomascus leucogenys":
            return SPECIES.Nomascus_leucogenys
        elif str == "Callithrix jacchus":
            return SPECIES.Callithrix_jacchus
        elif str == "Rattus norvegicus":
            return SPECIES.Rattus_norvegicus
        elif str == "Mus musculus":
            return SPECIES.Mus_musculus
        elif str == "Sus scrofa":
            return SPECIES.Sus_scrofa
        elif str == "Canis lupus familiaris":
            return SPECIES.Canis_lupus_familiaris
        elif str == "Gallus gallus":
            return SPECIES.Gallus_gallus
        elif str == "Danio rerio":
            return SPECIES.Danio_rerio
        elif str == "Drosophila melanogaster":
            return SPECIES.Drosophila_melanogaster
        elif str == "E coli":
            return SPECIES.E_coli
        elif str == "Arabidopsis thaliana":
            return SPECIES.Arabidopsis_thaliana
        elif str == "Saccharomyces cerevisiae":
            return SPECIES.Saccharomyces_cerevisiae
        else:
            if not str.__contains__('_'): assert False
            return SPECIES.from_string(str.replace('_', ' '))

    def shortest_name(self):
        return self.short_name()

    def taxon_ID(self):
        if self == SPECIES.Homo_sapiens:
            return "9606"
        elif self == SPECIES.Pan_troglodytes:
            return "9598"
        elif self == SPECIES.Pongo_abelii:
            return "9601"
        elif self == SPECIES.Papio_anubis:
            return "9555"
        elif self == SPECIES.Nomascus_leucogenys:
            return "61853"
        elif self == SPECIES.Callithrix_jacchus:
            return "9483"
        elif self == SPECIES.Rattus_norvegicus:
            return "10116"
        elif self == SPECIES.Mus_musculus:
            return "10090"
        elif self == SPECIES.Sus_scrofa:
            return "9823"
        elif self == SPECIES.Canis_lupus_familiaris:
            return "9615"
        elif self == SPECIES.Gallus_gallus:
            return "9031"
        elif self == SPECIES.Danio_rerio:
            return "7955"
        elif self == SPECIES.Drosophila_melanogaster:
            return "7227"
        elif self == SPECIES.E_coli:
            return "511145"
        elif self == SPECIES.Arabidopsis_thaliana:
            return "3702"
        elif self == SPECIES.Saccharomyces_cerevisiae:
            return "4932"
        else:
            assert False

    def short_name(self):
        if self == SPECIES.Homo_sapiens:
            return "Human"
        elif self == SPECIES.Pan_troglodytes:
            return "Chimpanzee"
        elif self == SPECIES.Pongo_abelii:
            return "Sumatran orangutan"
        elif self == SPECIES.Papio_anubis:
            return "Olive baboon"
        elif self == SPECIES.Nomascus_leucogenys:
            return "Gibbon"
        elif self == SPECIES.Callithrix_jacchus:
            return "Common marmoset"
        elif self == SPECIES.Rattus_norvegicus:
            return "Rat"
        elif self == SPECIES.Mus_musculus:
            return "Mouse"
        elif self == SPECIES.Sus_scrofa:
            return "Pig"
        elif self == SPECIES.Canis_lupus_familiaris:
            return "Dog"
        elif self == SPECIES.Gallus_gallus:
            return "Chicken"
        elif self == SPECIES.Danio_rerio:
            return "Zebrafish"
        elif self == SPECIES.Drosophila_melanogaster:
            return "D. melanogaster"
        elif self == SPECIES.E_coli:
            return "E. coli"
        elif self == SPECIES.Arabidopsis_thaliana:
            return "Flower (thale cress)"
        elif self == SPECIES.Saccharomyces_cerevisiae:
            return "Baker's yeast"
        else:
            assert False

    def __str__(self):
        return str(self.name).replace('_', ' ')


class TRANSCRIPT_CRITERIA(enum.Enum):
    NONE = 0
    LONGEST = 1
    LONGEST_CDS = 2
    LONGEST_CDS_AND_UTRs = 3  # similar to exons
    RANDOM = 4
