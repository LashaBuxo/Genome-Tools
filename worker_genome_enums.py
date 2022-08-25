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
    DIFF_NESTED_SHORTER = 3  # different strand + nested
    DIFF_NESTED_LONGER = 4  # different strand + nested
    DIVERGENT = 5
    SAME_NESTED = 6  # one of the gene entirely located into others boundaries on same strand
    TANDEM = 7  # genes located on same strand and neither is nested

    MULTI = 8  # for set of genes

    def short_name(self):
        if self == OVERLAP_TYPE.CONVERGENT:
            return "Convergent"
        elif self == OVERLAP_TYPE.DIFF_NESTED:
            return "Nested (diff strand)"
        elif self == OVERLAP_TYPE.DIFF_NESTED_SHORTER:
            return "Nested (diff strand, shorter)"
        elif self == OVERLAP_TYPE.DIFF_NESTED_LONGER:
            return "Nested (diff strand, longer)"
        elif self == OVERLAP_TYPE.DIVERGENT:
            return "Divergent"
        elif self == OVERLAP_TYPE.SAME_NESTED:
            return "Nested (same strand)"
        elif self == OVERLAP_TYPE.TANDEM:
            return "Tandem"
        elif self == OVERLAP_TYPE.PROMOTOR:
            return "PO"
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
    Rattus_norvegicus = 2
    Mus_musculus = 3
    Sus_scrofa = 4
    Xenopus_tropicalis = 5
    Danio_rerio = 6
    Tetraodon_nigroviridis = 7
    Caenorhabditis_elegans = 8
    Drosophila_melanogaster = 9
    Arabidopsis_thaliana = 10

    def shortest_name(self):
        return self.short_name()

    def taxon_ID(self):
        if self == SPECIES.Homo_sapiens:
            return "9606"
        elif self == SPECIES.Pan_troglodytes:
            return "9598"
        elif self == SPECIES.Rattus_norvegicus:
            return "10116"
        elif self == SPECIES.Mus_musculus:
            return "10090"
        elif self == SPECIES.Sus_scrofa:
            return "9823"
        elif self == SPECIES.Xenopus_tropicalis:
            return "8364"
        elif self == SPECIES.Danio_rerio:
            return "7955"
        elif self == SPECIES.Tetraodon_nigroviridis:
            return "99883"
        elif self == SPECIES.Caenorhabditis_elegans:
            return "6239"
        elif self == SPECIES.Drosophila_melanogaster:
            return "7227"
        elif self == SPECIES.Arabidopsis_thaliana:
            return "3702"
        else:
            assert False

    def short_name(self):
        if self == SPECIES.Homo_sapiens:
            return "Human"
        elif self == SPECIES.Pan_troglodytes:
            return "Chimpanzee"
        elif self == SPECIES.Rattus_norvegicus:
            return "Rat"
        elif self == SPECIES.Mus_musculus:
            return "Mouse"
        elif self == SPECIES.Sus_scrofa:
            return "Pig"
        elif self == SPECIES.Xenopus_tropicalis:
            return "Frog"
        elif self == SPECIES.Danio_rerio:
            return "ZebraFish"
        elif self == SPECIES.Tetraodon_nigroviridis:
            return "Tetraodon"
        elif self == SPECIES.Caenorhabditis_elegans:
            return "C. elegans"
        elif self == SPECIES.Drosophila_melanogaster:
            return "D. melanogaster"
        elif self == SPECIES.Arabidopsis_thaliana:
            return "Flower"
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
