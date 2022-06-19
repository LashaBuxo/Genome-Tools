## Genome Working Tools

Python project is a multi-functional BioInformatics tool, which simplifies calculations on genome annotations provided
by Ensembl or RefSeq (NCBI).

## Instruction to use

The project requires desired annotations
(GFF format), sequence files (FASTA format)
or other required files for calculations to be imported into the project and declared
in [worker_genome_values.py](worker_genome_values.py) file.

Project supports annotations and sequence files format provided
from [RefSeq (NCBI)](https://ftp.ncbi.nlm.nih.gov/genomes/) or [Ensembl](https://asia.ensembl.org/index.html). In
addition to adding appropriate files, specific organism name must be declared
in [worker_genome_values.py](worker_genome_values.py) as well, for which files are added.

_Note:_ NCBI .gff files lacks UTR 5' and UTR 3' features. In order to add these features run script add_ustrs_to_gff.py
after that use any converter tool to change the text format to the ANSII Unicode.
[add_ustrs_to_gff.py](used_data/genome_data/NCBI/add_utrs_to_gff.py)
is written by "David Managadze" (works at NCBI) and this script was suggested for use in
NCBI's [readme files](https://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt).

### Genome Load

```Python
from worker_genome import *

genome = GenomeWorker(SPECIES.Homo_sapiens,  # For which organism genome must be loaded?
                      ANNOTATIONS.ENSEMBL,  # Which annotation to use?
                      ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS,  # What type of features are required from annotation?
                      SEQUENCE_LOAD.NOT_LOAD)  # Should sequence also be loaded?
```

### Sequence Analyzer Usage

```Python
from worker_analyzer import *

analyzer = AnalyzerData()

analyzer.analyze_sequence_stats("CCAGCAGCAG",  # Nucleotide Sequence
                                1)  # ORF for sequence
print(analyzer.analyzed_peptide)  # Output: 'QQQ' 
analyzer.analyze_sequence_stats("ATGATG", 0)
print(analyzer.analyzed_peptide == 'QQQMM')  # Output: 'QQQMM'
print(analyzer.get_gc_content())  # Output: 0.64285714285
```

## Example Usage

Calculate protein-coding genes greater than 10000 (base pair) length:

```Python
total_genes = 0
genes_greater_1000000 = 0

for chr_id in range(1, genome.chromosomes_count() + 1):
    genes_cnt = genome.genes_count_on_chr(chr_id)
    total_genes += genes_cnt
    for i in range(0, genes_cnt):
        gene = genome.gene_by_indexes(chr_id, i)
        genes_greater_1000000 += 1 if gene.end - gene.start + 1 > 10000 else 0
print(f'{genes_greater_1000000} genes from {total_genes}')
```

```Console
14227 genes from 19176
```

Based on all protein-coding genes, it can be used to outline gene GC content, by calculating GC content in specific
regions of the gene and taking average value from all the gene. Regions are divided here into k=50 sub-regional parts.

Used organisms in analyses are: Caenorhabditis elegans, Drosophila melanogaster, Homo sapiens, Danio Rerio, Mus
Musculus, Rat norvegicus, Homo sapiens

#### When Ensembl annotated protein-coding genes are used:

![alt text for screen readers](Comparative%20Gene%20Outline%20(Ensembl,%20k=50,%20procc=1000).png "Text to show on mouseover")

#### When RefSeq (NCBI) annotated protein-coding genes are used:

![alt text for screen readers](Comparative%20Gene%20Outline%20(NCBI,%20k=50,%20procc=1000).png "Text to show on mouseover")

## 3rd Party Resources

* [gffutils](https://github.com/daler/gffutils) - We use gffutils for working with large GFF files
