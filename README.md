 

## Overlapping genes in Human Genome (general stats)
`overlapping_genes_stats.py` does general stats calculation on total genes, overlapping genes, overlapping gene clusters, etc... 

Genes are chosen from genome file features, if feature has `gene` **feature type** and if it has `protein_coding` attribute for **biotype**.

##### _Specifics:_
```text
Usage:
  overlapping_genes_stats.py <annotation>

Params (possible) to run:
  annotation: NCBI / Ensembl

Example:
  python overlapping_genes_stats.py NCBI

Output:
  prints stats in console
```
##### _Example Output:_
 
<table>
<tr>
<th> Using Ensembl </th>
<th> Using NCBI (Refseq) </th>
</tr>
<tr>
<td> 

```yaml
--------------------------clusters and their parameters-----------------------------------

cluster #1  |  size: 2 
	DNAJC16 = strand + | isoforms: 6 | expr. overall: 12.5[nTPM] | desc: DnaJ heat shock protein family (Hsp40) member C16 [Source:HGNC Symbol;Acc:HGNC:29157]
	AGMAT = strand - | isoforms: 1 | expr. overall: 8.2[nTPM] | desc: agmatinase [Source:HGNC Symbol;Acc:HGNC:18407]
	transcript:ENST00000375849 (DNAJC16) - transcript:ENST00000375826 (AGMAT):
		specifics: type - Convergent | len - 19 | exons - [5'-uuxxxxxxxxxxxxxXu-3'], [5'-uxXxxxxxu-3']
		flanks(+): [TTTAAGCTACCGGCTGCCCATGGATCTCAGTCTCTTTGCTTACTCACCAG]-[overlap]-[AACTTCGACAAGATCACAGCCCATCACGTTCAGGCCTTGACAACCCCTGA] 
		flanks(-): [CTGGTGAGTAAGCAAAGAGACTGAGATCCATGGGCAGCCGGTAGCTTAAA]-[overlap]-[TCAGGGGTTGTCAAGGCCTGAACGTGATGGGCTGTGATCTTGTCGAAGTT] 
		seq 1: 5'...A|AAGATCATACGGTGGTGA3'...
				GC - 42.11%
		seq 2: 5'...TCACCACCGTATGATCTTT3'...
				GC - 42.11%
		peptide 1: N...(646-652)KIIRW*...C
				TOP 3 aa: I (33.33%), R (16.67%), K (16.67%)
				degeneracy: medium (40.00%), low (40.00%), high (20.00%)
				aa types: non-polar (60.00%), + (40.00%), - (0.00%), polar (0.00%)
		peptide 2: N...(320-326)SPPYDL...C
				TOP 3 aa: P (33.33%), L (16.67%), S (16.67%)
				degeneracy: high (33.33%), medium (33.33%), low (33.33%)
				aa types: non-polar (50.00%), polar (33.33%), - (16.67%), + (0.00%)
		peptide hits 1:
				[6-652] -> 7 | PANTHER |  | PTHR44303 | 
				[6-652] -> 7 | PANTHER |  | PTHR44303:SF3 | 
		peptide hits 2:
				[72-344] -> 7 | TIGRfam | Agmatinase-related | TIGR01230 | IPR005925
				[1-352] -> 7 | alphafold |  | AF-Q9BSE5-F1.A | 
				[41-352] -> 7 | Gene3D |  | 3.40.800.10 | 
				[72-343] -> 7 | Pfam | Ureohydrolase | PF00491 | IPR006035
				[46-350] -> 7 | PIRSF | Ureohydrolase | PIRSF036979 | IPR006035
				[48-350] -> 7 | PANTHER | Ureohydrolase | PTHR11358 | IPR006035
				[61-352] -> 7 | Prosite_profiles | Ureohydrolase | PS51409 | IPR006035
				[55-343] -> 7 | CDD |  | cd11592 | 
				[53-346] -> 7 | SuperFamily | Ureohydrolase domain superfamily | SSF52768 | IPR023696

cluster #2  |  size: 2 
	HCRTR1 = strand + | isoforms: 4 | expr. overall: 1.7[nTPM] | desc: hypocretin receptor 1 [Source:HGNC Symbol;Acc:HGNC:4848]
	PEF1 = strand - | isoforms: 1 | expr. overall: 10.3[nTPM] | desc: penta-EF-hand domain containing 1 [Source:HGNC Symbol;Acc:HGNC:30009]
	transcript:ENST00000373705 (HCRTR1) - transcript:ENST00000373703 (PEF1):
		specifics: type - Convergent | len - 24 | exons - [5'-xxxxxxX-3'], [5'-uxXxxxu-3']
		flanks(+): [AGAAGAGTCTAGCTCTGTCCTGCCCATCGTGCCCCGGCCATGACCCGCAC]-[overlap]-[TGGAGCCCGAGCGGTCCCGGTCATACTGCTGGAAGAGGTTCTTCCACTGC] 
		flanks(-): [GTGCGGGTCATGGCCGGGGCACGATGGGCAGGACAGAGCTAGACTCTTCT]-[overlap]-[GCAGTGGAAGAACCTCTTCCAGCAGTATGACCGGGACCGCTCGGGCTCCA] 
		seq 1: 5'...CTTGCTGCAGCTCTGTGTAGCTAA3'...
				GC - 50.00%
		seq 2: 5'...TT|AGCTACACAGAGCTGCAGCAAG3'...
				GC - 50.00%
		peptide 1: N...(381-388)LAAALCS*...C
				TOP 3 aa: A (37.50%), L (25.00%), S (12.50%)
				degeneracy: high (42.86%), medium (42.86%), low (14.29%)
				aa types: non-polar (85.71%), polar (14.29%), + (0.00%), - (0.00%)
		peptide 2: N...(200-208)SYTELQQ...C
				TOP 3 aa: Q (28.57%), L (14.29%), S (14.29%)
				degeneracy: low (57.14%), high (28.57%), medium (14.29%)
				aa types: polar (71.43%), - (14.29%), non-polar (14.29%), + (0.00%)
		peptide hits 1:
		peptide hits 2:
				[118-283] -> 9 | CDD |  | cd16184 | 
				[181-216] -> 9 | Prosite_profiles | EF-hand domain | PS50222 | IPR002048
				[1-284] -> 9 | alphafold |  | AF-Q9UBV8-F1.A | 
				[114-279] -> 9 | SuperFamily | EF-hand domain pair | SSF47473 | IPR011992
				[186-212] -> 9 | Pfam | EF-hand domain | PF13405 | IPR002048
				[1-284] -> 9 | PANTHER |  | PTHR23064 | 
				[1-284] -> 9 | PANTHER |  | PTHR23064:SF45 | 
				[194-206] -> 7 | Prosite_patterns | EF-Hand 1, calcium-binding site | PS00018 | IPR018247
				[185-213] -> 9 | Smart | EF-hand domain | SM00054 | IPR002048
				[108-283] -> 9 | Gene3D |  | 1.10.238.10 | 

cluster #3  |  size: 2 
	ZBTB8A = strand + | isoforms: 2 | expr. overall: 11.0[nTPM] | desc: zinc finger and BTB domain containing 8A [Source:HGNC Symbol;Acc:HGNC:24172]
	ZBTB8OS = strand - | isoforms: 11 | expr. overall: 20.5[nTPM] | desc: zinc finger and BTB domain containing 8 opposite strand [Source:HGNC Symbol;Acc:HGNC:24094]
	transcript:ENST00000373510 (ZBTB8A) - transcript:ENST00000341885 (ZBTB8OS):
		specifics: type - Convergent | len - 98 | exons - [5'-uuuxxXu-3'], [5'-uXxu-3']
		flanks(+): [AGAATTTGGCATAGACAGCCTCCCCATTGACTTGGAAGCTGAACAACATC]-[overlap]-[CTGGTCATCCAACAGGTTGATGATAGTGAAGAAGAAGAAGAAAAAGAAAT] 
		flanks(-): [GATGTTGTTCAGCTTCCAAGTCAATGGGGAGGCTGTCTATGCCAAATTCT]-[overlap]-[ATTTCTTTTTCTTCTTCTTCTTCACTATCATCAACCTGTTGGATGACCAG] 
		seq 1: 5'...TT|ATGTCCCCATCAGATGGAGATAAGGATTCCAGATGGCACTTGAGTGAAGATGAGAATAGATCCTATGTGGAGATTGTAGAAGATGGGTCTGCTGAT3'...
				GC - 42.86%
		seq 2: 5'...AT|CAGCAGACCCATCTTCTACAATCTCCACATAGGATCTATTCTCATCTTCACTCAAGTGCCATCTGGAATCCTTATCTCCATCTGATGGGGACATAA3'...
				GC - 42.86%
		peptide 1: N...(386-419)MSPSDGDKDSRWHLSEDENRSYVEIVEDGSAD...C
				TOP 3 aa: S (18.75%), D (18.75%), E (12.50%)
				degeneracy: low (50.00%), high (28.12%), medium (21.88%)
				aa types: - (31.25%), non-polar (31.25%), polar (25.00%), + (12.50%)
		peptide 2: N...(33-65)QQTHLLQSPHRIYSHLHSSAIWNPYLHLMGT*...C
				TOP 3 aa: L (15.62%), H (15.62%), S (12.50%)
				degeneracy: low (41.94%), high (32.26%), medium (25.81%)
				aa types: non-polar (41.94%), polar (38.71%), + (19.35%), - (0.00%)
		peptide hits 1:
				[273-439] -> 34 | PANTHER |  | PTHR24414 | 
				[273-439] -> 34 | PANTHER |  | PTHR24414:SF25 | 
				[1-441] -> 34 | alphafold |  | AF-Q96BR9-F1.A | 
		peptide hits 2:
				[1-35] -> 3 | PANTHER | Archease | PTHR12682 | IPR002804
				[1-35] -> 3 | PANTHER |  | PTHR12682:SF12 | 

cluster #4  |  size: 2 
	PPIE = strand + | isoforms: 10 | expr. overall: 29.7[nTPM] | desc: peptidylprolyl isomerase E [Source:HGNC Symbol;Acc:HGNC:9258]
	BMP8B = strand - | isoforms: 1 | expr. overall: 5.9[nTPM] | desc: bone morphogenetic protein 8b [Source:HGNC Symbol;Acc:HGNC:1075]
	transcript:ENST00000356511 (PPIE) - transcript:ENST00000372827 (BMP8B):
		specifics: type - Convergent | len - 31 | exons - [5'-uxxxxxxxxxXu-3'], [5'-uxxXxxxxu-3']
		flanks(+): [ATCTACAGGTGGTGATTGTCATTTCAGAAACAAGAAGAGTCAGCAATTAC]-[overlap]-[AGCTCGTGCCGACGGCAGACCTGCCGGCCGTGGGAGCCGTGGACGTCATC] 
		flanks(-): [GTAATTGCTGACTCTTCTTGTTTCTGAAATGACAATCACCACCTGTAGAT]-[overlap]-[GATGACGTCCACGGCTCCCACGGCCGGCAGGTCTGCCGTCGGCACGAGCT] 
		seq 1: 5'...C|AGCCAGCCGAGGTCCTGGAAGCTGACGTAG3'...
				GC - 64.52%
		seq 2: 5'...C|TACGTCAGCTTCCAGGACCTCGGCTGGCTG3'...
				GC - 64.52%
		peptide 1: N...(284-294)SQPRSWKLT*...C
				TOP 3 aa: S (20.00%), L (10.00%), R (10.00%)
				degeneracy: high (44.44%), low (33.33%), medium (22.22%)
				aa types: polar (44.44%), non-polar (33.33%), + (22.22%), - (0.00%)
		peptide 2: N...(305-315)YVSFQDLGWL...C
				TOP 3 aa: L (20.00%), S (10.00%), G (10.00%)
				degeneracy: low (50.00%), high (30.00%), medium (20.00%)
				aa types: non-polar (60.00%), polar (30.00%), - (10.00%), + (0.00%)
		peptide hits 1:
				[131-292] -> 9 | Gene3D | Cyclophilin-like domain superfamily | 2.40.100.10 | IPR029000
				[143-288] -> 5 | Pfam | Cyclophilin-type peptidyl-prolyl cis-trans isomerase domain | PF00160 | IPR002130
				[1-291] -> 8 | PIRSF | Peptidyl-prolyl cis-trans isomerase E | PIRSF001475 | IPR016304
				[143-290] -> 7 | Prosite_profiles | Cyclophilin-type peptidyl-prolyl cis-trans isomerase domain | PS50072 | IPR002130
		peptide hits 2:
				[291-402] -> 11 | Gene3D | Cystine-knot cytokine | 2.10.90.10 | IPR029034
				[301-318] -> 11 | PRINTS |  | PR00669 | 
				[298-402] -> 11 | CDD |  | cd19398 | 
				[301-401] -> 11 | Pfam | Transforming growth factor-beta, C-terminal | PF00019 | IPR001839
				[1-402] -> 11 | alphafold |  | AF-P34820-F1.A | 
				[268-402] -> 11 | Prosite_profiles | Transforming growth factor-beta, C-terminal | PS51362 | IPR001839
				[293-401] -> 11 | SuperFamily | Cystine-knot cytokine | SSF57501 | IPR029034
				[301-402] -> 11 | Smart | Transforming growth factor-beta, C-terminal | SM00204 | IPR001839
				[4-402] -> 11 | PANTHER | Transforming growth factor-beta-related | PTHR11848 | IPR015615
				[4-402] -> 11 | PANTHER |  | PTHR11848:SF119 | 


```
</td>
<td>

```yaml
Organism: Homo sapiens
Annotation: NCBI
protein coding genes: 19235 (611 filtered out)

Genes on Positive(+) Strand: 9743
Genes on Negative(-) Strand: 9492

Overlapping Genes: 3860
Overlapping Gene clusters (>1 gene): 1692
     2-length clusters: 1408
     3-length clusters: 222
     4-length clusters: 33
     5-length clusters: 8
     6-length clusters: 7
     7-length clusters: 5
     8-length clusters: 1
     9-length clusters: 1
     11-length clusters: 1
     12-length clusters: 1
     14-length clusters: 1
     15-length clusters: 1
     17-length clusters: 1
     21-length clusters: 1
     22-length clusters: 1
```

</td>
</tr>
</table>
 
## Gene model in human genome by base occurrences 
 `base_occurrences.py` Builds 'Average' gene outline in genome (as .png).
  Genes is divided into 6 region: `UTR'5_Procession`, `UTR'5`, `inner CDSs`, `inner Introns`, `UTR'3`,` UTR'3_Procession`
  and each region is divided into _**k subregions**_, where _ATCG_ occurrences are calculated.
  
Based each nucleotide occurrences, graph is built from _**6 x k point**_
  where X axis corresponds to points, and Y axis corresponds  to nucleotide occurences.
##### _Specifics:_
```text 
Usage:
  base_occurrences.py <annotation> <subregions> <procession_length>

Params (possible) to run:
  annotation: NCBI / Ensembl
  subregions: k - number of subregions within each region
  procession_length: l - length of processions at the end of gene (before UTR5' or after UTR3')

Example:
  python base_occurrences.py Ensembl 50 1000

Output:
  created image in ./results/images folder

```
##### _Example Output:_

![](other/results/images/Gene%20Outline%20(Ensembl%2C%20k%3D50%2C%20procc%3D1000).png)
 
## Overlapping genes by specific regions 
 `overlapping_fragments.py`  Calculates overlapping gene pairs by
  only specific parts: CDS (coding sequence), EXON or whole gene
##### _Specifics:_
```text 
Usage:
  overlapping_fragments.py <annotation> <fragment_type> <strand_similarity> <ORF_similarity> <output_type>

Params (possible) to run:
  annotation: NCBI / Ensembl
  fragment_type: CDS / exon / gene
  strand_similarity: same_stranded / diff_stranded / both_stranded
  ORF_similarity: same_ORF / diff_ORF / both_ORF
  output_type: print_gene_pairs / print_gene_pairs_with_descriptions

Example:
  python overlapping_fragments.py NCBI CDS diff_stranded both_framed print_gene_pairs_with_descriptions
  Finds all overlapping genes by CDS fragments, which are on different strands
  and which have same ORF or different ORF

Output:
  creates .txt file in ./results folder
```
##### _Example Outputs:_
Folder [/results](/other/results) contains all the exported .txt files, where at the end of each file
there is statistic about composition of overlapped regions.


<hr style="border:2px solid gray"> </hr>

## Human Genome annotation/assembly versions
In project there are used specific (latest for 2022/04) annotations from 
`NCBI (Refseq)` and `Ensembl` with their corresponding genome builds: **Refseq** and **Genbank** accordingly.

Annotation and sequence files must be located in the `/genome_data/NCBI` or/and `/genome_data/Ensembl` folder and then 
values in `genome_lib_values.py` must be changed accordingly. 
##### _Specifics:_
 

<table>
<tr>
<th> Ensembl version  </th>
<th> NCBI (Refseq) version </th>
</tr>
<tr>
<td> 

```yaml
Genome Reference Assembly:
Genome Worker-build Genome Reference Consortium GRCh38.p13

Genome Primary Assembly:
GenBank Assembly ID
Genome Worker-build-accession GCA_000001405.28  

Genome Annotation:
Ensembl Release 106 (Apr 2022)
```
</td>
<td>
 

```yaml 
Genome Reference Assembly:
genome-build Genome Reference Consortium GRCh38.p14

Genome Primary Assembly:
NCBI Assembly ID 
genome-build-accession GCF_000001405.40

Genome Annotation:
NCBI Homo sapiens Annotation Release 110 (06/04/2022)

```
 

</td>
</tr>
</table>
