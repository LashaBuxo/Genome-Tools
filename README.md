 

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
Number of genes: 19957
Number of genes on Positive(+) Strand: 10092
Number of genes on Negative(-) Strand: 9865

Number of OG: 6795
Number of OG pairs: 4999
Number of OG clusters with >=2 gene: 2513

Genes by clusters length:
1-length clusters: 13162
2-length clusters: 1539
3-length clusters: 592
4-length clusters: 215
5-length clusters: 89
6-length clusters: 33
7-length clusters: 18
8-length clusters: 11
9-length clusters: 2
10-length clusters: 3
11-length clusters: 3
12-length clusters: 1
13-length clusters: 1
15-length clusters: 2
16-length clusters: 1
17-length clusters: 1
22-length clusters: 1
33-length clusters: 1
```
</td>
<td>

```yaml
Number of genes: 19846
Number of genes on Positive(+) Strand: 10029
Number of genes on Negative(-) Strand: 9817

Number of OG: 4533
Number of OG pairs: 3000
Number of OG clusters with >=2 gene: 1933

Genes by clusters length:
1-length clusters: 15313
2-length clusters: 1516
3-length clusters: 314
4-length clusters: 62
5-length clusters: 17
6-length clusters: 9
7-length clusters: 6
8-length clusters: 1
9-length clusters: 1
11-length clusters: 1
13-length clusters: 1
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

![](results/images/Gene%20Outline%20(Ensembl%2C%20k%3D50%2C%20procc%3D1000).png)
 
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
Folder [/results](/results) contains all the exported .txt files, where at the end of each file
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
genome-build Genome Reference Consortium GRCh38.p13

Genome Primary Assembly:
GenBank Assembly ID
genome-build-accession GCA_000001405.28  

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