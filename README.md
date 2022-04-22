# Overlapping genes in Human Genome (general stats)
In project there are used specific (latest for 2022/04) annotations from 
`NCBI (Refseq)` and `Ensembl` with their corresponding genome builds: **Refseq** and **Genbank** accordingly.
##### _Specifics:_
 

<table>
<tr>
<th> Ensembl version  </th>
<th> NCBI (Refseq) versuib </th>
</tr>
<tr>
<td> 

```json
Genome Reference Assembly:
	genome-build Genome Reference Consortium GRCh38.p13

Genome Primary Assembly:
	GenBank Assembly ID 
	genome-build-accession GCA_000001405.28  

Genome Annotation:
	Ensembl Release 105 (2021)
```
</td>
<td>
 

```json 
Genome Reference Assembly: 
	genome-build: Genome Reference Consortium GRCh38.p13

Genome Primary Assembly:
	NCBI Assembly ID 
	genome-build-accession: GCF_000001405.39

Genome Annotation:
	NCBI Homo sapiens Updated Annotation Release 109.20211119 
	11/19/2021
```
 

</td>
</tr>
</table>

# Overlapping genes in Human Genome (general stats)
`overlapping_genes_stats.py` does general stats calculation on total genes, overlapping genes, overlapping gene clusters, etc... 
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
Genes are chosen from genome file features,if feature has `gene` **feature type** and if it has `protein_coding` attribute for **biotype**.

<table>
<tr>
<th> Using Ensembl </th>
<th> Using NCBI (Refseq) </th>
</tr>
<tr>
<td> 

```yaml
Number of genes: 19951
Number of genes on Positive(+) Strand: 10086
Number of genes on Negative(-) Strand: 9865

Number of OG: 6758
Number of OG pairs: 4972
Number of OG clusters with >=2 gene: 2502

Genes by clusters length:
1-length clusters: 13193
2-length clusters: 1536
3-length clusters: 590
4-length clusters: 210
5-length clusters: 88
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
Number of genes: 19476
Number of genes on Positive(+) Strand: 9852
Number of genes on Negative(-) Strand: 9624

Number of OG: 4420
Number of OG pairs: 2917
Number of OG clusters with >=2 gene: 1893

Genes by clusters length:
1-length clusters: 15056
2-length clusters: 1489
3-length clusters: 307
4-length clusters: 61
5-length clusters: 15
6-length clusters: 6
7-length clusters: 6
8-length clusters: 2
9-length clusters: 1
11-length clusters: 1
13-length clusters: 1
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

![](results/images/average%20gene%20(Ensembl%2C%20k%3D50%2C%20procc%3D1000).png)
 
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
#### [overlapped genes (+descriptions) - CDS [diff_stranded, both_framed] (Ensembl).txt](/results/overlapped%20genes%20(+descriptions)%20-%20CDS%20[diff_stranded,%20both_framed]%20(Ensembl).txt)
#### [overlapped genes (+descriptions) - CDS [diff_stranded, both_framed] (NCBI)](/results/overlapped%20genes%20(+descriptions)%20-%20CDS%20[diff_stranded,%20both_framed]%20(NCBI).txt) 