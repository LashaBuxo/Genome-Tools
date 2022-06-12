from generator_writer import *

genome = GenomeWorker(SPECIES.Homo_sapiens, ANNOTATIONS.ENSEMBL, ANNOTATION_LOAD.GENES, SEQUENCE_LOAD.NOT_LOAD)
genome2 = GenomeWorker(SPECIES.Mus_musculus, ANNOTATIONS.ENSEMBL, ANNOTATION_LOAD.GENES, SEQUENCE_LOAD.NOT_LOAD)


def calc(gen: GenomeWorker):
    mito_genes = []
    mito_overlapped_genes = []
    mito_overlapped_to_mito_genes = []

    for chr_id in range(1, gen.chromosomes_count() + 1):
        genes_cnt = gen.genes_count_on_chr(chr_id)
        for ind in range(0, genes_cnt):
            gene = gen.gene_by_indexes(chr_id, ind)
            if gen.is_gene_MITO(gene.id):
                mito_genes.append(gene)
    for mito_gene in mito_genes:
        flag_overlapped = False
        flag_overlapped_to_mito = False
        for chr_id in range(1, gen.chromosomes_count() + 1):
            genes_cnt = gen.genes_count_on_chr(chr_id)
            for ind in range(0, genes_cnt):
                gene = gen.gene_by_indexes(chr_id, ind)
                if gene.chrom != mito_gene.chrom: continue
                if gene.id == mito_gene.id: continue
                if gen.are_segments_overlapped((mito_gene.start, mito_gene.end), (gene.start, gene.end)):
                    flag_overlapped = True
                    if gen.is_gene_MITO(gene.id):
                        if gen.get_gene_symbol(mito_gene) == 'ACOT7': print(f"lasha:{gen.get_gene_symbol(gene)}")
                        flag_overlapped_to_mito = True
        if flag_overlapped:
            mito_overlapped_genes.append(mito_gene)
        if flag_overlapped_to_mito: mito_overlapped_to_mito_genes.append(mito_gene)
    return mito_genes, mito_overlapped_genes, mito_overlapped_to_mito_genes


def find_common_set(arr1, arr2):
    ids = {}
    for gene in arr1:
        ids[genome.get_gene_symbol(gene)] = True
    common = []
    for gene in arr2:
        if ids.__contains__(genome.get_gene_symbol(gene)):
            common.append(genome.get_gene_symbol(gene))
    return common


mito_human, mito_human_overlapped, mito_human_overlapped_to_mito = calc(genome)
mito_mouse, mito_mouse_overlapped, mito_mouse_overlapped_to_mito = calc(genome2)

mito_common_genes = find_common_set(mito_human, mito_mouse)


mito_common_overlapped = find_common_set(mito_human_overlapped, mito_mouse_overlapped)
mito_common_overlapped_and_paired = []
for mito_sym in mito_common_overlapped:
    flag = False
    for chr_id in range(1, genome.chromosomes_count() + 1):
        genes_cnt = genome.genes_count_on_chr(chr_id)
        for ind in range(0, genes_cnt):
            gene = genome.gene_by_indexes(chr_id, ind)
            gene_sym = genome.get_gene_symbol(gene)
            if gene_sym == mito_sym: continue
            if genome.are_genes_overlapped(gene_sym, mito_sym) and genome2.are_genes_overlapped(gene_sym, mito_sym):
                flag = True
    if flag: mito_common_overlapped_and_paired.append(mito_sym)

mito_common_overlapped_to_mito = find_common_set(mito_human_overlapped_to_mito, mito_mouse_overlapped_to_mito)
mito_common_overlapped_to_mito_and_paired = []
for mito_sym in mito_common_overlapped_to_mito:
    flag = False
    for chr_id in range(1, genome.chromosomes_count() + 1):
        genes_cnt = genome.genes_count_on_chr(chr_id)
        for ind in range(0, genes_cnt):
            gene = genome.gene_by_indexes(chr_id, ind)
            gene_sym = genome.get_gene_symbol(gene)
            if gene_sym == mito_sym: continue
            if not genome.is_gene_MITO(gene.id): continue
            if genome.are_genes_overlapped(gene_sym, mito_sym) and genome2.are_genes_overlapped(gene_sym, mito_sym):
                flag = True
    if flag: mito_common_overlapped_to_mito_and_paired.append(mito_sym)


output_data_path = f'generated_data/texts/OGs mito stats.txt'

with open(output_data_path, 'w') as f:
    f.write(f'mito genes human - {len(mito_human)}; mouse - {len(mito_mouse)}\n')
    f.write(f'mito common genes - {len(mito_common_genes)}\n')

    f.write(f'mito overlapped genes: human - {len(mito_human_overlapped)}; mouse - {len(mito_mouse_overlapped)}\n')
    f.write(f'mito common overlapped genes - {len(mito_common_overlapped)}\n')
    f.write(f'mito common overlapped genes and paired - {len(mito_common_overlapped_and_paired)}\n\n')

    f.write(f'mito overlapped to mito genes: human - {len(mito_human_overlapped_to_mito)}; mouse - {len(mito_mouse_overlapped_to_mito)}\n')
    f.write(f'mito common overlapped to mito genes - {len(mito_common_overlapped_to_mito)}\n')
    f.write(f'mito common overlapped to mito genes and paired - {len(mito_common_overlapped_to_mito_and_paired)}\n')