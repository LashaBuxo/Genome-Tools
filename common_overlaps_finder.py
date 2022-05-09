file_human = "./generated_data/observed/OGs by CDS [Human, NCBI, diff_stranded, both_ORF].txt"
file_mouse = "./generated_data/observed/OGs by CDS [Mouse, NCBI, diff_stranded, both_ORF].txt"


def readgenes(path):
    file = open(path, 'r')
    lines = file.readlines()
    genes = []
    last_gene = ''
    for line in lines:
        if line.__contains__("Gene 1:") or line.__contains__("Gene 2:"):
            info = line.split('\t')[1].split(' ')[2]
            info = info.lower()
            if last_gene != '':
                genes.append((info, last_gene))
                last_gene = ''
            else:
                last_gene = info

    return genes


human_genes = readgenes(file_human)
mouse_genes = readgenes(file_mouse)
for gene in human_genes:
    flag = False
    for gene_m in mouse_genes:
        if (gene[0] == gene_m[0] or gene[0] == gene_m[1]) and (gene[1] == gene_m[0] or gene[1] == gene_m[1]):
            flag = True
    if flag == True: print(gene)
