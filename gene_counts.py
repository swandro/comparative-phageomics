

cutoffs = [60,65,70,75,80,85,90,92,95,97,99]

nt_file = open("gene_counts/nt_counts.tsv",'w')
aa_file = open("gene_counts/aa_counts.tsv",'w')

nt_file.write("number\t" + '\t'.join([str(x) for x in range(1,17)]) + "\ttotal\t" + '\n')
aa_file.write("number\t" + '\t'.join([str(x) for x in range(1,17)]) + "\ttotal\t" + '\n')

for cutoff in cutoffs:
    nt_file.write(str(cutoff) + '\t')
    aa_file.write(str(cutoff) + '\t')
    #Different cutoffs
    nt_blast_parsed_cutoff = parse_blast("blast_output/nt_blast.txt", cutoff, 50)
    aa_blast_parsed_cutoff = parse_blast("blast_output/aa_blast.txt", cutoff, 50)

    #Make master gene lists
    nt_master_genes_cutoff = make_master_gene_list(all_genes, nt_blast_parsed_cutoff)
    aa_master_genes_cutoff = make_master_gene_list(all_genes, aa_blast_parsed_cutoff)

    freq_nt = gene_frequency(nt_master_genes_cutoff)
    freq_aa = gene_frequency(aa_master_genes_cutoff)

    nt_result = []
    aa_result = []
    for i in range(1,17):
        nt_result.append(freq_nt[i])
        aa_result.append(freq_aa[i])

    nt_file.write('\t'.join([str(x) for x in nt_result]) + '\t' + str(sum(freq_nt.values())) + '\n')
    aa_file.write('\t'.join([str(x) for x in aa_result]) + '\t' + str(sum(freq_aa.values())) + '\n')

nt_file.close()
aa_file.close()


##Finding most conserved gene
aa_blast_parsed_cutoff = parse_blast("blast_output/aa_blast.txt", 95, 50)
aa_master_genes_cutoff = make_master_gene_list(all_genes, aa_blast_parsed_cutoff)

for gene in aa_master_genes_cutoff:
    if len(aa_master_genes_cutoff[gene].keys())==16:
        print(gene)

aa_master_genes_cutoff[21]
for gene in all_genes_dict["V12_RAST"]:
    if gene.id_only == "FLJIFKPE_00143":
        a=gene

str(a.seq)