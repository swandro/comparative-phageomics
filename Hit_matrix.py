from Bio import SeqIO
from Bio import AlignIO
from Bio import Seq
import glob
import subprocess



##############Save the parsed BLAST output as a hit matrix for each genome###########

#Helper function to write a hit table. ########################
def write_hit_matrix(genome, blast_output, folder):
    '''Function to save a hit matrix. Takes in a genome name, blast output, and destination folder
    In the hit matrix: Rows=genes, columns = other phages, fill = % identity
    '''
    with open( (folder + "/" + genome + "_hit_matrix.txt"), 'w') as outfile:
        #Get list of other genomes
        other_genome_names = [x for x in all_genes_dict.keys()]
        other_genome_names.remove(genome)

        #Write header
        firstline = ["gene"]
        firstline.extend(other_genome_names)
        outfile.write('\t'.join(firstline) + '\n')

        #loop through genes
        for gene in all_genes_dict[genome]:
            row_result = []
            gene_name = gene.name_id.split("__")[1]
            row_result.append(str(gene_name))

            for other_genome in other_genome_names:
                if gene_name not in blast_output[genome].keys() or other_genome not in blast_output[genome][gene_name].keys():
                    row_result.append('')
                else:
                    hit_result = 0.0
                    for hit_gene in blast_output[genome][gene_name][other_genome]:
                        if hit_gene[1] > hit_result:
                            hit_result = hit_gene[1]
                    '''
                    hit_result = ''
                    for hit_gene in blast_parsed[genome][gene_name][other_genome]:
                        hit_gene = [str(x) for x in hit_gene]
                        hit_result += '_'.join(hit_gene) + "|||"
                    '''
                    row_result.append(str(hit_result))
            #Write row
            outfile.write('\t'.join(row_result) + '\n')
##########################################################

###########Write hit matrices for each genome
#nt hit matrices
for genome in gbk_files_dict.keys():
    write_hit_matrix(genome, nt_blast_parsed, "hit_matrix/nt")
#aa hit matrices
for genome in gbk_files_dict.keys():
    write_hit_matrix(genome, aa_blast_parsed, "hit_matrix/aa")

######################################################################################