from Bio import SeqIO
from Bio import AlignIO
from Bio import Seq
from Bio.Blast.Applications import NcbiblastnCommandline
import sys
import os
import pathlib
import glob
import numpy
import subprocess

#######Code to run a blast#########
test = NcbiblastxCommandline(cmd="blastn.exe", out="test.txt", outfmt=5, query="seq.txt", db="db")
test()

'''
notes
a.seq = Seq object with sequence
a.features = list of SeqFeature
SeqFeature.qualifiers = dictionary of anotations
SeqFeature.qualifiers["db_xref"] or ["translation"] or ["product"]
'''



######Importing genbank files

#Get a list of all Genbank files
gbk_files = [] #make empty list
genome_file_list = glob.glob('genomes/*.gb')   ##Get names of all genomes

#Loop through and add genbank files to gbk_files
for genome in genome_file_list:
    with open(genome,   'r') as genbank_file:
        for record in SeqIO.parse(genbank_file, "genbank"):
            gbk_files.append(record)


#Get a list of all genes
all_genes = []
all_genes_dict = dict()
for genome in gbk_files:
    all_genes_dict[genome.name] = []
    for gene in genome.features:
        if gene.type in ["CDS", "tRNA"]:   ##only grabbing CDS and tRNA for now
            gene.seq = gene.extract(genome.seq)  ##Add the sequence to each gene
            qualifier_list = gene.qualifiers.keys()
            gene.id = str(genome.name) + '__'
            gene.origin = str(genome.name)
            if "db_xref" in qualifier_list:
                gene.id += str(gene.qualifiers["db_xref"][0])
                gene.id_only = str(gene.qualifiers["db_xref"][0])
            elif "locus_tag" in qualifier_list:
                gene.id += str(gene.qualifiers["locus_tag"][0])
                gene.id_only = str(gene.qualifiers["locus_tag"][0])
            elif "ID" in qualifier_list:
                gene.id += str(gene.qualifiers["ID"][0])
                gene.id_only = str(gene.qualifiers["ID"][0])
            else:
                raise Exception("Could not find the dictionary key for a unique gene identifier for genome: {G}, gene: {gen}".format(G= genome, gen=gene))
            all_genes.append(gene)  ##Add gene to list
            all_genes_dict[genome.name].append(gene)

#####Write genes to fasta file
db_fasta = ''
for gene in all_genes:
    db_fasta += ">" + gene.id + '\n' + str(gene.seq) + '\n'
outfile = open("db.fasta", 'w')
outfile.write(db_fasta)
outfile.close()
################

##Make BLASTdb of all genes
subprocess.call("makeblastdb.exe -in db.fasta -dbtype nucl -out blastDB", shell=True)
################


#Blast all genes
test = NcbiblastnCommandline(cmd="blastn.exe", task = 'dc-megablast', out="test_blast.txt", outfmt='"6 qseqid sseqid pident qlen length mismatch gapope evalue bitscore qcovhsp"' , query= "db.fasta", db="blastDB")
test()


def parse_blast(file:"filename of blast result", MIN_NUC_ID, MIN_Q_COV):
    '''Takes in a BLAST result in format 6 and generates a dictionary: dict("Phage": dict(gene:[hit organism, hit gene, hit stats])'''

    result = dict()   ##Set up the result dictionary
    with open(file, 'r') as blast_table:
        for line in blast_table:
            ###Parse the line#############
            linesplit = line.strip().split()
            query_genome = linesplit[0].split("__")[0]
            query_gene = linesplit[0].split("__")[1]
            hit_genome = linesplit[1].split("__")[0]

            hit_gene = linesplit[1].split("__")[1]
            hit_percent = float(linesplit[2])
            evalue = float(linesplit[-2])
            query_cov = float(linesplit[-1])
            #############################

            #Skip hit if does not meet provided cutoffs for nucleotide identity or query coverage
            if hit_percent < MIN_NUC_ID or query_cov < MIN_Q_COV:
                continue
            ###Add genome and gene to the dictionary if not already there#####
            if query_genome not in result.keys():
                result[query_genome] = dict()
            if query_gene not in result[query_genome].keys():
                result[query_genome][query_gene] = dict()
            if hit_genome not in result[query_genome][query_gene].keys():
                result[query_genome][query_gene][hit_genome]= []
            #Skip hit if hitting same genome##
            if hit_genome == query_genome:
                continue
            #Add passing hits to result
            result[query_genome][query_gene][hit_genome].append([hit_gene, hit_percent, query_cov, evalue])
        return(result)

filename = "test_blast.txt"
blast_parsed = parse_blast(filename, 40,40)

#Write output, one for each genome. Rows=genes, columns = other phages
for genome in all_genes_dict.keys():
    with open( (genome + "_hit_matrix.txt"), 'w') as outfile:
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
            gene_name = gene.id.split("__")[1]
            row_result.append(str(gene_name))

            for other_genome in other_genome_names:
                if gene_name not in blast_parsed[genome].keys() or other_genome not in blast_parsed[genome][gene_name].keys():
                    row_result.append('')
                else:
                    hit_result = 0.0
                    for hit_gene in blast_parsed[genome][gene_name][other_genome]:
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


###Calculate the most conserved genes####
master_gene_list = dict()
i=1
for gene in all_genes:
    gene_found = False
    ##See if gene already exists in master_gene_list
    for master_gene in master_gene_list.keys():
        if gene.origin in master_gene_list[master_gene].keys():
            if gene.id_only in master_gene_list[master_gene][gene.origin]:
                gene_found = True
                break
    if gene_found:
        continue
    same_genes = dict({gene.origin:[gene.id_only]})
    if gene.id_only in blast_parsed[gene.origin].keys():
        for phage in blast_parsed[gene.origin][gene.id_only].keys():
            if len(blast_parsed[gene.origin][gene.id_only][phage]) != 0:
                if phage not in same_genes.keys():
                    same_genes[phage] = []
                for hit in blast_parsed[gene.origin][gene.id_only][phage]:
                    same_genes[phage].append(hit[0])
    master_gene_list[i] = same_genes
    i += 1

def find_gene(name, organism):
    '''Finds a gene from a phage in the master gene list. Returns an index from the master gene list'''
    result = []
    for master_gene_i in master_gene_list.keys():
        if organism in master_gene_list[master_gene_i].keys():
            if name in master_gene_list[master_gene_i][organism]:
                result.append(master_gene_i)
    return result

lens = []
find_gene("SEED:fig|6666666.336390.peg.202", "Car")
master_gene_list[369]
for G in master_gene_list.keys():
    lens.append(len(master_gene_list[G]))


zeros = 0
ones = 0
more = 0
for GENE in all_genes:
    num = 0
    for master_gene in master_gene_list.keys():
        for gene_L in master_gene_list[master_gene].values():
            if GENE.id_only in gene_L:
                num += 1
    if num == 1:
        ones += 1
    elif num == 0:
        zeros += 1
    else:
        more += 1
        print(GENE)

