from Bio import SeqIO
from Bio import AlignIO
from Bio import Seq
import glob
from collections import Counter

###Need to first run ImportGBK_BLAST.py to get variables

#####Parsing BLAST results######################################

#Function to do the work
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

##############################Parse the BLAST results###########
nt_blast_parsed = parse_blast("blast_output/nt_blast.txt", 50,50)
aa_blast_parsed = parse_blast("blast_output/aa_blast.txt", 50,50)
##################################################


########Create a master gene list####################

###Create a master gene list in which each entry is a gene and which gene it corresponds to in each genome####

###Helper functions
def make_master_gene_list(all_genes, parsed_blast_output):
    '''Function to make a master gene list'''
    result = dict()
    i=1
    for gene in all_genes:
        gene_found = False
        ##See if gene already exists in master_gene_list
        for master_gene in result.keys():
            if gene.origin in result[master_gene].keys():
                if gene.id_only in result[master_gene][gene.origin]:
                    gene_found = True
                    break
        if gene_found:
            continue
        same_genes = dict({gene.origin:[gene.id_only]})
        if gene.id_only in parsed_blast_output[gene.origin].keys():
            for phage in parsed_blast_output[gene.origin][gene.id_only].keys():
                if len(parsed_blast_output[gene.origin][gene.id_only][phage]) != 0:
                    if phage not in same_genes.keys():
                        same_genes[phage] = []
                    for hit in parsed_blast_output[gene.origin][gene.id_only][phage]:
                        same_genes[phage].append(hit[0])
        result[i] = same_genes
        i += 1
    return result
    #############################################################################


def find_gene(name, organism, master_gene_list):
    '''function to find a gene in the master gene list'''
    '''Finds a gene from a phage in the master gene list. Returns an index from the master gene list'''
    result = []
    for master_gene_i in master_gene_list.keys():
        if organism in master_gene_list[master_gene_i].keys():
            if name in master_gene_list[master_gene_i][organism]:
                result.append(master_gene_i)
    return result
    ####################################################

#Make aa and nt master gene list
nt_master_genes = make_master_gene_list(all_genes, nt_blast_parsed)
aa_master_genes = make_master_gene_list(all_genes, aa_blast_parsed)

##############Histogram

def gene_frequency(master_gene_list):
    counts = []
    for gene in master_gene_list:
        count = 0
        for genome in master_gene_list[gene]:
            count += len(master_gene_list[gene][genome])
        counts.append(count)
    counts.sort()
    return(Counter(counts))


sum(gene_frequency(nt_master_genes).values())