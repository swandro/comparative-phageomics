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
for genome in gbk_files:
    for gene in genome.features:
        if gene.type in ["CDS", "tRNA"]:   ##only grabbing CDS and tRNA for now
            gene.seq = gene.extract(genome.seq)  ##Add the sequence to each gene
            qualifier_list = gene.qualifiers.keys()
            gene.id = str(genome.name) + '__'
            if "db_xref" in qualifier_list:
                gene.id += str(gene.qualifiers["db_xref"][0])
            elif "locus_tag" in qualifier_list:
                gene.id += str(gene.qualifiers["locus_tag"][0])
            elif "ID" in qualifier_list:
                gene.id += str(gene.qualifiers["ID"][0])
            else:
                raise Exception("Could not find the dictionary key for a unique gene identifier for genome: {G}, gene: {gen}".format(G= genome, gen=gene))
            all_genes.append(gene)  ##Add gene to list

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
test = NcbiblastnCommandline(cmd="blastn.exe", out="test_blast.txt", outfmt=6 , query= "db.fasta", db="blastDB")
test()

def parse_blast(file:"filename of blast result"):
    '''Takes in a BLAST result in format 6 and generates a dictionary: dict("Phage": dict(gene:[hit organism, hit gene, hit stats])'''

    result = dict()   ##Set up the result dictionary
    with open(file, 'r') as blast_table:
        for line in blast_table:
            ###Parse the line#############
            linesplit = line.strip().split()
            query_genome = linesplit[0].split("__")[0]
            query_gene = linesplit[0].split("__")[1]
            hit_genome = linesplit[1].split("__")[0]
            if hit_genome == query_genome:
                continue
            hit_gene = linesplit[1].split("__")[1]
            hit_percent = float(linesplit[2])
            evalue = linesplit[-1]
            #############################

            ###Add results to the dictionary#####
            if query_genome not in result.keys():
                result[query_genome] = dict()
            if query_gene not in result[query_genome].keys():
                result[query_genome][query_gene] = []
            result[query_genome][query_gene].append([hit_genome, hit_gene, hit_percent, evalue])
        return(result)

filename = "test_blast.txt"
blast_parsed = parse_blast(filename)
stats = dict()
for key in blast_parsed.keys():
    stats[key] = []
    for gene in blast_parsed[key].keys():
        stats[key].append(len(blast_parsed[key][gene]))

for key in stats.keys():
    print(key)
    print("mean: {}".format(sum(stats[key])/len(stats[key])))

sims = []
for key in blast_parsed.keys():
    for gene in blast_parsed[key].keys():
        for hit in blast_parsed[key][gene]:
            sims.append(hit[3])

sims = [float(x) for x in sims]