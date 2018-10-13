from Bio import SeqIO
from Bio import AlignIO
from Bio import Seq
from Bio.Blast.Applications import NcbiblastxCommandline
import sys
import os
import pathlib
import glob
import numpy

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
        for redord in SeqIO.parse(genbank_file, "genbank"):
            gbk_files.append(redord)


#Get a list of all genes
all_genes = []
for genome in gbk_files:
    for gene in genome.features:
        if gene.type in ["CDS", "tRNA"]:   ##only grabbing CDS and tRNA for now
            gene.seq = gene.extract(genome.seq)  ##Add the sequence to each gene
            all_genes.append(gene)  ##Add gene to list

