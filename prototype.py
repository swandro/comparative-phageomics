from Bio import SeqIO
from Bio import AlignIO
from Bio.Blast.Applications import NcbiblastxCommandline
import sys
import os
import pathlib
import glob
import numpy

###Code to run a blast
test = NcbiblastxCommandline(cmd="blastn.exe", out="test.txt", outfmt=5, query="seq.txt", db="db")
test()


#Importing genbank files
gbk_files = []
genome_file_list = glob.glob('genomes/*.gb')
for genome in genome_file_list:
    with open(genome,   'r') as genbank_file:
        temp = SeqIO.parse(genbank_file, "genbank")

for record in gbk_files:
    print(record)

    with open("cor6_6.gb", "rU") as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            print(record)