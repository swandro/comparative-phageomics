from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import Seq
import glob
import subprocess
from Bio import SeqIO
import os


#Search through all phage genomes and see how many are similar to ours.
genome_paths = glob.glob("all_phages/phage_genebank.gb")
#10,000 genomes

with open("all_phages/phage_genbank.gb",  'r') as genbank_file:
    for record in SeqIO.parse(genbank_file, "genbank"):
        fasta = open("temp/temp.fasta",'w')
        #######Add sequence to genes############
        i=1
        for gene in record.features:
            if gene.type in ["CDS", "tRNA"]:  ##only grabbing CDS and tRNA for now
                gene.seq = gene.extract(record.seq)  ##Add the sequence to each gene
                fasta.write(">{}_{}\n".format(record.name,str(i)))
                fasta.write(str(gene.seq) + '\n')
                i += 1
        fasta.close()
        nt_blast = NcbiblastnCommandline(cmd="blastn.exe", task='dc-megablast', out="temp/{}_blast.txt".format(record.name),
                                         outfmt='"6 qseqid sseqid pident qlen length mismatch gapope evalue bitscore qcovhsp"',
                                         query="temp/temp.fasta", db="blast_db/nt_db")
        nt_blast()




#Search through phagesdb genomes and see how many are similar to ours.
#1,900 genomes

with open("phagesdb/phages_db.gb",  'r') as genbank_file:
    for record in SeqIO.parse(genbank_file, "genbank"):
        fasta = open("phagesdb/temp.fasta",'w')
        #######Add sequence to genes############
        i=1
        for gene in record.features:
            if gene.type in ["CDS", "tRNA"]:  ##only grabbing CDS and tRNA for now
                gene.seq = gene.extract(record.seq)  ##Add the sequence to each gene
                fasta.write(">{}_{}\n".format(record.name,str(i)))
                fasta.write(str(gene.seq) + '\n')
                i += 1
        fasta.close()
        blastpath = "phagesdb/blast/{}_blast.txt".format(record.name)
        nt_blast = NcbiblastnCommandline(cmd="blastn.exe", task='dc-megablast', out=blastpath,
                                         outfmt='"6 qseqid sseqid pident qlen length mismatch gapope evalue bitscore qcovhsp"',
                                         query="phagesdb/temp.fasta", db="blast_db/nt_db")
        nt_blast()
        #see if it has hits
        empty = False
        with open("phagesdb/blast/" + str(record.name) + "_blast.txt", 'r') as hitfile:
            if len(hitfile.readlines()) < 1:
                empty = True
            else:
                gbpath = "phagesdb/genomes/{}.gb".format(record.name)
                SeqIO.write(record,gbpath, "genbank" )
        if empty:
            os.remove(blastpath)
        os.remove("phagesdb/temp.fasta")
