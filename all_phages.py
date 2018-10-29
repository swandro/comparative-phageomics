from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import Seq
import glob
import subprocess
from Bio import SeqIO
import os


#Search through all phage genomes in Genbank (all_phages/phage_genbank.gb) and see how many are similar to ours.
genome_paths = glob.glob("all_phages/phage_genebank.gb")
#10,000 genomes
#Go through each genome and BLAST
all_gb_genes = []
with open("all_phages/phage_genbank.gb",  'r') as genbank_file:
    for record in SeqIO.parse(genbank_file, "genbank"):
        #fasta = open("temp/temp.fasta",'w')
        #######Add sequence to genes############
        i=1
        for gene in record.features:
            if gene.type in ["CDS", "tRNA"]:  ##only grabbing CDS and tRNA for now
                gene.seq = gene.extract(record.seq)  ##Add the sequence to each gene
                gene.uniq_id = ">{}__{}\n".format(record.name,str(i))
                all_gb_genes.append(gene) #Add genes to all gene list
                #fasta.write(">{}_{}\n".format(record.name,str(i)))
                #fasta.write(str(gene.seq) + '\n')
                i += 1
        #fasta.close()
        #nt_blast = NcbiblastnCommandline(cmd="blastn.exe", task='dc-megablast', out="temp/{}_blast.txt".format(record.name),
        #                                 outfmt='"6 qseqid sseqid pident qlen length mismatch gapope evalue bitscore qcovhsp"',
        #                                 query="temp/temp.fasta", db="blast_db/nt_db")
        #nt_blast()

#Parse through metadatafile
gb_phage_hits = []
with open("all_phages/phage_hits.txt",'r') as hits:
    for line in hits:
        gb_phage_hits.append(line.strip())
#Parse through all genbank files and keep ones that were on my hit list
gb_phages = dict()
with open("all_phages/phage_genbank.gb",  'r') as genbank_file:
    for record in SeqIO.parse(genbank_file, "genbank"):
        id = str(record.id).split('.')[0]
        if id in gb_phage_hits:
            gb_phages[id] = record

#Parse blast
gb_nt_blast_parsed = parse_blast("all_phages/blast/merged_blast.txt", 50,50)

#Merge the blast results
all_blast_parsed = dict()
all_blast_parsed.update(gb_nt_blast_parsed)
all_blast_parsed.update(nt_blast_parsed)

#get list of all genes in new phages
for gb in gb_phages:
    for gene in gb.features
        gb_phages["NC_029026"].features[0].__dict__

#make new master gene list
relevant_gb_genes = []
for gene in all_gb_genes:
    name = gene.uniq_id.split("__")[0][1:]
    if name in gb_phage_hits:
        relevant_gb_genes.append(gene)
relevant_gb_genes.append(all_genes)

extended_nt_master_genes = make_master_gene_list(relevant_gb_genes, all_blast_parsed)
for genome in gb_phages:
    print(genome)
    break



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
