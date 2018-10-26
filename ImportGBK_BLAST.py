from Bio import SeqIO
from Bio import AlignIO
from Bio import Seq
from Bio.Blast.Applications import NcbiblastnCommandline
import glob
import subprocess

######################Importing and formatting data####################################

##Get filenames of genomes in the genbank file format
genome_file_list = glob.glob('prokka_genomes/*.gbf')   ##Get names of all genomes

#Create dictionary of imported Genbank files. Genome name is Key
gbk_files_dict = dict()
for genome in genome_file_list:
    with open(genome,   'r') as genbank_file:
        for record in SeqIO.parse(genbank_file, "genbank"):
            gbk_files_dict[record.name] = record


#Create list of all genes and dictionary of all genes by genome
##Only keep CDS and tRNA
all_genes = []
all_genes_dict = dict()

for genome_name in gbk_files_dict.keys():
    all_genes_dict[genome_name] = []
    genome = gbk_files_dict[genome_name]
    for gene in genome.features:
        if gene.type in ["CDS", "tRNA"]:   ##only grabbing CDS and tRNA for now
            gene.seq = gene.extract(genome.seq)  ##Add the sequence to each gene
            qualifier_list = gene.qualifiers.keys()
            gene.name_id = str(genome.name) + '__'
            gene.origin = str(genome.name)
            if "locus_tag" in qualifier_list:
                gene.name_id += str(gene.qualifiers["locus_tag"][0])
                gene.id_only = str(gene.qualifiers["locus_tag"][0])
            else:
                raise Exception("Could not find the dictionary key for a unique gene identifier for genome: {G}, gene: {gen}".format(G= genome, gen=gene))
            all_genes.append(gene)  ##Add gene to list
            all_genes_dict[genome_name].append(gene)

#################################################################

###############BLASTING genes##################################

##################Make blast DB###########
### Only needs to be run once
#Make nucleotide db
db_fna = ''
#Make amino acid db
db_faa = ''
#Make nucleotide db of just trna genes
db_trna = ''

for gene in all_genes:
    if gene.type == "CDS":
        db_faa += ">" + gene.name_id + '\n' + str(gene.seq.translate()) + '\n'
        db_fna += ">" + gene.name_id + '\n' + str(gene.seq) + '\n'
    elif gene.type == "tRNA":
        db_fna += ">" + gene.name_id + '\n' + str(gene.seq) + '\n'
        db_trna += ">" + gene.name_id + '\n' + str(gene.seq) + '\n'

with open("blast_db/db.fna", 'w') as outfile_fna:
    outfile_fna.write(db_fna)
with open("blast_db/db.faa", 'w') as outfile_faa:
    outfile_faa.write(db_faa)
with open("blast_db/db_tRNA.fna", 'w') as outfile_trna:
    outfile_trna.write(db_trna)

#Call BLAST to make blastDBs
subprocess.call("makeblastdb.exe -in blast_db/db.fna -dbtype nucl -out blast_db/nt_db", shell=True)
subprocess.call("makeblastdb.exe -in blast_db/db.faa -dbtype prot -out blast_db/aa_db", shell=True)
subprocess.call("makeblastdb.exe -in blast_db/db_tRNA.fna -dbtype nucl -out blast_db/trna_db", shell=True)

#Names of BLAST dbs
#blast_db/nt_db
#blast_db/aa_db
#blast_db/trna_db

#########################################

#BLAST the databases against themselves
nt_blast = NcbiblastnCommandline(cmd="blastn.exe", task = 'dc-megablast', out="blast_output/nt_blast.txt", outfmt='"6 qseqid sseqid pident qlen length mismatch gapope evalue bitscore qcovhsp"' , query= "blast_db/db.fna", db="blast_db/nt_db")
aa_blast = NcbiblastnCommandline(cmd="blastp.exe",  out="blast_output/aa_blast.txt", outfmt='"6 qseqid sseqid pident qlen length mismatch gapope evalue bitscore qcovhsp"' , query= "blast_db/db.faa", db="blast_db/aa_db")
trna_blast = NcbiblastnCommandline(cmd="blastn.exe", task = 'dc-megablast', out="blast_output/trna_blast.txt", outfmt='"6 qseqid sseqid pident qlen length mismatch gapope evalue bitscore qcovhsp"' , query= "blast_db/db_tRNA.fna", db="blast_db/trna_db")
nt_blast()
aa_blast()
trna_blast()

################################################################