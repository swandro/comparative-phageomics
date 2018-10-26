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

####Create a master gene list######