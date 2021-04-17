
#Python script to generate a presence absence matrix for each strain and the homologous genes that they either have or do not have 

from Bio import SeqIO
import csv
import os


def getList(): #getting the information about each strain from the input csv file
    strain_information = []
    with open("./strain_symptom_samples.csv", "r") as in_file: #grabbing the strain # and symptom from the csv file
        reader = csv.reader(in_file, delimiter = ',')
        i = False
        for row in in_file:
            if i != False: #avoiding the first "Strain # and Participant Symtom" heading
                strain_information.append(row.strip("\n").split(",")) #gives ["accession number", "strain_number", "Symtom status"] 
            i = True
    
    return strain_information
 

def getStrainGenes(strain_info_in):
    gene_dict = {}
    for i in range(len(strain_info_in)):
        strain_acc = strain_info_in[i][0]
        strain_num = strain_info_in[i][1]
        x = [] #list to append genes to
        records = list(SeqIO.parse("./ProteinSeqs/"+strain_acc+"_protein.faa","fasta"))
        for record in records:
            x.append(record.id)
        gene_dict[strain_num]=x #making a dictionary with the strain number as the key and the list of gene accession numbers as the value

    return gene_dict
 

def getClusterGenes(cluster_filename_in):
    cluster_genes = []
    with open(os.path.join("./cluster_dir/" + cluster_filename_in)) as cluster:
        for gene in SeqIO.parse(cluster, "fasta"):
             gene_name = gene.id
             cluster_genes.append(gene_name) #get a list of the genes in this one cluster
             #for computational purposes, this list is only created once and then compared to the genes in each strain 
    return cluster_genes
         
 
        
strain_info = getList() #getting the strain information from the "strain_symptom_samples" csv file

strain_gene_dict = getStrainGenes(strain_info) #dictionary containing the strains as keys and the gene accession numbers for each strain as values


cluster_matrix = open("cluster_pres_abs_matrix.csv", "w", newline = "") 
writer = csv.writer(cluster_matrix)

header = []
header.append("Cluster Number")
for strain in strain_info:
    header.append(strain[1])
    
writer.writerow(header)

for cluster_filename in os.listdir("./cluster_dir"):
    cluster_values = [] #setting up a list for the cluster row presence absence values
    cluster_values.append(cluster_filename) #getting a list of cluster names for the csv header
    cluster_genes = getClusterGenes(cluster_filename)   #generate the list of gene accession numbers present in this cluster 
    for strain in strain_gene_dict: #then go through each list of strain genes in order to see if there is a gene from that strain present in this cluster
        present_value = 0 #int value to determine if pres or abs
        for gene in strain_gene_dict[strain]: #for each gene within that strain's gene list 
            if gene in cluster_genes: #if that strain gene is in the cluster gene list
                present_value = 1 #make this "1" for present
                break #break from the for loop to reduce computational burden once gene found
        cluster_values.append(present_value)
    #write out the cluster value row to the csv file
    writer.writerow(cluster_values)
    
cluster_matrix.close() #this matrix will have to be transposed so that the strain numbers are the first column and each column is the cluster presence absence values

        
        
                
    