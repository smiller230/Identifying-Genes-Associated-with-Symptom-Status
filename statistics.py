from Bio import SeqIO
import csv
import os
import pandas as pd
from scipy.stats import kendalltau

def transpose(file):
    m = pd.read_csv(file)
    rez = m.transpose()
    new_header = rez.iloc[0]
    rez =rez[1:]
    rez.columns =new_header
    transposed = open("transposed.csv","w")

    writer = csv.writer(transposed)
    rez.to_csv(transposed,index=False,header=True)



def kendall(names):
    matrix = pd.read_csv("transposed.csv") #open pres_abs matrix
    correlations = {}
    results1 = open("correlation_results.txt","w")
    results2 = open("correlation_matrix.csv","w")
    for i in names:
        x = matrix[i].to_list()
        y = matrix["Strain Symptom"].to_list()
        corr,_= kendalltau(x,y)
        if str(corr) == "nan":
            correlations[i]= 0
        else:
            correlations[i]= corr
    v=list(correlations.values())
    k=list(correlations.keys())
    positive = k[v.index(max(v))]
    negative = k[v.index(min(v))]
    results1.write("The cluster with the highest positive correlation to symptom status is " + positive + "\n")
    results1.write("The cluster with the highest negative correlation to symptom status is " + negative)
    writer = csv.writer(results2)
    for key,value in correlations.items():
        writer.writerow([key,value])
    results2.close()
    return correlations
cluster_values = []
for cluster_filename in os.listdir("./cluster_dir"):
    cluster_values.append(cluster_filename) #getting a list o
    
transpose("cluster_pres_abs_matrix.csv")
kendall(cluster_values)
