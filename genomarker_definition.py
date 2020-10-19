#New genomarker definition

import os
import sys
import json 
import random 
from multiprocessing import Process
import numpy as np 
from pyspark.sql import SparkSession
from pyspark import SparkConf, SparkContext
import pyspark.sql.functions as f


def process(file,filemarkers,fileout):
    spark = SparkSession.builder.appName("GEL data pipeline").getOrCreate()
    spark.conf.set("spark.driver.memory", "300g")
    filemarkers = spark.read.load(filemarkers+'2', format="csv", sep="\t", inferSchema="true", header="true")
    annot = spark.read.load(file+'_toDF', format="csv", sep="\t", inferSchema="true", header="true")
    annot_dbnsfp_gene=spark.read.load('dbNSFP4.1_gene.complete', format="csv", sep="\t", inferSchema="true", header="true")
    annot=annot.join(annot_dbnsfp_gene,annot_dbnsfp_gene.Ensembl_gene==annot.Gene,'left')
    df = filemarkers.join(annot, ['Location'], 'left')
    dfcolumns2=df.columns
    dfcolumns2.remove("index")
    dfcolumns2.remove("id")
    dfcolumns2.remove("cn")
    dfcolumns2.remove("#Uploaded_variation")
    dfcolumns2.remove("Location")
    dfcolumns2.remove("Allele")
    
    dfcolumns=dfcolumns2.copy()
    dfcolumns.remove("Gene")


    df=df.groupBy("index","id","cn","#Uploaded_variation","Location","Allele","Gene").agg(*(f.concat_ws('|',f.collect_set(c)).alias(c) for c in dfcolumns))
    df=df.groupBy("index","id","cn","#Uploaded_variation","Location","Allele").agg(*(f.concat_ws(';',f.collect_list(c)).alias(c) for c in dfcolumns2))
    df.write.json(fileout)
    #pandas_df=df.toPandas()
    #pandas_df.to_json(fileout)

def preprocess(file,filemarkers):
    #Sacamos los headers 
    commandcp = "grep -E \"^#\" " + file + " | tail -n 1 > " + file+"_headers"
    os.system(commandcp)
    commandcp = "grep -Ev \"^#\" " + file + "  > " + file+"_wo_headers"
    os.system(commandcp)
    commandcp = "cat " + file + "_headers > " + file+"_toDF" 
    os.system(commandcp)
    commandcp = "cat " + file + "_wo_headers  >> " + file+"_toDF"
    os.system(commandcp)

    #Leemos el archivo de markers y añadimos index 
    filePosToWrite = open(filemarkers+"2","w")
    filePosToWrite.write("Location\tid\tindex\tcn\n")
    chrpos={}
    with open(filemarkers) as f:
        print("Lectura de markers file")
        lines=f.readlines()
        for line in lines:
            lineparts=line.strip().split(':')
            chr=lineparts[0]
            index=0
            if(chr in chrpos):
                index=chrpos[chr]+1
            else:
                index=1
            chrpos[chr]=index
            pos=lineparts[1]
            id=lineparts[2]
            filePosToWrite.write(str(chr).replace('chr','')+":"+str(pos)+"\t"+str(id)+"\t"+str(index)+"\tzzg_m_"+str(chr).replace("chr","")+"_"+str(index)+"\n")
    


def main(argv):
    #Leemos fichero de marcadores
    # #index=2 and id=446699
    # Añadir el index por chr al archivo 
    file=argv[0]
    filemarkers=argv[1]
    fileout=argv[2]
    preprocess(file,filemarkers)
    process(file,filemarkers,fileout)


if __name__ == "__main__":
   main(sys.argv[1:])

