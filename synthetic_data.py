#1. create a synthetic genomaker definition with 5M and 10M variants 
#(only properties needed is the chromosome, location, position, index and cn). 
#cn has the following structure: zzg_chr_index so for example for marker in index 7 of chr 10 you would have cn = zzg_m_10_7
#2. Create bit array for each chromosomes with UKB participant IDs (ask Manos for the UKB participants). 
#The collections are named chrID e.g. chr1
#3. Ingest in new demo server.
import os
import sys
import json 
import random 
from multiprocessing import Process
import numpy as np 

def createGMdef(numbervars):
    data={}
    data['genomarkers']=[]
    chrsList = range(1,25)
    chrindex={}
    with open('data_'+str(numbervars)+'.json', 'w') as outfile:

        for i in range(1, numbervars):
            if i%100000==0:
                print(i)
            chr=random.choice(chrsList)
            position=random.choice(range(1,250000000))
            location=str(chr)+":"+str(position)

            while(location in data.values()):
                position=random.choice(range(1,250000000))
                location=str(chr)+":"+str(position)
            index=0
            if chr in chrindex:
                index=chrindex[chr]+1
            else:
                index=1
            chrindex[chr]=index
            cn="zzg_m_"+str(chr)+"_"+str(index)
            data['genomarkers'].append({'index':str(index),'Chromosome':str(chr),'Location':location,'cn':cn})
            outfile.write('{\"index\":\"'+str(index)+'\",\"Chromosome\":\"'+str(chr)+'\",\"Location\":\"'+location+'\",\"cn\":\"'+cn+'\"}\n')
            #json.dump(data, outfile)
    return data


def createBA(data,chr,id_conversion,subseti):



    print("Empezamos chr "+str(chr))
    filePosToWrite = open(str(chr)+"_"+str(len(data['genomarkers']))+"_"+str(subseti)+"_toBA.txt","w")
    markerchrpos={}
    for datagm in data['genomarkers']:
        if datagm['Chromosome']==str(chr):
            markerchrpos[datagm['Location']]=datagm['index']
    for id in id_conversion:
        
        string= ''.join(str(random.choice(range(0,2))) for _ in range(len(markerchrpos)+1))
        #for i in range(1,len(markerchrpos)+1):
        #    string=string+str(random.choice(range(0,2)))
        filePosToWrite.write(id+"\t"+string+"\n")









def main(argv):
    
    data5M=createGMdef(5000001)
    data10M=createGMdef(10000001)

    idconversion=[]
    with open('participants_id.csv') as f:
        print("Lectura de participants")
        lines=f.readlines()
        for line in lines:
            lineparts=line.strip().split('\t')
            idconversion.append(lineparts[0])
    subsets=np.array_split(idconversion,31)
    thread=[]


    for i in range(1,4):
        subseti=0
        for subset in subsets:
            subseti=subseti+1
            process=Process(target=createBA,args=(data5M,i,subset,subseti))
            thread.append(process)
        #process=Process(target=createBA,args=(data10M,i))
        #thread.append(process)
          
    for process in thread:
        process.start()
    for process in thread:
        process.join()



if __name__ == "__main__":
   main(sys.argv[1:])