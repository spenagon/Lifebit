import hail as hl 
from threading import Thread  
import os
import sys
from collections import deque
import time 
import csv 
from multiprocessing import Process,Condition,Queue
from hail.utils.java import FatalError







def an_item_is_available(queue):
    if queue.empty():
        print("an item is NOT available "+str(queue))
        return False
    else:
        print("an item is available "+str(queue))
        return True
    

def get_an_available_item(queue):
    print("get an available item "+str(queue))
    return queue.get()
    

def make_an_item_available(queue,item):
    print("make an item available "+str(queue)+" "+item)
    queue.put(item)
    
def copyS3ToEc2(cond1,q,myListFile,inputDir,outputDir,qaws_size):
    file1 = open(myListFile, 'r') 
    Lines = file1.readlines() 
    for line in Lines: 
        file=line.strip()
        if not line.startswith("#"):
            commandcp = "aws s3 cp " + file + " " + inputDir+" --no-progress"
            
            os.system(commandcp)
            #q.put(file)
            cond1.acquire()
            print("Thread aws make item available "+file)
            make_an_item_available(q,file)
            qaws_size=qaws_size+1
            cond1.notify_all()
            cond1.release()
            if qaws_size>3:
                time.sleep(1200)
    time.sleep(300)
    cond1.acquire()
    print("Thread aws make END available")
    make_an_item_available(q,"END")
    cond1.notify_all()
    cond1.release() 

def cmdsubprocess(fileName,subpart,inputDir,outputDir):

    idconversion={}
    with open('id_conversion') as f:
        #print("Lectura de id_conversion")
        lines=f.readlines()
        for line in lines:
            
            lineparts=line.strip().split('\t')
            #print(lineparts[1]+"\t"+lineparts[0])
            idconversion[lineparts[1]]=lineparts[0]
            
    
    chrName=fileName.split("_")[-3]
    #Load markers file 
    #chr10:81984:AA:22
    print("Thread cmd: Loading markers file")
    markerspos={}
    markersall={}
    with open('800k_to_extract_indexed.txt') as csvfile:
        markersreader = csv.reader(csvfile, delimiter=':')
        i=0
        for row in markersreader:
            if row[0]==chrName:
                #print("Thread cmd: markerspos "+row[0]+":"+row[1])
                markerspos[row[0]+":"+row[1]]=i
                i=i+1

    with open("genomarkers_gt_7.csv") as csvfile:
        markersreader = csv.reader(csvfile, delimiter=':')
        i=0
        for row in markersreader:
            if row[0]==chrName:
                #print("Thread cmd: markersall "+row[0]+":"+row[1]+":"+row[2])
                markersall[row[0]+":"+row[1]+":"+row[2]]=i
                i=i+1
    

    print("Thread cmdsubprocess "+str(subpart)+ " analysis file :"+fileName+"_"+str(subpart)+" of chr "+chrName)

    tablefile_subpart=outputDir+"/"+fileName+"_"+str(subpart)+".tsv"
    
    toCreateFile=False
    samplesStrings={}
    samplesAllStrings={}
    samples={}

    #First check if export files exist
    if os.path.isfile(outputDir+"/"+"extracted_"+chrName):
    #If exists load strings in dict
        print("Cargamos archivos de strings")
        with open(outputDir+"/"+"extracted_"+chrName) as f:
            lines=f.readlines()
            for line in lines:
                lineparts=line.strip().split("\t")
                #print("\tCarga "+lineparts[0]+" "+lineparts[1])
                samplesStrings[lineparts[0]]=lineparts[1]
        with open(outputDir+"/"+"extracted_alleles_"+chrName) as f:
            lines=f.readlines()
            for line in lines:
                lineparts=line.strip().split("\t")
                #print("\tCarga All "+lineparts[0]+" "+lineparts[1])
                samplesAllStrings[lineparts[0]]=lineparts[1]

    else:
        #print("Es necesario crear archivos de strings para este cromosoma")
        toCreateFile=True
        
    #Create files for position and alleles
        #Create string for all samples with length markers and idconversion
    with open(tablefile_subpart) as f:
        #print("Comenzamos lectura de tsv")
        lines=f.readlines()
        count=0
        for line in lines:
        #First line is sample names
            lineparts=line.strip().split("\t")
            if count==0:
                #print("Leemos primera linea")
                for i in range(1,len(lineparts)):
                    sampleName=lineparts[i]
                    #print("sampleName: "+sampleName)
                    sampleNameConv=idconversion[sampleName]
                    #print("sampleNameConv: "+sampleNameConv)
                    
                    samples[i]=sampleNameConv
                    if toCreateFile:
                        #print("Creamos strings")
                        string_pos = "0" * len(markerspos)
                        string_all="0" * len(markersall)
                        #print(string_pos)
                        #print(string_all)
                        samplesStrings[sampleNameConv]=string_pos
                        samplesAllStrings[sampleNameConv]=string_all
                    
                    
                count=count+1
            else: 
                count=count+1
                if count%50==0:
                    print("Subpart "+str(subpart)+" analyzed "+str(count))
                #print("Leemos linea genotipos")               
                        #chr10:129040748:A:G
                varLocParts=lineparts[0].split(":")
                                    
                #print(varLocParts)
                #print("Posicion a cambiar en string pos")
                if varLocParts[0]+":"+varLocParts[1]+":"+varLocParts[2]+varLocParts[3] in markersall:
                    #print(markerspos[varLocParts[0]+":"+varLocParts[1]])

                    for i in range(1,len(lineparts)):
                        gt=lineparts[i]
                        #print("Muestra "+samples[i])
                        #print("Genotipo "+gt)
                        #print("String pos antes ")
                        #print(samplesStrings[samples[i]])

                        #Position
                        #if gt=="1" or gt=="2":
                        if gt=="0/1" or gt=="1/1" or gt=="1" or gt=="1/.":
                            ss=samplesStrings[samples[i]]
                            posToChange=markerspos[varLocParts[0]+":"+varLocParts[1]]
                            samplesStrings[samples[i]]=ss[:posToChange]+"1"+ss[posToChange+1:]
                        #Alleles
                        #chr10:129040748:A:G
                        ss=samplesAllStrings[samples[i]]
                        #if gt=="1":
                        if gt=="0/1":
                            posToChange=markersall[varLocParts[0]+":"+varLocParts[1]+":"+varLocParts[2]+varLocParts[3]]
                            samplesAllStrings[samples[i]]=ss[:posToChange]+"1"+ss[posToChange+1:]
                        #elif gt=="2":
                        elif gt=="1/1":
                            posToChange=markersall[varLocParts[0]+":"+varLocParts[1]+":"+varLocParts[3]+varLocParts[3]]
                            samplesAllStrings[samples[i]]=ss[:posToChange]+"1"+ss[posToChange+1:]
                        #elif gt=="0":
                        elif gt=="0/0":
                            posToChange=markersall[varLocParts[0]+":"+varLocParts[1]+":"+varLocParts[2]+varLocParts[2]]
                            samplesAllStrings[samples[i]]=ss[:posToChange]+"1"+ss[posToChange+1:]
                        elif gt=="1":
                            if varLocParts[0]+":"+varLocParts[1]+":"+varLocParts[3] in markersall:
                                posToChange=markersall[varLocParts[0]+":"+varLocParts[1]+":"+varLocParts[3]]
                                samplesAllStrings[samples[i]]=ss[:posToChange]+"1"+ss[posToChange+1:]
                        elif gt=="0":
                            if varLocParts[0]+":"+varLocParts[1]+":"+varLocParts[2] in markersall:
                                posToChange=markersall[varLocParts[0]+":"+varLocParts[1]+":"+varLocParts[2]]
                                samplesAllStrings[samples[i]]=ss[:posToChange]+"1"+ss[posToChange+1:]
                        elif gt=="0/.":
                            if varLocParts[0]+":"+varLocParts[1]+":"+varLocParts[2]+"." in markersall:
                                posToChange=markersall[varLocParts[0]+":"+varLocParts[2]+":"+"."]
                                samplesAllStrings[samples[i]]=ss[:posToChange]+"1"+ss[posToChange+1:]
                        elif gt=="1/.":
                            if varLocParts[0]+":"+varLocParts[1]+":"+varLocParts[3]+"." in markersall:
                                posToChange=markersall[varLocParts[0]+":"+varLocParts[3]+":"+"."]
                                samplesAllStrings[samples[i]]=ss[:posToChange]+"1"+ss[posToChange+1:]


                        #print("String pos despues")
                        #print(samplesStrings[samples[i]])
        #Escribimos ficheros a disco
        #print("Escribimos archivos a disco")
        filePosToWrite = open(outputDir+"/"+"extracted_"+chrName+"_"+str(subpart),"w")
        fileAllToWrite = open(outputDir+"/"+"extracted_alleles_"+chrName+"_"+str(subpart),"w")
        for key in samples:
            sampleNameConv=samples[key]
            sampleString=samplesStrings[sampleNameConv]
            sampleAllString=samplesAllStrings[sampleNameConv]
            filePosToWrite.write(sampleNameConv+"\t"+sampleString+"\n")
            fileAllToWrite.write(sampleNameConv+"\t"+sampleAllString+"\n")
        filePosToWrite.close()
        fileAllToWrite.close()

def hailthread(cond1,q,cond2,qcm,inputDir,outputDir,qaws_size):

    #Load id_conversion file
    #table_idconv=hl.import_table('id_conversion')

    #Load markers files
    #table_makers_pos=hl.import_table('800k_to_extract_indexed2.txt',delimiter=':',no_header=True,impute=True)
    #table_markers_all=hl.import_table('800k_to_extract_indexed_alleles_gt2.txt',delimiter=':',no_header=True,impute=True)
    
    #cut -f 1 -d',' 800k_to_extract_indexed2.txt > interval_table
    #awk -F':' '{print $1"\t"$2"\t"$2}' interval_table > interval_table2

    hl.init()
    cond1.acquire()
    while not an_item_is_available(q):
        #print("Thread hail to sleep")
        #time.sleep(300)
        print("Thread hail to wait")

        cond1.wait()
    
    file=get_an_available_item(q)
    print("Thread hail get item "+file)
    qaws_size=qaws_size-1
    cond1.release()

    interval_table=hl.import_locus_intervals('interval_table2',reference_genome='GRCh38')

    while file != "END":
        
        
        fileParts=file.split("/")[-1]
        fileName=fileParts.replace(".vcf.gz","").replace(".gvcf.gz","")
        chrName=fileName.split("_")[-3]
        #myFNAL=fileName.split("\\.")
        #myTempId=myFNAL[0]
        #Load gVCF file
        #data=hl.import_vcf("/mnt/vol1/java/gel_test.vcf",force_bgz=True,reference_genome='GRCh38')
        #data=hl.import_vcf("/mnt/vol1/java/gel_mainProgramme_aggV2_chr10_129040437_131178399.vcf.gz",force_bgz=True,reference_genome='GRCh38')
        try:
            data=hl.import_vcf(inputDir+"/"+fileParts,force_bgz=True,reference_genome='GRCh38')        
            #data=hl.import_vcf(file.replace("s3://","s3a://"),force_bgz=True,reference_genome='GRCh38')
            #Filters PASS
            if chrName!="chrY":
                data = data.filter_rows(data.filters.size()>0, keep=False)
            #Multiallelic
            data=hl.split_multi_hts(data)
            #Join with markers
            data_filtered=data.filter_rows(hl.is_defined(interval_table[data.locus]))
            #Replace with 0s and 1s
            #Export
            #data_filtered_annot=data_filtered.annotate_entries(output=(data_filtered.GT.is_het()|data_filtered.GT.is_hom_var()))
            #Cambiamos la key para que contenga los alelos y no simplemente la posición
            variant_key=data_filtered.key_rows_by(variant=hl.variant_str(data_filtered.locus,data_filtered.alleles))
            #Exportamos el campo creado anteriormente y la cuenta del número de alelos ALT 
            #variant_key.GT.n_alt_alleles().export(outputDir+"/"+fileName+".tsv")
            #Change to export GT and not only the number of alt alleles because sex chromosomes need to check GT
            variant_key.GT.export(outputDir+"/"+fileName+".tsv")

            #Extract INFO fields

            data=hl.import_vcf(inputDir+"/"+fileParts,force_bgz=True,reference_genome='GRCh38',drop_samples=True)	
            #Filters PASS
            if chrName!="chrY":
                data = data.filter_rows(data.filters.size()>0, keep=False)
            #Multiallelic
            data=hl.split_multi_hts(data)
            #Join with markers
            data_filtered=data.filter_rows(hl.is_defined(interval_table[data.locus]))
            
            if chrName!="chrY":
                data_sr=data_filtered.select_rows(data_filtered.info.medianDepthAll,
                data_filtered.info.medianDepthNonMiss,
                data_filtered.info.medianGQ,
                data_filtered.info.missingness,
                data_filtered.info.completeGTRatio,
                data_filtered.info.ABratio,
                data_filtered.info.MendelSite,
                data_filtered.info.AN,
                data_filtered.info.AC,
                data_filtered.info.AC_Hom,
                data_filtered.info.AC_Het)

            else:
                data_sr=data_filtered.select_rows(
                data_filtered.info.AN,
                data_filtered.info.AC,
                data_filtered.info.AC_Hom,
                data_filtered.info.AC_Het)

            ht=data_sr.make_table()
            ht.export(outputDir+"/"+fileName+"_INFO.tsv")
            os.system("sed -i 's/\[//g' "+outputDir+"/"+fileName+"_INFO.tsv")
            os.system("sed -i 's/]//g' "+outputDir+"/"+fileName+"_INFO.tsv")
            os.system("cat "+outputDir+"/"+fileName+"_INFO.tsv | grep -v locus "+" >> "+outputDir+"/INFO_"+chrName)
            
            cond2.acquire()
            print("Thread hail make item available "+fileName)
            make_an_item_available(qcm,file)
            cond2.notify_all()
            cond2.release() 
        except FatalError as e:
            print("Exception2 in file:"+file)
            os.system("rm "+inputDir+"/"+fileParts)
            
        except AssertionError as e:
            print("Exception3 in file:"+file)
            os.system("rm "+inputDir+"/"+fileParts)
            
        except Exception as e:
            print("Exception in file:"+file)
            os.system("rm "+inputDir+"/"+fileParts)
            
            #raise Exception 
        cond1.acquire()
        while not an_item_is_available(q):
            #print("Thread hail to sleep")
            #time.sleep(300)
            print("Thread hail to wait")
            cond1.wait()
        
        file=get_an_available_item(q)
        print("Thread hail get item "+file)
        qaws_size=qaws_size-1
        cond1.release()
    time.sleep(300)
    cond2.acquire()
    print("Thread hail make END available")
    make_an_item_available(qcm,"END")
    cond2.notify_all()
    cond2.release() 


def rmThread(cond3,qrm,inputDir,outputDir):
    cond3.acquire()
    while not an_item_is_available(qrm):
        #print("Thread rm to sleep")
        #time.sleep(300)
        print("Thread rm to wait")
        cond3.wait()
    
    file=get_an_available_item(qrm)
    print("Thread rm get item "+file)
    cond3.release()
    while file != "END":
        fileParts=file.split("/")[-1]
        fileName=fileParts.replace(".vcf.gz","").replace(".gvcf.gz","")
        
        os.system("rm "+inputDir+"/"+fileParts)
        os.system("rm "+outputDir+"/"+fileName+"*")
        cond3.acquire()
        while not an_item_is_available(qrm):
            #print("Thread rm to sleep")
            #time.sleep(300)
            print("Thread rm to wait")
            cond3.wait()
        file=get_an_available_item(qrm)
        print("Thread rm get item "+file)
        cond3.release()

def cmdThread(cond2,qcm,cond3,qrm,inputDir,outputDir):
    cond2.acquire()
    while not an_item_is_available(qcm):
        print("Thread cmd to wait")
        cond2.wait()
    file=get_an_available_item(qcm)
    print("Thread cm get item "+file)
    cond2.release()
    while file != "END":
        fileParts=file.split("/")[-1]
        fileName=fileParts.replace(".vcf.gz","").replace(".gvcf.gz","")
        chrName=fileName.split("_")[-3]
        tablefile=outputDir+"/"+fileName+".tsv"

        tablefile1=outputDir+"/"+fileName+"_1.tsv"
        tablefile2=outputDir+"/"+fileName+"_2.tsv"
        tablefile3=outputDir+"/"+fileName+"_3.tsv"
        tablefile4=outputDir+"/"+fileName+"_4.tsv"
        tablefile5=outputDir+"/"+fileName+"_5.tsv"
        tablefile6=outputDir+"/"+fileName+"_6.tsv"
        tablefile7=outputDir+"/"+fileName+"_7.tsv"
        tablefile8=outputDir+"/"+fileName+"_8.tsv"
        os.system("cut -f 1,2-9999 "+tablefile+" > "+tablefile1)
        os.system("cut -f 1,10000-19999 "+tablefile+" > "+tablefile2)
        os.system("cut -f 1,20000-29999 "+tablefile+" > "+tablefile3)
        os.system("cut -f 1,30000-39999 "+tablefile+" > "+tablefile4)
        os.system("cut -f 1,40000-49999 "+tablefile+" > "+tablefile5)
        os.system("cut -f 1,50000-59999 "+tablefile+" > "+tablefile6)
        os.system("cut -f 1,60000-69999 "+tablefile+" > "+tablefile7)
        os.system("cut -f 1,70000- "+tablefile+" > "+tablefile8)



        p1=Process(target=cmdsubprocess,args=(fileName,1,inputDir,outputDir))
        p1.start()
        p2=Process(target=cmdsubprocess,args=(fileName,2,inputDir,outputDir))
        p2.start()
        p3=Process(target=cmdsubprocess,args=(fileName,3,inputDir,outputDir))
        p3.start()
        p4=Process(target=cmdsubprocess,args=(fileName,4,inputDir,outputDir))
        p4.start()
        p5=Process(target=cmdsubprocess,args=(fileName,5,inputDir,outputDir))
        p5.start()
        p6=Process(target=cmdsubprocess,args=(fileName,6,inputDir,outputDir))
        p6.start()
        p7=Process(target=cmdsubprocess,args=(fileName,7,inputDir,outputDir))
        p7.start()
        p8=Process(target=cmdsubprocess,args=(fileName,8,inputDir,outputDir))
        p8.start()

        p1.join()
        p2.join()
        p3.join()
        p4.join()
        p5.join()
        p6.join()
        p7.join()
        p8.join()


        #Unimos las partes
        print("Concat") 
        result=os.system("cat "+outputDir+"/"+"extracted_"+chrName+"_1 > "+ outputDir+"/"+"extracted_"+chrName)
        print(result)
        result=os.system("cat "+outputDir+"/"+"extracted_"+chrName+"_2 >> "+ outputDir+"/"+"extracted_"+chrName)
        print(result)
        result=os.system("cat "+outputDir+"/"+"extracted_"+chrName+"_3 >> "+ outputDir+"/"+"extracted_"+chrName)
        print(result)
        result=os.system("cat "+outputDir+"/"+"extracted_"+chrName+"_4 >> "+ outputDir+"/"+"extracted_"+chrName)
        print(result)
        result=os.system("cat "+outputDir+"/"+"extracted_"+chrName+"_5 >> "+ outputDir+"/"+"extracted_"+chrName)
        print(result)
        result=os.system("cat "+outputDir+"/"+"extracted_"+chrName+"_6 >> "+ outputDir+"/"+"extracted_"+chrName)
        print(result)
        result=os.system("cat "+outputDir+"/"+"extracted_"+chrName+"_7 >> "+ outputDir+"/"+"extracted_"+chrName)
        print(result)
        result=os.system("cat "+outputDir+"/"+"extracted_"+chrName+"_8 >> "+ outputDir+"/"+"extracted_"+chrName)
        print(result)
        print("Concat") 
        result=os.system("cat "+outputDir+"/"+"extracted_alleles_"+chrName+"_1 > "+ outputDir+"/"+"extracted_alleles_"+chrName)
        print(result)
        result=os.system("cat "+outputDir+"/"+"extracted_alleles_"+chrName+"_2 >> "+ outputDir+"/"+"extracted_alleles_"+chrName)
        print(result)
        result=os.system("cat "+outputDir+"/"+"extracted_alleles_"+chrName+"_3 >> "+ outputDir+"/"+"extracted_alleles_"+chrName)
        print(result)
        result=os.system("cat "+outputDir+"/"+"extracted_alleles_"+chrName+"_4 >> "+ outputDir+"/"+"extracted_alleles_"+chrName)
        print(result)
        result=os.system("cat "+outputDir+"/"+"extracted_alleles_"+chrName+"_5 >> "+ outputDir+"/"+"extracted_alleles_"+chrName)
        print(result)
        result=os.system("cat "+outputDir+"/"+"extracted_alleles_"+chrName+"_6 >> "+ outputDir+"/"+"extracted_alleles_"+chrName)
        print(result)
        result=os.system("cat "+outputDir+"/"+"extracted_alleles_"+chrName+"_7 >> "+ outputDir+"/"+"extracted_alleles_"+chrName)
        print(result)
        result=os.system("cat "+outputDir+"/"+"extracted_alleles_"+chrName+"_8 >> "+ outputDir+"/"+"extracted_alleles_"+chrName)
        print(result)
        cond3.acquire()
        print("Thread cm make item available "+file)
        make_an_item_available(qrm,file)
        cond3.notify_all()
        cond3.release() 
        cond2.acquire()
        while not an_item_is_available(qcm):
            print("Thread cm to wait")

            cond2.wait()
        
        file=get_an_available_item(qcm)
        print("Thread cm get item "+file)
        cond2.release()
    cond3.acquire()
    print("Thread cm make END available")
    make_an_item_available(qrm,"END")
    cond3.notify_all()
    cond3.release() 
    


def main(argv):
    
    myListFile=argv[0]
    print(myListFile)
    inputDir=argv[1]
    print(inputDir)
    outputDir=argv[2]
    print(outputDir)


    cond1=Condition()
    cond2=Condition()
    cond3=Condition()

    q=Queue()
    qcm=Queue()
    qrm=Queue()
    qaws_size=0
#    Thread(target=copyS3ToEc2,args=(cond1,)).start()
#    Thread(target=hailthread,args=(cond1,cond2,)).start()
    paws=Process(target=copyS3ToEc2,args=(cond1,q,myListFile,inputDir,outputDir,qaws_size,))
    paws.start()
    phail=Process(target=hailthread,args=(cond1,q,cond2,qcm,inputDir,outputDir,qaws_size,))
    phail.start()

    pcm=Process(target=cmdThread,args=(cond2,qcm,cond3,qrm,inputDir,outputDir,))
    pcm.start()
    prm=Process(target=rmThread,args=(cond3,qrm,inputDir,outputDir,))
    prm.start()

    paws.join()
    print("paws joined")
    phail.join()
    print("phail joined")
    pcm.join()
    print("pcm joined")
    prm.join()
    print("prm joined")


if __name__ == "__main__":
   main(sys.argv[1:])