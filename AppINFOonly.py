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
            #if qaws_size>3:
            #    time.sleep(1200)
    time.sleep(300)
    cond1.acquire()
    print("Thread aws make END available")
    make_an_item_available(q,"END")
    cond1.notify_all()
    cond1.release() 



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


            #Extract INFO fields
            
            

            data=hl.import_vcf(inputDir+"/"+fileParts,force_bgz=True,reference_genome='GRCh38',drop_samples=True)	
            #Filters PASS
            if chrName!="chrY":
                data = data.filter_rows(data.filters.size()>0, keep=False)
            #Multiallelic
            data=hl.split_multi_hts(data)
            #Join with markers
            data_filtered=data.filter_rows(hl.is_defined(interval_table[data.locus]))
            
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

            ht=data_sr.make_table()
            ht.export(outputDir+"/"+fileName+"_INFO.tsv")
            os.system("sed -i 's/\[//g' "+outputDir+"/"+fileName+"_INFO.tsv")
            os.system("sed -i 's/]//g' "+outputDir+"/"+fileName+"_INFO.tsv")
            os.system("cat "+outputDir+"/"+fileName+"_INFO.tsv | grep -v locus "+" >> "+outputDir+"/INFO_"+chrName)
            os.system("rm "+inputDir+"/"+fileParts)
            
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


    paws.join()
    print("paws joined")
    phail.join()
    print("phail joined")


if __name__ == "__main__":
   main(sys.argv[1:])