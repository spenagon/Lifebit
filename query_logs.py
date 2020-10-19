
import glob
import sys 
import codecs
import re
import collections  

logdict={}


def parselog(file):
    with codecs.open(file,"r","utf-8") as f:
        lines=f.readlines()
        theme=""
        theme_string=""
        for line in lines:
            line=line.strip()
            if line.startswith("--"):
                if theme != "" and theme in logdict:
                    logdict[theme]=logdict[theme]+"\n"+file+"\n"+theme_string
                elif theme !="":
                    logdict[theme]=file+"\n"+theme_string
                
                theme=""
                theme_string=""
            elif line !="": 
                if theme=="":
                    theme=line
                else:
                    theme_string=theme_string+"\n"+line
        if theme != "" and theme in logdict:
            logdict[theme]=logdict[theme]+"\n"+file+"\n"+theme_string
        elif theme !="":
            logdict[theme]=file+"\n"+theme_string

def querylog():
    contbool=True
    logdict_s=collections.OrderedDict(sorted(logdict.items()))
    while contbool:
        for key in logdict_s.keys():
            if re.match("n?[1-9]", key):
                #print("Entra aqui")
                print(key)
        print("Opciones: exacto, start o end")
        myinput=input()
        if myinput=="end":
            contbool=False
        elif myinput=="exacto":
            myinput=input()
            if myinput in logdict:
                print(logdict[myinput])
            else:
                print("key not in logdict")
        
            print("\n continue")
            input()
        elif myinput=="start":
            myinput=input()
            for key in logdict_s.keys():
                #print(myinput+"\t"+key)
                if key.startswith(myinput):
                    print(key)
                    print(logdict_s[key])
                    print("**********")
            print("\ncontinue")
            input()

def main(argv):
    myDir=argv[0]

    mylist = [f for f in glob.glob(myDir+"/2020*")]
    for myfile in mylist:
        if "." not in myfile: 
            print(myfile)
            parselog(myfile)
    querylog()


if __name__ == "__main__":
   main(sys.argv[1:])
