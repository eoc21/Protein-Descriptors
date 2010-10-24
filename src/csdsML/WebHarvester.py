'''
Created on May 12, 2010
Harvests all PDB structures from PDB database
@author: ed
'''
import urllib, sys, os, random, math

otherProteins = open(sys.argv[1],'r')
otherProts = otherProteins.readlines()
notPDZDomain = []
pdzDomain =[]
pdzIds = []
notPDZIds = []
pdzActives = open(sys.argv[2],'w')
pdzInactives = open(sys.argv[3],'w')
for i in range(len(otherProts)):
    anotherProtein = str(otherProts[i]).strip()
    if ';'  in anotherProtein and 'PDZ' in anotherProtein:
        proteinIds = anotherProtein.split(';')
        pdzIds.append(str(proteinIds[0]).strip())
        result =  str(proteinIds[0]).strip()+";"+str(proteinIds[1]).strip()+";"+str(proteinIds[2]).strip()+"\n"
    elif ';' in anotherProtein and 'PDZ' not in anotherProtein and anotherProtein !="":
        proteinIdsX =  anotherProtein.split(';')
        notPDZIds.append(str(proteinIdsX[0]).strip())
        result = str(proteinIdsX[0]).strip()+";"+str(proteinIdsX[1]).strip()+";"+str(proteinIdsX[2]).strip()+"\n"


random.shuffle(notPDZIds)
sizeofInactives = int(math.floor(len(notPDZIds)/10))
tenPercent = notPDZIds[:sizeofInactives]
#Now for harvesting
outputDirectory="/InactiveProteins/"
for i in range(len(tenPercent)):
    id = str(tenPercent[i]).strip()
    PDBfile = "http://www.rcsb.org/pdb/files/"+id+".pdb"
    datasource = urllib.urlopen(PDBfile)
    DS = datasource.readlines()
    pathname = os.path.dirname(sys.argv[4])
    val=pathname+outputDirectory+id+".pdb"
    f=open(val, 'w')
    for i in range(len(DS)):
        A =  DS[i]
        if(A[0:6] !='ANISOU'):
            B = str(A)
            f.write(B)    
        datasource.close()