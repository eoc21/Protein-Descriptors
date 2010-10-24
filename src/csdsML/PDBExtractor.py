#A script for CSDS to download proteins from the PDB
# using a list of protein IDs stored in a file on her LINUX machine
#This program was written on a Windows machine!
#04/11/08

import urllib, sys, os

def extractPDBFiles(fileName):
    pdbFiles = open(fileName,'r')
    values = str(fileName).split("/")
    outputDirectory = ""
    if "1" in values[len(values)-1]:
        outputDirectory = "/PDZCluster1/"
    elif "2" in values[len(values)-1]:
        outputDirectory = "/PDZCluster2/"
    elif "3" in values[len(values)-1]:
        outputDirectory = "/PDZCluster3/"
    else:
        outputDirectory = "/PDZCluster4/"
    pdbIds =  pdbFiles.readlines()
    for i in range(len(pdbIds)):
        apdb = str(pdbIds[i]).strip()
        PDBfile = "http://www.rcsb.org/pdb/files/"+apdb+".pdb"
        datasource = urllib.urlopen(PDBfile)
        DS = datasource.readlines()
        pathname = os.path.dirname(sys.argv[1])
        val=pathname+outputDirectory+apdb+".pdb"
        f=open(val, 'w')
        for i in range(len(DS)):
            A =  DS[i]
            if(A[0:6] !='ANISOU'):
                B = str(A)
                f.write(B)    
        datasource.close()
extractPDBFiles(sys.argv[1])
extractPDBFiles(sys.argv[2])
extractPDBFiles(sys.argv[3])
extractPDBFiles(sys.argv[4])