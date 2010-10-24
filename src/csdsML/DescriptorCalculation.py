'''
Created on May 13, 2010
Central  class to  calculate all  descriptors
@author: ed
'''
import os, sys
from PDBReader import *
from USRDescriptor import *
def stringFormatting(aList):
    unformatted = str(aList)
    unformatted1 = unformatted.replace('[', '')
    unformatted2 = unformatted1.replace(']','')
    unformatted3 = unformatted2.replace("'","").replace("\"","")
    return unformatted3      

def getUSRDescriptor(fileName):
    usrDescriptor = []
    pdbFile = open(fileName,'r')
    print fileName
    protein = pdbFile.readlines()
    aProtein = Protein()
    for i in range(len(protein)):
        aLine = str(protein[i]).strip().split()
        if aLine[0] == 'ATOM':
            if len(aLine) == 12:
                a = Atom(aLine[1],aLine[2],aLine[3],aLine[4],aLine[5],aLine[6],aLine[7],aLine[8],aLine[9],aLine[10],aLine[11])
                aProtein.addAtom(a)
            elif len(aLine) == 11:
                a = Atom(aLine[1],aLine[2],aLine[3],aLine[4],aLine[0],aLine[5],aLine[6],aLine[7],aLine[8],aLine[9],aLine[10])
                aProtein.addAtom(a)
    usr = USR()
    centroidPoint = usr.CalculateCentroid(aProtein)
    euclideanDistance = usr.EuclideanDistanceMeasure(aProtein,centroidPoint,0)
    v = usr.ClosestAndFurthestAtomToCentroid(aProtein)
    ffaIndex = usr.FurthestAtomFromFurthestAtom(aProtein,v[1])
    #From centroid
    meanFromCentroid = usr.MomentToCentroid(1,aProtein,centroidPoint)
    varianceFromCentroid = usr.MomentToCentroid(2,aProtein,centroidPoint)
    skewnessFromCentroid = usr.MomentToCentroid(3,aProtein,centroidPoint)
    #From closest atom
    m1CA = usr.MomentToX(1,aProtein,v[0].index)
    m2CA = usr.MomentToX(2,aProtein,v[0].index)
    m3CA = usr.MomentToX(3,aProtein,v[0].index)
    #FurthestAtom
    m1FA = usr.MomentToX(1,aProtein,v[1].index)
    m2FA = usr.MomentToX(2,aProtein,v[1].index)
    m3FA = usr.MomentToX(3,aProtein,v[1].index)
    #Furthest to Furthest atom
    m1FFA = usr.MomentToX(1,aProtein,ffaIndex)
    m2FFA = usr.MomentToX(2,aProtein,ffaIndex)
    m3FFA = usr.MomentToX(3,aProtein,ffaIndex)
    usrDescriptor.append(meanFromCentroid)
    usrDescriptor.append(varianceFromCentroid)
    usrDescriptor.append(skewnessFromCentroid)
    usrDescriptor.append(m1CA)
    usrDescriptor.append(m2CA)
    usrDescriptor.append(m3CA)
    usrDescriptor.append(m1FA)
    usrDescriptor.append(m2FA)
    usrDescriptor.append(m3FA)
    usrDescriptor.append(m1FFA)
    usrDescriptor.append(m2FFA)
    usrDescriptor.append(m3FFA)
    return usrDescriptor

#list all files in a directory
path = sys.argv[1]
output = open(sys.argv[2],'w')
errorFile = open(sys.argv[3],'w')
files = os.listdir(path)
counter = 0
errors=[]
for file in files:
    if file !='.pdb':
        errorCounter=0
        pdbFile = open(path+"/"+file,'r')
        try:
            USRDescriptorVector = getUSRDescriptor(path+"/"+file)
        except(ValueError,IndexError):
            errorCounter=1
            errors.append(file)
            #errorString = file+"\n"
            #errorFile.write(errorString)
            print "error in: "+file
 #       print USRDescriptorVector
        pdb = pdbFile.readlines()
        pdbIdentification = str(file).replace(".pdb", "").replace("'", "")
        protein = Protein()
        for i in range(len(pdb)):
            pdbLine = str(pdb[i]).strip()
            if "HEADER" in pdbLine:
                protein.header = pdbLine
            elif "TITLE" in pdbLine:
                protein.title =  pdbLine
            elif pdbLine.startswith("AUTHOR"):
                protein.author.append(pdbLine)
            elif pdbLine.startswith("REMARK"):
                protein.remark.append(pdbLine)
            elif pdbLine.startswith("SEQRES"):
                separatecomponents =  pdbLine.split()
                for i in range(len(separatecomponents)):
                    if i>3:
                        protein.sequence.append(separatecomponents[i])
        
        descriptor = ProteinDescriptors()
        aminoAcidComposition = stringFormatting(descriptor.calculateAminoAcidComposition(protein))
        dipeptideComposition = stringFormatting(descriptor.calculateDipeptideComposition(protein))
        moreauBrotoAutoCorrelation = stringFormatting(descriptor.calculateMoreauBrotoAutoCorrelationDescriptor(protein))
        moranAutoCorrelation = stringFormatting(descriptor.calculateMoranAutoCorrelation(protein))
        gearyAutoCorrelation = stringFormatting(descriptor.calculateGearyAutoCorrelation(protein))
        usrDescriptor = stringFormatting(USRDescriptorVector)
        #Add MD descriptors here.
        outputString =file+","+aminoAcidComposition+","+dipeptideComposition+","+moreauBrotoAutoCorrelation+","+moranAutoCorrelation+","+gearyAutoCorrelation+","+usrDescriptor+","+"InActive"+"\n"
        val = outputString.split(",")
        labels =[]
        for i in range(len(val)):
            labels.append("V"+str(i))
        formattedLabels = stringFormatting(labels)
        if counter == 0:
            output.write(formattedLabels+"\n")
        if errorCounter==0:
            output.write(outputString)
        print counter
        counter = counter+1
#Write errors
for i in range(len(errors)):
    anError =  str(errors[i]).strip()
    errorString = anError+"\n"
    errorFile.write(errorString)    
errorFile.close()   