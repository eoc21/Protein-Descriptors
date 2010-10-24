'''
Created on Nov 28, 2009
Script to calculate theoretical protein descriptors-
Amino acid composition, dipeptide composition,moreauBrotoAutoCorrelation,
moranAutoCorrelation & gearyAutoCorrelation.
2pzd or 2pdz? for one of the proteins.
@author: ed
'''
import os, sys
from PDBReader import *
from USRDescriptor import *

class MDDescriptor():
    def __init__(self,proteinId,descriptorVector):
        self.id = proteinId
        self.descriptors = descriptorVector

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
#Add molecular dynamic descriptors
mdDescriptorFile = open(sys.argv[3],'r')
molecularDynamicDescriptors = mdDescriptorFile.readlines()
mdDescriptors = []
for i in range(len(molecularDynamicDescriptors)):
    idAndDescriptor = str(molecularDynamicDescriptors[i]).split(",")
    descriptor = []
    for j in range(len(idAndDescriptor)):
        d = str(idAndDescriptor[j]).strip()
        descriptor.append(d) 
    md = MDDescriptor(str(descriptor[0]),str(descriptor[1:]))
    print md.id,len(str(md.descriptors).split(","))
    mdDescriptors.append(md)
output = open(sys.argv[2],'w')
#list all files in a directory
path = sys.argv[1]
files = os.listdir(path)
counter = 0
for file in files:
    if file !='.pdb':
        pdbFile = open(path+"/"+file,'r')
        USRDescriptorVector = getUSRDescriptor(path+"/"+file)
        print USRDescriptorVector
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
        #Add MD descriptors here.
        outputString =""
        for j in range(len(mdDescriptors)):
            if str(mdDescriptors[j].id) == pdbIdentification:
                mdDescriptor = str(mdDescriptors[j].descriptors).split(",")
                mdDescriptor = stringFormatting(mdDescriptor)
                outputString = aminoAcidComposition+","+dipeptideComposition+","+moreauBrotoAutoCorrelation+","+moranAutoCorrelation+","+gearyAutoCorrelation+","+"PDZ4"+"\n"
        val = outputString.split(",")
        labels =[]
        for i in range(len(val)):
            labels.append("V"+str(i))
        formattedLabels = stringFormatting(labels)
        if counter == 0:
            output.write(formattedLabels+"\n")
        output.write(outputString)
        counter = counter+1
        