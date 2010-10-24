'''
Created on Dec 7, 2009
Simplo script to concatenate molecular dynmaic simulation descriptors into 1 file
@author: ed
'''
import os, sys, glob
class MDDescriptor():
    def __init__(self,proteinId,descriptorArray):
        self.id =  proteinId
        self.descriptors = descriptorArray
        self.concatenatedDescriptors =  []
        
class BFactorEnergy():
    def __init__(self,proteinId,bFactor):
        self.id = proteinId
        self.bFactor = bFactor

class BindingEnergy():
    def __init__(self,proteinId,ligandCrystalEnergy,etwvEnergy,gifvEnergy):
        self.id = proteinId
        self.ligandCrystalEnergy = ligandCrystalEnergy
        self.etwv = etwvEnergy
        self.gifvEnergy = gifvEnergy
        
bfactorEnergies = []
md157 = open(sys.argv[1],'r')
energyDescriptors = open(sys.argv[2],'r')
path = sys.argv[3]
outputFile = open(sys.argv[4],'w')
#Processes all the .all files
for infile in glob.glob( os.path.join(path, '*.all') ):
  data = open(infile,'r')
  descriptor = data.readlines()
  bFactor = str(descriptor[0]).split()
  splitUpInfile = str(infile).split("/")
  fileId = splitUpInfile[len(splitUpInfile)-1]
  processedId = fileId.replace("af1.all","")
  processedId1 = processedId.replace("bf-","")
  bfactorEnergy = BFactorEnergy(processedId1,bFactor[1])
  bfactorEnergies.append(bfactorEnergy)

#process energies
bindingEnergyDescriptors = []
bindingEnergies = energyDescriptors.readlines()
for i in range(len(bindingEnergies)):
    val = str(bindingEnergies[i]).split()
    if val[2] == 'NA':
        aBindingEnergy = BindingEnergy(val[0],val[1],"0",val[3])
        bindingEnergyDescriptors.append(aBindingEnergy)
    else:
        aBindingEnergy = BindingEnergy(val[0],val[1],val[2],val[3])
        bindingEnergyDescriptors.append(aBindingEnergy)

#processes md157 file
molecularDynamicDescriptors = md157.readlines()
MDDescriptors = []
concatenatedDescriptors = []
for i in range(len(molecularDynamicDescriptors)):
    descriptor=[]
    individualDescriptor = str(molecularDynamicDescriptors[i]).split(",")
    print individualDescriptor,len(individualDescriptor)
    for j in range(len(individualDescriptor)):
        aDescriptorValue = str(individualDescriptor[j]).strip()
        descriptor.append(aDescriptorValue)
    anMDDescriptor = MDDescriptor(str(individualDescriptor[0]),str(descriptor[1:]))
    MDDescriptors.append(anMDDescriptor)
#Concatenate the descriptors
mdAndBindingEnergyDescriptors =[]
for i in range(len(MDDescriptors)):
    for j in range(len(bindingEnergyDescriptors)):
        if MDDescriptors[i].id == bindingEnergyDescriptors[j].id:
            mdBE = str(MDDescriptors[i].id)+","+str(MDDescriptors[i].descriptors)+","+str(bindingEnergyDescriptors[j].ligandCrystalEnergy)+","+str(bindingEnergyDescriptors[j].etwv)+","+str(bindingEnergyDescriptors[j].gifvEnergy)
            mdBEProcessed = mdBE.replace("[","")
            mdBEProcessed1 = mdBEProcessed.replace("]","")
            mdAndBindingEnergyDescriptors.append(mdBEProcessed1)
           
#Concatenate B-Factor
for i in range(len(mdAndBindingEnergyDescriptors)):
    id = str(mdAndBindingEnergyDescriptors[i]).split(",")
    for j in range(len(bfactorEnergies)):
        if str(id[0]) == str(bfactorEnergies[j].id).lower():
            descriptorWithBFactor = str(mdAndBindingEnergyDescriptors[i])+","+str(bfactorEnergies[j].bFactor)
            descriptorWithBFactor = descriptorWithBFactor.replace("'","")
            outputFile.write(descriptorWithBFactor+"\n")

     
     
     