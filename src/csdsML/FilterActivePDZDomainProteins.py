'''
Created on Oct 23, 2010
Simply filter actives and add labels and remove first column
@author: ed
'''
import os,sys
class ProteinWithClassValue():
    def __init__(self,ID,classValue):
        self.id = ID
        self.classVal = classValue
        
proteins=[]
dataLabels = open(sys.argv[1],'r')
dataFileDescriptors = open(sys.argv[2],'r')
outputFile = open(sys.argv[3],'w')
dataDescriptors = dataFileDescriptors.readlines()
dataFile = dataLabels.readlines()
for i in range(len(dataFile)):
    result =  str(dataFile[i]).strip().split(",")
    p = ProteinWithClassValue(result[0],result[1])
    proteins.append(p)
counter = 0
for i in range(len(dataDescriptors)):
   dataDescriptorVal = str(dataDescriptors[i]).strip().split(",")
   if i==0:
       v = str(dataDescriptorVal[1:-1])
       v1 = v.replace("[", "")
       v2 = v1.replace("]","")
       v3 = v2.replace("'","")
       outputFile.write(v3+"\n")
   #dataDescriptorVal[1:-1]
   descriptorResult =  str(dataDescriptors[i]).strip()
   for j in range(len(proteins)):
       if proteins[j].id  in descriptorResult:
           d = descriptorResult.split(",")
           descriptorR1 =  d[1:-1]
           descriptorR2 = str(descriptorR1).replace("[", "")
           descriptorR3 = str(descriptorR2).replace("]","")
           descriptorR4 = str(descriptorR3).replace("'", "")
           val = descriptorR4+","+proteins[j].classVal
           outputFile.write(val+"\n")
           counter = counter+1
