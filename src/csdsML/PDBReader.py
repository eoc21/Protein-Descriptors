'''
Created on Nov 28, 2009
Module to read protein databank files
@author: ed
'''
from math import *
class Atom():
    def __init__(self,index,atomType,aminoAcid,chain,val, x,y,z,A,anotherVal, atomT):
        self.index = index
        self.specificAtomType = atomType
        self.aminoAcid = aminoAcid
        self.chain = chain
        self.extraValue1 = val
        self.xCoordinate = x
        self.yCoordinate = y
        self.zCoordinate = z
        self.extraValue2 = A
        self.extraValue3 = anotherVal
        self.atomType = atomT
    
    def getXCoordinate(self):
        return self.xCoordinate
    
    def getYCoordinate(self):
        return self.yCoordinate
    
    def getZCoordinate(self):
        return self.zCoordinate

class Protein():
    def __init__(self):
        self.sequence = []
        self.header = ""
        self.title = ""
        self.models = []
        self.author =[]
        self.remark =[]
        self.singleAASequence = []
        self.aaDictionary ={}
        self.geneDataBankId = ""
        self.atoms = []
    def addAtom(self, atom):
        self.atoms.append(atom)
    
    def getAtom(self, i):
        return self.atoms[i]
    
    def getAtoms(self):
        return self.atoms
    
    def convertAAToSingleLetters(self):
        self.aaDictionary = {'GLY': 'G', 'ALA': 'A','LEU':'L','MET':'M','PHE':'F','TRP':'W','LYS':'K',
              'GLN':'Q','GLU':'E','SER':'S','PRO':'P','VAL':'V','ILE':'I','CYS':'C','TYR':'Y',
              'HIS':'H','ARG':'R','ASN':'N','ASP':'D','THR':'T'}
        for i in range(len(self.sequence)):
            if self.sequence[i] == 'TPO' or self.sequence[i]=='A' or self.sequence[i]=='C' or self.sequence[i]=='G' or self.sequence[i]=='SIN' or self.sequence[i]=='CH2' or self.sequence[i]=='BMT' or self.sequence[i]=='ABA' or self.sequence[i]=='KPH' or self.sequence[i]=='DPN' or self.sequence[i]=='U' or self.sequence[i]=="NH2" or self.sequence[i]=="FME" or self.sequence[i]=='UNK' or self.sequence[i]=='E' or self.sequence[i]=='DT'  or self.sequence[i] =='DA' or self.sequence[i]=='DG' or self.sequence[i]=='DC' or self.sequence[i]=='PCA' or self.sequence[i] =='MLY' or self.sequence[i] == 'MSE' or self.sequence[i] == 'MIS' or self.sequence[i] == 'SEP' or self.sequence[i] == 'TYI' or self.sequence[i] == 'PTR' or self.sequence[i] == 'SNN' or self.sequence[i]== 'ACE':
                pass
            if self.sequence[i] not in self.aaDictionary:
                pass
            else:
                aminoAcidSingleLetter = self.aaDictionary[self.sequence[i]]
                self.singleAASequence.append(aminoAcidSingleLetter)

class Model():
    def __init__(self,atoms):
        self.atoms = atoms

class ProteinDescriptors():
    def __init__(self):
        self.aaComposition = []
        self.dipeptideComposition = []
        self.hydrophobicityDictionary = {'A':0.02,'R':-0.42,'N':-0.77,'D':-1.04,'C': 0.77,'Q':-1.10, 
                                         'E':-1.14,'G':-0.80,'H':0.26,'I':1.81,'L':1.14,'K':-0.41,'M':1.00,
                                          'F':1.35,'P':-0.09,'S':-0.97,'T':-0.77,'W':1.71,'Y':1.11,'V':1.13}
    
        self.averageFlexibilityDictionary = {'A':0.357,'R':0.529,'N':0.463,'D':0.511,'C': 0.346,'Q':0.493, 
                                         'E':0.497,'G':0.544,'H':0.323,'I':0.462,'L':0.365,'K':0.466,'M':0.295,
                                          'F':0.314,'P':0.509,'S':0.507,'T':0.444,'W':0.305,'Y':0.420,'V':0.386}
    
        
        self.polarizabilityDictionary= {'A':0.046,'R':0.291,'N':0.134,'D':0.105,'C': 0.128,'Q':0.180, 
                                         'E':0.151,'G':0.000,'H':0.230,'I':0.186,'L':0.186,'K':0.219,'M':0.221,
                                          'F':0.290,'P':0.131,'S':0.062,'T':0.108,'W':0.409,'Y':0.298,'V':0.140}
        
        self.freeEnergyOfSolutionInWaterDictionary = {'A':-0.368,'R':-1.03,'N':0,'D':2.06,'C': 4.53,'Q':0.731, 
                                         'E':1.77,'G':-0.525,'H':0,'I':0.791,'L':1.07,'K':0,'M':0.656,
                                          'F':1.06,'P':-2.24,'S':-0.524,'T':0,'W':1.60,'Y':4.91,'V':0.401}
        
        self.residueASAInTripeptideDictionary = {'A':115,'R':225,'N':160,'D':150,'C':135,'Q':180, 
                                         'E':190,'G':75,'H':195,'I':175,'L':170,'K':200,'M':185,
                                          'F':210,'P':145,'S':115,'T':140,'W':255,'Y':230,'V':155}
        
        self.residueVolumeDictionary = {'A':52.6,'R':109.1,'N':75.7,'D':68.4,'C':68.3,'Q':89.7, 
                                         'E':84.7,'G':36.3,'H':91.9,'I':102.0,'L':102.0,'K':105.1,'M':97.7,
                                          'F':113.9,'P':73.6,'S':54.9,'T':71.2,'W':135.4,'Y':116.2,'V':85.1}
        
        self.stericParameterDictionary = {'A':0.52,'R':0.68,'N':0.76,'D':0.76,'C':0.62,'Q':0.68, 
                                         'E':0.68,'G':0.00,'H':0.70,'I':1.02,'L':0.98,'K':0.68,'M':0.78,
                                          'F':0.70,'P':0.36,'S':0.53,'T':0.50,'W':0.70,'Y':0.70,'V':0.76}
        
        self.relativeMutabilityDictionary = {'A':100,'R':65,'N':134,'D':106,'C':20,'Q':93, 
                                         'E':102,'G':49,'H':66,'I':96,'L':40,'K':56,'M':94,
                                          'F':41,'P':56,'S':120,'T':97,'W':18,'Y':41,'V':74}
        
        self.MoreauBrotoDescriptors = []
        self.MoranAutoCorrelation = []
        self.GearyAutoCorrelation = []
        
    def calculateAminoAcidComposition(self,protein):
        protein.convertAAToSingleLetters()
        sequenceSize = len(protein.sequence)
        alanineCounter = 0
        arginineCounter = 0
        asparagineCounter = 0
        asparticAcidCounter = 0
        cysteineCounter = 0
        glutamicAcidCounter = 0
        glutamineCounter = 0
        glycineCounter = 0
        histidineCounter = 0
        isoleucineCounter = 0
        leucineCounter = 0
        lysineCounter = 0
        methioneCounter = 0
        phenylalanineCounter = 0
        prolineCounter =  0
        serineCounter = 0
        threonineCounter = 0
        tryptophanCounter = 0
        tyrosineCounter = 0
        valineCounter = 0
        for i in range(len(protein.singleAASequence)):
            if protein.singleAASequence[i] == 'A':
                alanineCounter = alanineCounter+1
            elif protein.singleAASequence[i] == 'R':
                arginineCounter = arginineCounter+1
            elif protein.singleAASequence[i] == 'N':
                asparagineCounter = asparagineCounter+1
            elif protein.singleAASequence[i] == 'D':
                asparticAcidCounter = asparticAcidCounter+1
            elif protein.singleAASequence[i] == 'C':
                cysteineCounter = cysteineCounter+1
            elif protein.singleAASequence[i] == 'E':
                glutamicAcidCounter = glutamicAcidCounter+1
            elif protein.singleAASequence[i] == 'Q':
                glutamineCounter = glutamineCounter+1
            elif protein.singleAASequence[i] == 'G':
                glycineCounter = glycineCounter+1
            elif protein.singleAASequence[i] == 'H':
                histidineCounter =  histidineCounter+1
            elif protein.singleAASequence[i] == 'I':
                isoleucineCounter = isoleucineCounter+1
            elif protein.singleAASequence[i] == 'L':
                leucineCounter = leucineCounter+1
            elif protein.singleAASequence[i] == 'K':
                lysineCounter = lysineCounter+1
            elif protein.singleAASequence[i] == 'M':
                methioneCounter =  methioneCounter+1
            elif protein.singleAASequence[i] == 'F':
                phenylalanineCounter = phenylalanineCounter+1
            elif protein.singleAASequence[i] == 'P':
                prolineCounter = prolineCounter+1
            elif protein.singleAASequence[i] == 'S':
                serineCounter = serineCounter+1
            elif protein.singleAASequence[i] == 'T':
                threonineCounter = threonineCounter+1
            elif protein.singleAASequence[i] == 'W':
                tryptophanCounter = tryptophanCounter+1
            elif protein.singleAASequence[i] == 'Y':
                tyrosineCounter = tyrosineCounter+1
            elif protein.singleAASequence[i] == 'V':
                valineCounter = valineCounter+1
        
        if alanineCounter == 0:
            self.aaComposition.append(0)
        else:
            self.aaComposition.append(float(alanineCounter)/sequenceSize)
        if arginineCounter == 0:
            self.aaComposition.append(0)
        else:
            self.aaComposition.append(float(arginineCounter)/sequenceSize)
        if asparagineCounter == 0:
            self.aaComposition.append(0)
        else:
            self.aaComposition.append(float(asparagineCounter)/sequenceSize)
        if asparticAcidCounter == 0:
            self.aaComposition.append(0)
        else:
            self.aaComposition.append(float(asparticAcidCounter)/sequenceSize)
        if cysteineCounter == 0:
            self.aaComposition.append(0)
        else:
            self.aaComposition.append(float(cysteineCounter)/sequenceSize)
        if glutamicAcidCounter == 0:
            self.aaComposition.append(0)
        else:
            self.aaComposition.append(float(glutamicAcidCounter)/sequenceSize)
        if glutamineCounter == 0:
            self.aaComposition.append(0)
        else:
            self.aaComposition.append(float(glutamineCounter)/sequenceSize)
        if glycineCounter == 0:
            self.aaComposition.append(0)
        else:
            self.aaComposition.append(float(glycineCounter)/sequenceSize)
        if histidineCounter == 0:
            self.aaComposition.append(0)
        else:
            self.aaComposition.append(float(histidineCounter)/sequenceSize)
        if isoleucineCounter == 0:
            self.aaComposition.append(0)
        else:
            self.aaComposition.append(float(isoleucineCounter)/sequenceSize)
        if leucineCounter == 0:
            self.aaComposition.append(0)
        else:
            self.aaComposition.append(float(leucineCounter)/sequenceSize)
        if lysineCounter == 0:
            self.aaComposition.append(0)
        else:
            self.aaComposition.append(float(lysineCounter)/sequenceSize)
        if methioneCounter == 0:
            self.aaComposition.append(0)
        else:
            self.aaComposition.append(float(methioneCounter)/sequenceSize)
        if phenylalanineCounter == 0:
            self.aaComposition.append(0)
        else:
            self.aaComposition.append(float(phenylalanineCounter)/sequenceSize)
        if prolineCounter == 0:
            self.aaComposition.append(0)
        else:
            self.aaComposition.append(float(prolineCounter)/sequenceSize)
        if serineCounter == 0:
            self.aaComposition.append(0)
        else:
            self.aaComposition.append(float(serineCounter)/sequenceSize)
        if threonineCounter == 0:
            self.aaComposition.append(0)
        else:
            self.aaComposition.append(float(threonineCounter)/sequenceSize)
        if tryptophanCounter == 0:
            self.aaComposition.append(0)
        else:
            self.aaComposition.append(float(tryptophanCounter)/sequenceSize)
        if tyrosineCounter == 0:
            self.aaComposition.append(0)
        else:
            self.aaComposition.append(float(tyrosineCounter)/sequenceSize)
        if valineCounter == 0:
            self.aaComposition.append(0)
        else:
            self.aaComposition.append(float(valineCounter)/sequenceSize)   
        return self.aaComposition
    
    def calculateDipeptideComposition(self,protein):
        aminoAcidSingleLetters = protein.aaDictionary.values()
        aminoAcidSingleLetters.sort()
        proteinSequenceSize = len(protein.singleAASequence)
        dipeptideComposition =[]
        for i in range(len(aminoAcidSingleLetters)):
            for j in range(len(aminoAcidSingleLetters)):
                dipeptideCountValue = 0
                for k in range(len(protein.singleAASequence)):
                    if k < len(protein.singleAASequence)-1:
                        if protein.singleAASequence[k] == aminoAcidSingleLetters[i] and protein.singleAASequence[k+1] == aminoAcidSingleLetters[j]:
                           dipeptideCountValue = dipeptideCountValue+1
                if proteinSequenceSize == 0 or int(dipeptideCountValue) ==0:
                    dipeptideComposition.append(0)
                else:
                    dipeptideComposition.append(float(dipeptideCountValue)/proteinSequenceSize) 
        return dipeptideComposition

    def calculateMoreauBrotoAutoCorrelationDescriptor(self,protein):
        '''
        For each descriptor calculate the mean and standard deviation,
        then z-scale.
        '''
        scaledAverageFlexibilityDictionary = self.ZScale(self.averageFlexibilityDictionary)
        scaledFreeEnergyOfSolutionInWaterDictionary = self.ZScale(self.freeEnergyOfSolutionInWaterDictionary)
        scaledHydrophobicityDictionary = self.ZScale(self.hydrophobicityDictionary)
        scaledPolarizabilityDictionary = self.ZScale(self.polarizabilityDictionary)
        scaledRelativeMutabilityDictionary = self.ZScale(self.relativeMutabilityDictionary)
        scaledResidueASAInTripeptideDictionary = self.ZScale(self.residueASAInTripeptideDictionary)
        scaledResidueVolumeDictionary = self.ZScale(self.residueVolumeDictionary)
        scaledStericParameterDictionary = self.ZScale(self.stericParameterDictionary)
        #loop through sequence vector at i and i+d
        autocorrelationDescriptors = []
        def calculateIndividualVector(aDictionary):
            autocorrelationValue = []
            for i in range(1,31):
                descriptorValue = 0
                for k in range(len(protein.singleAASequence)):
                    if k < len(protein.singleAASequence)-i:
                        val = aDictionary[protein.singleAASequence[k]]* aDictionary[protein.singleAASequence[k+i]]
                        descriptorValue = descriptorValue+val
                autocorrelationValue.append(descriptorValue)
            return autocorrelationValue
        
        hydrophobicityAutoCorrelationDescriptors = calculateIndividualVector(scaledHydrophobicityDictionary)
        scaledAverageFlexibilityAutoCorrelationDescriptors = calculateIndividualVector(scaledAverageFlexibilityDictionary)
        scaledFreeEnergyOfSolutionInWaterAutoCorrelationDescriptors = calculateIndividualVector(scaledFreeEnergyOfSolutionInWaterDictionary)
        scaledPolarizabilityAutoCorrelationDescriptors = calculateIndividualVector(scaledPolarizabilityDictionary)
        scaledRelativeMutabilityAutoCorrelationDescriptors = calculateIndividualVector(scaledRelativeMutabilityDictionary)
        scaledResidueASAInTripeptideAutoCorrelationDescriptors = calculateIndividualVector(scaledResidueASAInTripeptideDictionary)
        scaledResidueVolumeAutoCorrelationDescriptors = calculateIndividualVector(scaledResidueVolumeDictionary)
        scaledStericParameterAutoCorrelationDescriptors = calculateIndividualVector(scaledStericParameterDictionary)
        autocorrelationDescriptors.append(hydrophobicityAutoCorrelationDescriptors)
        autocorrelationDescriptors.append(scaledAverageFlexibilityAutoCorrelationDescriptors)
        autocorrelationDescriptors.append(scaledFreeEnergyOfSolutionInWaterAutoCorrelationDescriptors)
        autocorrelationDescriptors.append(scaledPolarizabilityAutoCorrelationDescriptors)
        autocorrelationDescriptors.append(scaledRelativeMutabilityAutoCorrelationDescriptors)
        autocorrelationDescriptors.append(scaledResidueASAInTripeptideAutoCorrelationDescriptors)
        autocorrelationDescriptors.append(scaledResidueVolumeAutoCorrelationDescriptors)
        autocorrelationDescriptors.append(scaledStericParameterAutoCorrelationDescriptors)
        self.MoreauBrotoDescriptors = autocorrelationDescriptors
        return self.MoreauBrotoDescriptors

    def ZScale(self,dictionary):
        '''
        Z-scales dictionary values.
        '''
        meanValue = sum(dictionary.values())/20
        values = dictionary.values()
        summedVal = 0
        for i in range(len(values)):
            individualValue = pow(values[i]-meanValue,2)
            summedVal = summedVal+individualValue
        standardDeviation = sqrt(summedVal/20)
        keys = dictionary.keys()
        zscaledValues = []
        for i in range(len(keys)):
            zscaledValues.append((dictionary[keys[i]]-meanValue)/standardDeviation)
        standardizedDictionary = dict(zip(keys, zscaledValues))
        return standardizedDictionary
    
    def calculateMoranAutoCorrelation(self,protein):
        scaledAverageFlexibilityDictionary = self.ZScale(self.averageFlexibilityDictionary)
        scaledFreeEnergyOfSolutionInWaterDictionary = self.ZScale(self.freeEnergyOfSolutionInWaterDictionary)
        scaledHydrophobicityDictionary = self.ZScale(self.hydrophobicityDictionary)
        scaledPolarizabilityDictionary = self.ZScale(self.polarizabilityDictionary)
        scaledRelativeMutabilityDictionary = self.ZScale(self.relativeMutabilityDictionary)
        scaledResidueASAInTripeptideDictionary = self.ZScale(self.residueASAInTripeptideDictionary)
        scaledResidueVolumeDictionary = self.ZScale(self.residueVolumeDictionary)
        scaledStericParameterDictionary = self.ZScale(self.stericParameterDictionary)
        #Calculate average property along the sequence for each property value
        
        def calculateAverageProperty(dictionary):
            val = 0
            for i in range(len(protein.singleAASequence)):
                val = val+dictionary[protein.singleAASequence[i]]
            if val  == 0 or len(protein.singleAASequence) == 0:
                meanVal = 0
            else:
                meanVal = val/len(protein.singleAASequence)
            return meanVal
        meanScaledAverageFlexibilityValue = calculateAverageProperty(scaledAverageFlexibilityDictionary)
        meanScaledFreeEnergyOfSolutionInWaterValue = calculateAverageProperty(scaledFreeEnergyOfSolutionInWaterDictionary)
        meanScaledHydrophobicityValue = calculateAverageProperty(scaledHydrophobicityDictionary)
        meanScaledPolarizabilityValue = calculateAverageProperty(scaledPolarizabilityDictionary)
        meanScaledRelativeMutabilityValue = calculateAverageProperty(scaledRelativeMutabilityDictionary)
        meanScaledResidueASAInTripeptideValue = calculateAverageProperty(scaledResidueASAInTripeptideDictionary)
        meanScaledResidueVolumeValue = calculateAverageProperty(scaledResidueVolumeDictionary)
        meanScaledStericParameterValue = calculateAverageProperty(scaledStericParameterDictionary)
        autocorrelationDescriptors = []
        
        def calculateIndividualVector(aDictionary,aMeanValue):
            autocorrelationValue = []
            for i in range(1,31):
                descriptorValue = 0
                denominator = 0
                for k in range(len(protein.singleAASequence)):
                    if k < len(protein.singleAASequence)-i:
                        val = (aDictionary[protein.singleAASequence[k]]-aMeanValue)*(aDictionary[protein.singleAASequence[k+i]]-aMeanValue)
                        descriptorValue = descriptorValue+val
                        denomin = pow((aDictionary[protein.singleAASequence[k]]-aMeanValue),2)
                        denominator = denominator+denomin
                if len(protein.singleAASequence)==0 or descriptorValue == 0:
                    numerator = 1
                else:
                    numerator = descriptorValue/(len(protein.singleAASequence)-i)
                if denominator == 0 or len(protein.singleAASequence)==0:
                    denominator = 1
                else:
                    denominator = denominator/(len(protein.singleAASequence))
                autocorrelationValue.append(numerator/denominator)
            return autocorrelationValue
        
        autocorrelationDescriptors.append(calculateIndividualVector(scaledAverageFlexibilityDictionary,meanScaledAverageFlexibilityValue))
        autocorrelationDescriptors.append(calculateIndividualVector(scaledFreeEnergyOfSolutionInWaterDictionary,meanScaledFreeEnergyOfSolutionInWaterValue))
        autocorrelationDescriptors.append(calculateIndividualVector(scaledHydrophobicityDictionary,meanScaledHydrophobicityValue))
        autocorrelationDescriptors.append(calculateIndividualVector(scaledPolarizabilityDictionary,meanScaledPolarizabilityValue))
        autocorrelationDescriptors.append(calculateIndividualVector(scaledRelativeMutabilityDictionary,meanScaledRelativeMutabilityValue))
        autocorrelationDescriptors.append(calculateIndividualVector(scaledResidueASAInTripeptideDictionary,meanScaledResidueASAInTripeptideValue))
        autocorrelationDescriptors.append(calculateIndividualVector(scaledResidueVolumeDictionary,meanScaledResidueVolumeValue))
        autocorrelationDescriptors.append(calculateIndividualVector(scaledStericParameterDictionary,meanScaledStericParameterValue)) 
        self.MoranAutoCorrelation = autocorrelationDescriptors
        return self.MoranAutoCorrelation
    
    def calculateGearyAutoCorrelation(self,protein):
        scaledAverageFlexibilityDictionary = self.ZScale(self.averageFlexibilityDictionary)
        scaledFreeEnergyOfSolutionInWaterDictionary = self.ZScale(self.freeEnergyOfSolutionInWaterDictionary)
        scaledHydrophobicityDictionary = self.ZScale(self.hydrophobicityDictionary)
        scaledPolarizabilityDictionary = self.ZScale(self.polarizabilityDictionary)
        scaledRelativeMutabilityDictionary = self.ZScale(self.relativeMutabilityDictionary)
        scaledResidueASAInTripeptideDictionary = self.ZScale(self.residueASAInTripeptideDictionary)
        scaledResidueVolumeDictionary = self.ZScale(self.residueVolumeDictionary)
        scaledStericParameterDictionary = self.ZScale(self.stericParameterDictionary)
        #Calculate average property along the sequence for each property value
        def calculateAverageProperty(dictionary):
            val = 0
            for i in range(len(protein.singleAASequence)):
                val = val+dictionary[protein.singleAASequence[i]]
            if len(protein.singleAASequence)==0:
                meanVal=0
            else:
                meanVal = val/len(protein.singleAASequence)
            return meanVal
        meanScaledAverageFlexibilityValue = calculateAverageProperty(scaledAverageFlexibilityDictionary)
        meanScaledFreeEnergyOfSolutionInWaterValue = calculateAverageProperty(scaledFreeEnergyOfSolutionInWaterDictionary)
        meanScaledHydrophobicityValue = calculateAverageProperty(scaledHydrophobicityDictionary)
        meanScaledPolarizabilityValue = calculateAverageProperty(scaledPolarizabilityDictionary)
        meanScaledRelativeMutabilityValue = calculateAverageProperty(scaledRelativeMutabilityDictionary)
        meanScaledResidueASAInTripeptideValue = calculateAverageProperty(scaledResidueASAInTripeptideDictionary)
        meanScaledResidueVolumeValue = calculateAverageProperty(scaledResidueVolumeDictionary)
        meanScaledStericParameterValue = calculateAverageProperty(scaledStericParameterDictionary)
        autocorrelationDescriptors = []
        
        def calculateIndividualVector(aDictionary,aMeanValue):
            autocorrelationValue = []
            for i in range(1,31):
                descriptorValue = 0
                denominator = 0
                for k in range(len(protein.singleAASequence)):
                    if k < len(protein.singleAASequence)-i:
                        val = pow((aDictionary[protein.singleAASequence[k]]-aDictionary[protein.singleAASequence[k+i]]),2)
                        descriptorValue = descriptorValue+val
                        denomin = pow((aDictionary[protein.singleAASequence[k]]-aMeanValue),2)
                        denominator = denominator+denomin
                if descriptorValue==0 or len(protein.singleAASequence)==0:
                    numerator=1
                else:
                    numerator = descriptorValue/(2*(len(protein.singleAASequence)-i))
                if len(protein.singleAASequence)==0:
                    denominator=1
                else:
                    denominator = denominator/(len(protein.singleAASequence)-1)
                if numerator == 0 or denominator==0:
                    autocorrelationValue.append(1)
                else:
                    autocorrelationValue.append(numerator/denominator)
            return autocorrelationValue
    
        autocorrelationDescriptors.append(calculateIndividualVector(scaledAverageFlexibilityDictionary,meanScaledAverageFlexibilityValue))
        autocorrelationDescriptors.append(calculateIndividualVector(scaledFreeEnergyOfSolutionInWaterDictionary,meanScaledFreeEnergyOfSolutionInWaterValue))
        autocorrelationDescriptors.append(calculateIndividualVector(scaledHydrophobicityDictionary,meanScaledHydrophobicityValue))
        autocorrelationDescriptors.append(calculateIndividualVector(scaledPolarizabilityDictionary,meanScaledPolarizabilityValue))
        autocorrelationDescriptors.append(calculateIndividualVector(scaledRelativeMutabilityDictionary,meanScaledRelativeMutabilityValue))
        autocorrelationDescriptors.append(calculateIndividualVector(scaledResidueASAInTripeptideDictionary,meanScaledResidueASAInTripeptideValue))
        autocorrelationDescriptors.append(calculateIndividualVector(scaledResidueVolumeDictionary,meanScaledResidueVolumeValue))
        autocorrelationDescriptors.append(calculateIndividualVector(scaledStericParameterDictionary,meanScaledStericParameterValue)) 
        self.GearyAutoCorrelation = autocorrelationDescriptors
        return self.GearyAutoCorrelation
    