'''
Created on Oct 24, 2010
Utility class for atoms and proteins.
@author: ed
'''
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
