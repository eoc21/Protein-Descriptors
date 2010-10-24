'''
Created on May 11, 2010
Descriptor calculates Ultrafast shape descriptors for proteins.
@author: ed
'''
import os, sys, math, operator
from PDBReader import *

class Point():
    def __init__(self,x=0,y=0,z=0):
        self.xCoordinate = x
        self.yCoordinate = y
        self.zCoordinate = z

class AtomIndexAndDistance():
    def __init__(self,index,distance):
        self.index = index
        self.distance = distance
    
class USR():
    def __init__(self):
        pass
    def CalculateCentroid(self,aProtein):
        centroid = Point()
        summedX = 0.0
        summedY = 0.0
        summedZ = 0.0
        atomsInProtein = len(aProtein.getAtoms())
        for i in range(len(aProtein.getAtoms())):
            summedX = summedX+float(aProtein.getAtom(i).getXCoordinate())
            summedY = summedY+float(aProtein.getAtom(i).getYCoordinate())
            summedZ = summedZ+float(aProtein.getAtom(i).getZCoordinate()) 
           
        if atomsInProtein == 0 or summedX == 0:
            centroid.xCoordinate = 0
        if atomsInProtein == 0 or  summedY == 0:
            centroid.yCoordinate = 0
        if atomsInProtein == 0 or summedZ == 0:
            centroid.zCoordinate = 0
        else:
            centroid.xCoordinate = summedX/atomsInProtein
            centroid.yCoordinate = summedY/atomsInProtein
            centroid.zCoordinate = summedZ/atomsInProtein
        return centroid

    def ClosestAndFurthestAtomToCentroid(self,aProtein):
        centroid = self.CalculateCentroid(aProtein)
        atomIndexWithDistance=[]
        closestAndFurthestAtom=[]
        for i in range(len(aProtein.getAtoms())):
            dist = self.EuclideanDistanceMeasure(aProtein,centroid,i)
            result = AtomIndexAndDistance(i,dist)
            atomIndexWithDistance.append(result)
        atomIndexWithDistance.sort(key = operator.attrgetter('distance'))
        closestAndFurthestAtom.append(atomIndexWithDistance[0])
        closestAndFurthestAtom.append(atomIndexWithDistance[len(atomIndexWithDistance)-1])
        return closestAndFurthestAtom

    def FurthestAtomFromFurthestAtom(self,aProtein,furthestAtom):
        atomIndexWithDistance=[]
        furthestAtomIndex = furthestAtom.index
        x = float(aProtein.getAtom(furthestAtomIndex).getXCoordinate())
        y = float(aProtein.getAtom(furthestAtomIndex).getYCoordinate())
        z = float(aProtein.getAtom(furthestAtomIndex).getZCoordinate())
    
        def Distance():
            for i in range(len(aProtein.getAtoms())):
                xvalue = math.pow(float(aProtein.getAtom(i).getXCoordinate()) -x,2)
                yvalue = math.pow(float(aProtein.getAtom(i).getYCoordinate())-y,2)
                zvalue = math.pow(float(aProtein.getAtom(i).getZCoordinate())-z,2)
                distance = math.sqrt(xvalue+yvalue+zvalue)
                result = AtomIndexAndDistance(i,distance)
                atomIndexWithDistance.append(result)
            atomIndexWithDistance.sort(key = operator.attrgetter('distance'))
            return atomIndexWithDistance
        atomsWithDistance = Distance()
        furthestFromFurthestAtomIndex = atomsWithDistance[len(atomsWithDistance)-1].index
        return furthestFromFurthestAtomIndex

    def EuclideanDistanceMeasure(self,protein,b,i):
        xvalue = math.pow((float(protein.getAtom(i).getXCoordinate()))-float(b.xCoordinate),2)
        yvalue = math.pow((float(protein.getAtom(i).getYCoordinate()))-float(b.yCoordinate),2)
        zvalue = math.pow((float(protein.getAtom(i).getZCoordinate()))-float(b.zCoordinate),2)
        distance = math.sqrt((xvalue+yvalue+zvalue))
        return distance

    def MomentToCentroid(self,momentId,aProtein,centroidPoint):
    
        def MeanFromCentroid():
            distance = 0
            for i in range(len(aProtein.getAtoms())):
                       xval = math.pow(math.fabs(float(aProtein.getAtom(i).getXCoordinate())-centroidPoint.xCoordinate),2)
                       yval = math.pow(math.fabs(float(aProtein.getAtom(i).getYCoordinate())-centroidPoint.yCoordinate),2)
                       zval = math.pow(math.fabs(float(aProtein.getAtom(i).getZCoordinate())-centroidPoint.zCoordinate),2)
                       dist =  math.sqrt(xval+yval+zval)
                       distance = distance+dist
            mean = distance/len(aProtein.getAtoms())
            return mean
        if momentId == 1:
            meanValue = MeanFromCentroid()
            return meanValue
    
        if momentId == 2:
            summedDistance = 0
            meanDistanceToCentroid = MeanFromCentroid()
            for i in  range(len(aProtein.getAtoms())):
                xval = math.pow(math.fabs(float(aProtein.getAtom(i).getXCoordinate())-centroidPoint.xCoordinate),2)
                yval = math.pow(math.fabs(float(aProtein.getAtom(i).getYCoordinate())-centroidPoint.yCoordinate),2)
                zval = math.pow(math.fabs(float(aProtein.getAtom(i).getZCoordinate())-centroidPoint.zCoordinate),2)
                dist =  math.sqrt(xval+yval+zval)
                summedDistance = summedDistance+math.pow((dist-meanDistanceToCentroid),2)
            variance = summedDistance/len(aProtein.getAtoms())
            return variance
    
        if momentId == 3:
            summedDistance = 0
            numeratorSum =  0
            denominatorSum =  0
            meanDistanceToCentroid = MeanFromCentroid()
            for i in  range(len(aProtein.getAtoms())):
                xval = math.pow(math.fabs(float(aProtein.getAtom(i).getXCoordinate())-centroidPoint.xCoordinate),2)
                yval = math.pow(math.fabs(float(aProtein.getAtom(i).getYCoordinate())-centroidPoint.yCoordinate),2)
                zval = math.pow(math.fabs(float(aProtein.getAtom(i).getZCoordinate())-centroidPoint.zCoordinate),2)
                dist =  math.sqrt(xval+yval+zval)
                numeratorSum = numeratorSum+math.pow(math.fabs((dist-meanDistanceToCentroid)),3)
                denominatorSum = denominatorSum+math.pow(math.fabs((dist-meanDistanceToCentroid)),2)
            numerator = numeratorSum/len(aProtein.getAtoms())
            denominator = math.pow(denominatorSum/len(aProtein.getAtoms()),1.5)
            skewness = numerator/denominator
            return skewness

    def MomentToX(self,MomentId,aProtein,atomIndexValue):
        def MeanFromPoint():
            distance = 0
            for i in range(len(aProtein.getAtoms())):
                xval = math.pow(math.fabs(float(aProtein.getAtom(i).getXCoordinate())-float(aProtein.getAtom(atomIndexValue).getXCoordinate())),2)
                yval = math.pow(math.fabs(float(aProtein.getAtom(i).getYCoordinate())-float(aProtein.getAtom(atomIndexValue).getYCoordinate())),2)
                zval = math.pow(math.fabs(float(aProtein.getAtom(i).getZCoordinate())-float(aProtein.getAtom(atomIndexValue).getZCoordinate())),2)
                dist =  math.sqrt(xval+yval+zval)
                distance = distance+dist
            mean = distance/len(aProtein.getAtoms())
            return mean
        if MomentId == 1:
            meanValue = MeanFromPoint()
            return meanValue
    
        if MomentId == 2:
            summedDistance = 0
            meanDistance = MeanFromPoint()
            for i in  range(len(aProtein.getAtoms())):
                xval = math.pow(math.fabs(float(aProtein.getAtom(i).getXCoordinate())-float(aProtein.getAtom(atomIndexValue).getXCoordinate())),2)
                yval = math.pow(math.fabs(float(aProtein.getAtom(i).getYCoordinate())-float(aProtein.getAtom(atomIndexValue).getYCoordinate())),2)
                zval = math.pow(math.fabs(float(aProtein.getAtom(i).getZCoordinate())-float(aProtein.getAtom(atomIndexValue).getZCoordinate())),2)
                dist =  math.sqrt(xval+yval+zval)
                summedDistance = summedDistance+math.pow((dist-meanDistance),2)
            variance = summedDistance/len(aProtein.getAtoms())
            return variance

        if MomentId == 3:
            summedDistance = 0
            numeratorSum =  0
            denominatorSum =  0
            meanDistance = MeanFromPoint()
            for i in range(len(aProtein.getAtoms())):
                xval = math.pow(math.fabs(float(aProtein.getAtom(i).getXCoordinate())-float(aProtein.getAtom(atomIndexValue).getXCoordinate())),2)
                yval = math.pow(math.fabs(float(aProtein.getAtom(i).getYCoordinate())-float(aProtein.getAtom(atomIndexValue).getYCoordinate())),2)
                zval = math.pow(math.fabs(float(aProtein.getAtom(i).getZCoordinate())-float(aProtein.getAtom(atomIndexValue).getZCoordinate())),2)
                dist =  math.sqrt(xval+yval+zval)
                numeratorSum = numeratorSum+math.pow(math.fabs((dist-meanDistance)),3)
                denominatorSum = denominatorSum+math.pow(math.fabs((dist-meanDistance)),2)
            numerator = numeratorSum/len(aProtein.getAtoms())
            denominator = math.pow(denominatorSum/len(aProtein.getAtoms()),1.5)
            skewness = numerator/denominator
            return skewness