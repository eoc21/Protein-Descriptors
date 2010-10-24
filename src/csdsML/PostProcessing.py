'''
Created on May 15, 2010
Simple script to remove file ids so descriptors can be put into weka
@author: ed
'''
import os, sys

datafile = open(sys.argv[1],'r')
formattedOut = open(sys.argv[2],'w')

data = datafile.readlines()
for i in range(len(data)):
    s = str(data[i]).split(",")
   # s1 = str(s[1:])
   # s1 =s1.replace('[','')
   # s2 = s1.replace(']','')
   # s3 = s2.replace('\'','')
    formattedOut.write(str(s[1:]))
