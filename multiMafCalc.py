# -*- coding: utf-8 -*-

"""
Title: multiMafCalc.py
Created on Mon Mar  7 13:34:31 2022
Author: Roni Odai
Description:
    This script takes .ped files as input, cross referces them to input .map files for SNP IDs.
    Minor allele frequency (MAF) is calculated for all SNPs and stored in a dictionary.
    MAF is also calculated for successive .ped file inputs.
    Time conditions between .ped files and SNPs of interest are procured as user input.
    MAF vs. time is plotted for all SNPs of interest.
List of functions:
    mafCalc(ped, dotMap)
    Parses .ped file, chunks nucleotides into alleles, cross referces corresponding .map file for IDs.
    Calculates MAF and appends into dicitonary with its IDs as keys.
List of non-standard modules
    Pyplot from Matpplotlib: plots Time Vs. MAF with error bars.
Outline:
    - Parse .ped files and compute allele counts
    - Cross reference to .map files and extract SNP IDs
    - Calculate allele frequencies
    - Get time and SNP IDs of interest from user
    - Plot MAF for SNP IDs of interest
Usage:
    python multiMafCalc.py test1.map test1.ped test2.map test2.ped test3.map test3.ped test4.map test4.ped
"""


import sys
import numpy as np
from matplotlib import pyplot as plt


snpmaf={}
#Define function with ped and map files as input
def mafCalc(ped, dotMap):
    indalleles=[]
    for line in ped:
        #Empty allele list once a line/individual has been iterated through
        alleles=[]
        #Extract genotypes from ped file
        snps=line.split()[6:len(line.split())]
        #Chunk genotypes into entries of two to make an allele
        for snp in range(0,len(snps), 2):
            allele=''.join(snps[snp:snp+2])
            alleles.append(allele)
        #Append all alleles in order of occurence to retain correspondence to an individual
        indalleles.append(alleles)
    global samplesize
    samplesize=len(indalleles)
    #Iterate through map file lines with an iterator
    for i, bline in enumerate(dotMap):
        #Extract snpid
        snpid=bline.split()[1]
        #Empty dictionary for snps once it has been iterated
        snpdict=dict()
        #Iterate through indalleles listed list, list is the individual, listed list the individuals alleles
        for j in range(0,len(indalleles)):
            #Add alleles and increment their values once per occurrence
            if indalleles[j][i] not in snpdict:
                snpdict[indalleles[j][i]]=0
            if indalleles[j][i] in snpdict:
                snpdict[indalleles[j][i]]+=1
        major=0
        #Set the value of keys/alleles containing '0' to nothing
        for k in snpdict:
            if '0' in k:
                snpdict[k]=0
            if k[0]!=k[1]:
               major += snpdict[k]
        #Calculate the allele freqeuncy of the most abudant allele if the sum of values isnt zero
        if sum(snpdict.values()) !=0:
            major+=max(snpdict.values())*2
            allelefreq=major/(sum(snpdict.values())*2)
            #calculate maf, which is 1 substracted the major allele frequency
            maf=1-allelefreq
        else:
            maf=0
        #add to the global variable dict having snpids as keys and MAFs as values
        #set.default adds the snp id as key if it isnt already there, if it is the corresponding MAF will be appended to the list
        snpmaf.setdefault(snpid,[]).append(maf)

time=[0]
for index, arg in enumerate(sys.argv):
    if index%2 ==0 and index+2<len(sys.argv) and index!=0:
        #Get values for time between datasets
        try:
            time.append(float(input('Time between datasets ' + sys.argv[index]+ ' and ' + sys.argv[index+2] + ': ' )))
        except ValueError:
            print('Input needs to be numerical')
    if index%2 !=0:
        #Run mafCalc for
        with open(sys.argv[index], 'r') as dotMap, open(sys.argv[index+1], 'r') as ped:
            mafCalc(ped, dotMap)

time=list(np.cumsum(time))

try:
    snpOfInterest=str(input('Provide SNP(s) of interest: '))


    snpOfInterest=snpOfInterest.split()

    for snip in snpOfInterest:
        plt.errorbar(time, snpmaf[snip], marker='.',yerr=np.sqrt((snpmaf[snip]*(1-np.array(snpmaf[snip])))/(2*samplesize)))
    plt.legend(snpOfInterest)
    plt.xlabel('Time')
    plt.ylabel('MAF')
    plt.savefig('MAF over time.pdf')
    plt.close()
except KeyError:
    print('SNP ' + snip + ' not found in the datasets')
