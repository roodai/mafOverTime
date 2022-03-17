# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 13:34:31 2022

@author: inf-29-2021
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
    #Amount of individuals needed outside of function for standard error
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
            #If allele is heterozygotic add one to the count of major
            if k[0]!=k[1]:
               major += snpdict[k]
        #Calculate the allele freqeuncy of the most abudant allele if the sum of values isnt zero
        if sum(snpdict.values()) !=0:
            #The major allele in homozygotic alleles adds two to the count of major
            major+=max(snpdict.values())*2
            #divide the count of major by the 2*sum of all alleles
            allelefreq=major/(sum(snpdict.values())*2)
            #calculate maf, which is 1 substracted the major allele frequency
            maf=1-allelefreq
        #if there is only one type of allele there is not minor allele
        else:
            maf=0
        #add to the global variable dict having snpids as keys and MAFs as values
        #set.default adds the snp id as key if it isnt already there, if it is the corresponding MAF will be appended to the list
        snpmaf.setdefault(snpid,[]).append(maf)
        
#Set initial time value to 0
time=[0]
for index, arg in enumerate(sys.argv):
    if index%2 ==0 and index+2<len(sys.argv) and index!=0:
        #Get values for time between datasets
        try:
            time.append(float(input('Time between datasets ' + sys.argv[index]+ ' and ' + sys.argv[index+2] + ': ' )))
        except ValueError:
            print('Input needs to be numerical')
    if index%2 !=0:
        #Run mafCalc for inputs
        with open(sys.argv[index], 'r') as dotMap, open(sys.argv[index+1], 'r') as ped:
            mafCalc(ped, dotMap)

#Time provided is between datasets, cumulative sum needed to depict proper time
time=list(np.cumsum(time))

#Get snps, or other IDs of interest
try:
    snpOfInterest=str(input('Provide SNP(s) of interest: '))
    snpOfInterest=snpOfInterest.split()
    
    #Iterate through snps of interest for plotting purposes
    for snip in snpOfInterest:
        #plot time versus maf of snp interest with error bars
        #standard deviation of allele freqeucies used to determine error bars
        plt.errorbar(time, snpmaf[snip], marker='.',yerr=np.sqrt((snpmaf[snip]*(1-np.array(snpmaf[snip])))/(2*samplesize)))
    #Add descriptions to plot and close it
    plt.legend(snpOfInterest)
    plt.savefig('MAF over time')
    plt.close()
except KeyError:
    print('SNP ' + snip + ' not found in the datasets')