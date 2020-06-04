#!/usr/bin/env python3

#####################################################################################################
# File: compareGenomes.py 
# Executable: python compareGenomes.py 
# Pupose: Program calculates the composition statistics of two genomes: GC content, codon usage, etc.
#
# Student: Carlos Arevalo(caeareva)
# Group: None 
#
#####################################################################################################
'''
This program calculates and compares two genomes gc content, codon usage and codon composition. 
The method used to calculate the codon usage and codon composition was by dividing the specific 
values in each function of both genomes and taking the log base 2. I think the program needs improvement
and I will appreacite your feedback how to develop it more efficiently. 
'''

import sequenceAnalysis as sa
import math

def main(): # comparedGenomes is developed inside the main() 
    '''
    Funtion reads over the input fasta file and executes functions. Function calculates the 
    sequence length, GC content, and amino acids compositions.
    '''
    # method reads first genome using FastAreader
    myReaderA =  sa.FastAreader('testGenome.fa') # make sure to change this to use stdin
    nucParamA = sa.NucParams()
    for header, sequence in myReaderA.readFasta():
        dnaCountA = nucParamA.nucCount()
        nucParamA.addSequence(sequence)
        pass

    # method reads second genome using FastAreader
    myReaderB =  sa.FastAreader('haloVolc1_1-genes.fa') # make sure to change this to use stdin
    nucParamB = sa.NucParams()
    for header, sequence in myReaderB.readFasta():
        nucParamB.addSequence(sequence)
        dnaCountB = nucParamB.nucCount()
        pass
        
    # calculates sequence lengths of genomes A and B and returns them in Mb multpiplying by (1/1000000)
    print("Sequence length A, B = {:0.2f} Mb, {:0.2f} Mb\n".format(dnaCountA*(1/10**6), (dnaCountB)*(1/10**6)))

    # access nucComposition() in nucParms class
    codonCompA = nucParamA.nucComposition()
    codonCompB = nucParamB.nucComposition()
    # access nucCount() in nucParms class
    dnaCountA = nucParamA.nucCount()
    dnaCountB = nucParamB.nucCount()
    # calculates the number of Gs and Cs in sequence for both genomes
    gcContentA = codonCompA['G'] + codonCompA['C']
    gcContentB = codonCompB['G'] + codonCompB['C']
    # calculates GC content for genome A
    gcContentA = (gcContentA/dnaCountA)*100
    # calculates GC content for genome B
    gcContentB = (gcContentB/dnaCountB)*100
    # compares both genomes GC contents
    bothGC = (math.log(gcContentA/gcContentB))
    pass

    # prints the GC content of genome A and genome B in order
    print("GC content A, B = {:0.2f}%, {:0.2f}%\n".format(gcContentA, gcContentB))
    # prints both genomes GC content comparison 
    print("GC content ratio of both genomes = {:0.2f}%\n".format(bothGC)) 
    
    # acccess codonComposition() function in nucParams class  
    codonCountA = nucParamA.codonComposition()
    codonCountB = nucParamB.codonComposition()
    # iterate over codon counts and amino acids in dictionary for both genomes
    for codons in sorted(codonCountA) or sorted(codonCountB):
        # access and sorts items in AA dictionary genome A
        charA = nucParamA.rnaCodonTable[codons]
        # access and sorts items in AA dictionary genome B
        charB = nucParamB.rnaCodonTable[codons]
        pass
        
        # calculates relative codon usage for each codon in both genomes
        compositionA = codonCountA[codons]
        compositionB = codonCountB[codons]
        # iterates over amino acids in both genomes
        for letter in charA or charB: 
            aminoCountA = nucParamA.aaComp[letter]
            aminoCountB = nucParamB.aaComp[letter] 
            pass
        
            if compositionA or compositionB != 0: # if composition is not zero, calculates the codon frequency 
                # calculates the relative composition of both genomes using the log function base 2
                relativeComp = math.log(((compositionA / aminoCountA) / (compositionB / aminoCountB)), 2)
                # calculates the amino acid composition of both genomes using the log function base 2
                aaCount = math.log(compositionA / compositionB, 2)
                print('{} : {} {:5.1f} ({:6.2f})'.format(codons, letter, relativeComp, aaCount))
        
        
if __name__ == "__main__":
    main()

