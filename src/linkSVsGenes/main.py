"""
	Set of scripts intended to link SVs to genes by rules based on gains and losses of regulatory elements caused by TAD disruptions

	The idea is that we first make a 'neighborhood' in which all regulatory elements are assigned to genes (regulator set).
	These are the regulatory elements within the TAD of the gene, and specifically what is linked to the gene in the data for eQTLs, promoters and enhancers.
	Then we look at the TADs that are disrupted by SVs, and determine how the interactions between genes and regulatory elements are re-wired.

	The setup of these scripts is:
	
	- The main script where the input (genes) are parsed and the relevant scripts are called
	- The neighborhoodDefiner, which takes the genes as input and assigns their regulator set. It also accepts SVs as input, and then makes a call to link them to the genes.
	- The derivativeTADMaker, which for each SV determines how the TAD is disrupted, and then assigns to each gene what the gains and losses in the regulator set are.
	- The geneRanking script, which then takes the genes and SVs that are linked to these, and outputs a list of SV-gene pairs and which regulatory elements are gained and lost for that gene.
	This used to provide a ranking, but not anymore. It also outputs the bags used for mulitple instance learning.

	The main script has 4 input parameters:
	- the uuid; this will be the output name of the folder that the data is written to.
	- permutationYN; if set to True, will randomly shuffle the SVs.
	- permutationRound; if running multiple permutations, this can be assigned a different value to get different output files for each run within the uuid folder.
	- path; path to where the settings (with all data locations) can be found.


"""
#############
###Imports###

import sys
import numpy as np
import random
import pickle as pkl
import os
import re
import time

path = sys.argv[4]
sys.path.insert(1, path) #path to the settings file

from neighborhoodDefiner import NeighborhoodDefiner
from geneRanking import GeneRanking
from inputParser import InputParser
from genomicShuffler import GenomicShuffler
import settings

###############
###Main code###

startTime = time.time() #Keep run time of the program

#0. Collect all the relevant parameters here for clarity
uuid = sys.argv[1] #It is the folder name that will be created within the output folder specified in the settings, and where the output will be written to.
permutationYN = sys.argv[2] #True or False depending on if we want to permute or not (keep in mind True should be a string)

#1. Read and parse the causal genes and the nonCausal genes. For now, keep them separate to test on causal/non-causal genes separately
causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.

#Combine the genes into one set.
causalGenes = np.concatenate((causalGenes, nonCausalGenes), axis=0)

#2. Read the SVs, depending on the data source. Currently a TCGA format is supported (see data preprocessing), PCAWG and HMF
variantData = []

if settings.general['source'] == 'PCAWG':
	print("Reading SV data PCAWG")
	svDir = settings.files['svDir']
	#get the cancer type from the settings
	cancerType = settings.general['cancerType']
	svData = InputParser().getSVsFromFile_pcawg(svDir, cancerType)
	print(svData) #check if the SVs are read in correctly
if settings.general['source'] == 'HMF':
	print("Reading SV data")
	svDir = settings.files['svDir']
	cancerType = settings.general['cancerType']
	svData = InputParser().getSVs_hmf(svDir, cancerType)
	print(svData) #check if the SVs are read in correctly
if settings.general['source'] == 'HMF_simple':
	print("Reading SV data")
	svDir = settings.files['svDir']
	cancerType = settings.general['cancerType']
	svData = InputParser().getSVs_hmf_simple(svDir)
	print(svData) #check if the SVs are read in correctly
if settings.general['source'] == 'Nunes':
	print("Reading SV data Nunes et al.")
	svDir = settings.files['svDir']
	#get the cancer type from the settings
	cancerType = settings.general['cancerType']
	svData = InputParser().getSVs_nunes(svDir, cancerType)
	print(svData) #check if the SVs are read in correctly

#Check the SV distribution
delCount = 0
dupCount = 0
invCount = 0
itxCount = 0
for sv in svData:
	
	if sv[8].svType == 'DEL':
		delCount += 1
	elif sv[8].svType == 'DUP':
		dupCount += 1
	elif sv[8].svType == 'INV':
		invCount += 1
	elif sv[8].svType == 'ITX':
		itxCount += 1
		
print(delCount, dupCount, invCount, itxCount)

#3. If this is a permutation run, we wish to shuffle these SVs
if permutationYN == "True":
	print("Shuffling variants")
	genomicShuffler = GenomicShuffler()
	#Shuffle the variants
	svData = genomicShuffler.shuffleSVs(svData)
	
permutationRound = ""
if permutationYN == "True":
	permutationRound = sys.argv[3] #will be appended to the output file name. 

#2. Get the neighborhood for these genes based on the SVs
print("Defining the neighborhood for the genes and the SVs")
NeighborhoodDefiner(causalGenes, svData)

#3. Get the SV-gene pairs and their gains and losses, output this to a file, and create the bags for MIL.
print("Gathering gains and losses for SV-gene pairs and output to a file")
geneRanking = GeneRanking(causalGenes[:,3], svData, sys.argv[1], permutationRound)

endTime = time.time()
print("The program took ", endTime-startTime, " seconds to complete")

	
