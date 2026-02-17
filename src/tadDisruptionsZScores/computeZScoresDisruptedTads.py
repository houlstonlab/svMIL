"""
	
	For each TAD, determine if these are disrupted by ANY SV.
	If yes, get all the genes in those TADs, and compute the z-score for those genes from the patient with the SV to all patients without any SV disrupting the boundary
	Then filter out the genes in both this set and the null distribution that have SNVs, CNVs or coding SVs, because that effect is likely not from the SV.

	We also report p-values, but those are not used elsewhere anymore.

"""
import sys

path = sys.argv[1]
sys.path.insert(1, path)

#this code depends on the input parser from the linking part. This is quick and dirty, is there a better solution?
sys.path.insert(1, 'linkSVsGenes/')

import settings
from inputParser import InputParser
import numpy as np
import os

import glob
import re
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy import stats
import matplotlib.pyplot as plt
from scipy.stats import rankdata
from genomicShuffler import GenomicShuffler

###parameters
geneNameConversionFile = settings.files['geneNameConversionFile']
expressionFile = settings.files['normalizedExpressionFile']
outDir = sys.argv[2]
randomize = sys.argv[3] #shuffle expression to get random z-scores?

specificOutDir = outDir + '/tadDisruptionsZScores/'

if not os.path.exists(specificOutDir):
	os.makedirs(specificOutDir)

#For each TAD, determine which genes are there

#first get all genes and their positions
causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.

#Combine the genes into one set.
allGenes = np.concatenate((causalGenes, nonCausalGenes), axis=0)

genes = []
for gene in allGenes:
	genes.append([gene[0], gene[1], gene[2], gene[3].name])

genes = np.array(genes, dtype='object')

#also use a map for the gene names
geneNameConversionMap = dict()
with open(geneNameConversionFile, 'r') as inF:

	lineCount = 0
	for line in inF:

		if lineCount < 1:
			lineCount += 1
			continue
		line = line.strip()
		splitLine = line.split("\t")
		ensgId = splitLine[3]
		splitEnsgId = ensgId.split('.') #we only keep everything before the dot
		geneName = splitLine[4]
		geneNameConversionMap[splitEnsgId[0]] = geneName

#Get all SVs
if settings.general['source'] == 'HMF':
	svDir = settings.files['svDir']
	svData = InputParser().getSVs_hmf(svDir, settings.general['cancerType'])
elif settings.general['source'] == 'PCAWG':
	svDir = settings.files['svDir']
	svData = InputParser().getSVsFromFile_pcawg(svDir, settings.general['cancerType'])
elif settings.general['source'] == 'Nunes':
	svDir = settings.files['svDir']
	svData = InputParser().getSVs_nunes(svDir, settings.general['cancerType'])
else:
	print('Other data sources not supported')
	exit(1)

#fix this
filteredSVs = svData

print(len(np.unique(filteredSVs[:,7])))

#load the expression data from the pre-normalized file.


#get the gene expression
expressionData = []
samples = []
with open(expressionFile, 'r') as inF:
	lineCount = 0
	for line in inF:
		line = line.strip()
		if lineCount == 0:

			samples = ['']
			samples += line.split("\t")

			lineCount += 1
			continue
		splitLine = line.split("\t")
		geneName = splitLine[0]

		#no need to check for genes that we did not check SVs for
		if geneName not in genes[:,3]:
			continue

		data = splitLine[1:len(splitLine)]

		fixedData = [geneName]
		fixedData += data

		expressionData.append(fixedData)

expressionData = np.array(expressionData, dtype="object")

#remove the columns that do not have SVs, and are thus not for this cancer type. 
columnsToRemove = []
filteredSamples = []
patientsWithSVs = np.unique(filteredSVs[:,7])
for sampleCol in range(1, len(samples)):
	if samples[sampleCol] not in patientsWithSVs:
		columnsToRemove.append(sampleCol)
	else:
		filteredSamples.append(samples[sampleCol])

expressionData = np.delete(expressionData, columnsToRemove, axis=1)
samples = [''] + filteredSamples


#Get the TADs.
#For each TAD, if there is an SV disrupting either of the boundaries, it goes into the disrupted class for that patient.
#All other patients that do not have SVs disrupting this TAD go into the negative group.
tadFile = settings.files['tadFile']

print("Getting TADs")
tadData = InputParser().getTADsFromFile(tadFile)

tadDisruptions = dict() #keep each TAD, and a list of which patients have an SV disrupting that TAD.
for tad in tadData:

	tadStr = tad[0] + '_' + str(tad[1]) + '_' + str(tad[2])
	tadDisruptions[tadStr] = []

	#Check if there is any SV overlapping the TAD on either side.

	#for intra-chromosomal SVs, we check if these overlap the boundary of the TAD.
	#for inter-chromosomal SVs, we check if these are inside a TAD.

	#get a subset of SVs on the right chromosome
	svChr1Subset = filteredSVs[filteredSVs[:,0] == tad[0]]
	intraSVSubset = svChr1Subset[svChr1Subset[:,0] == svChr1Subset[:,3]]

	#In the intra set, check if there are SVs that start before the TAD, but end within the TAD.

	#This TAD is disrupted if an SV starts or ends in it.
	#so an SV either starts before the end, and ends after the start
	#or the sv ends after the start, but starts before the end. 

	#SVs that start before the TAD end, and end after the TAD start. 
	startMatches = (intraSVSubset[:,1] >= tad[1]) * (intraSVSubset[:,1] <= tad[2]) * (intraSVSubset[:,5] >= tad[2])
	endMatches = (intraSVSubset[:,5] >= tad[1]) * (intraSVSubset[:,1] <= tad[1]) * (intraSVSubset[:,5] <= tad[2])
	
	#startMatches = intraSVSubset[:,1] <= tad[2]
	#endMatches = intraSVSubset[:,5] >= tad[1]
			
	#either of these must be true for the TAD to be disrupted by an intra-chromosomal SV.
	allMatches = intraSVSubset[startMatches + endMatches]
	
	for match in allMatches:
		
		svStr = match[0] + '_' + str(match[1]) + '_' + str(match[2]) + '_' + match[3] + '_' + str(match[4]) + '_' + str(match[5]) + '_' + match[7]

		if match[7] not in tadDisruptions[tadStr]:
			
			tadDisruptions[tadStr].append([match[7], match])

	#then repeat for interchromosomal SVs.
	#check if there is any bp in this TAD to count it as disrupted.
	
	#get a subset of all SVs that are interchromosomal, then subset for the ones that are either on chr1 or chr2
	
	interSVSubset = filteredSVs[filteredSVs[:,0] != filteredSVs[:,3]]
	interSVsChr1 = interSVSubset[interSVSubset[:,0] == tad[0]]
	interSVsChr2 = interSVSubset[interSVSubset[:,3] == tad[0]]
	
	#which bps of the chr1 set are in the TAD
	chr1Matches = (interSVsChr1[:,1] >= tad[1]) * (interSVsChr1[:,1] <= tad[2])
	chr2Matches = (interSVsChr2[:,5] >= tad[1]) * (interSVsChr2[:,5] <= tad[2])
	
	allChr1Matches = interSVsChr1[chr1Matches]
	allChr2Matches = interSVsChr2[chr2Matches]
	
	
	for match in allChr1Matches:
		
		svStr = match[0] + '_' + str(match[1]) + '_' + str(match[2]) + '_' + match[3] + '_' + str(match[4]) + '_' + str(match[5]) + '_' + match[7]
		
		if match[7] not in tadDisruptions[tadStr]:

			tadDisruptions[tadStr].append([match[7], match])
			
	for match in allChr2Matches:
		
		svStr = match[0] + '_' + str(match[1]) + '_' + str(match[2]) + '_' + match[3] + '_' + str(match[4]) + '_' + str(match[5]) + '_' + match[7]
		
		if match[7] not in tadDisruptions[tadStr]:

			tadDisruptions[tadStr].append([match[7], match])
	
#check which TADs we have with disruptions, validate

disrCount = 0
nonDisrCount = 0
for tad in tadData:
	
	tadStr = tad[0] + '_' + str(tad[1]) + '_' + str(tad[2])
	
	if len(tadDisruptions[tadStr]) > 0:
		disrCount += 1
	else:
		nonDisrCount += 1
		
print('disrupted tads: ', disrCount)
print('non-disrupted tads: ', nonDisrCount)

allPatientsWithDisruptions = dict()
for tad in tadDisruptions:
	for match in tadDisruptions[tad]:
		allPatientsWithDisruptions[match[0]] = 0

#check how many unique SVs disrupt a TAD pair
disruptingSVs = dict()
typeDistribution = dict()
for tad in tadDisruptions:
	
	for matchList in tadDisruptions[tad]:
		
		match = matchList[1]
		svStr = match[0] + '_' + str(match[1]) + '_' + str(match[2]) + '_' + match[3] + '_' + str(match[4]) + '_' + str(match[5]) + '_' + match[7]
		svType = match[8].svType
		
		
		if svType not in typeDistribution:
			typeDistribution[svType] = 0
		
		if svStr not in disruptingSVs: 
			typeDistribution[svType] += 1

			disruptingSVs[svStr] = 0


print('disrupting SVs: ', len(disruptingSVs))
print('per type: ', typeDistribution)

#Shuffle expression if requested
if randomize == 'True':

	genes = expressionData[:,0]
	expression = expressionData[:,1:]
	np.random.shuffle(expression)

	shuffledExpressionData = np.empty(expressionData.shape, dtype='object')
	shuffledExpressionData[:,0] = genes
	shuffledExpressionData[:,1:] = expression

	expressionData = shuffledExpressionData

#get all the mutation information for each patient, so that we can easily filter out mutated genes
mutDir = outDir + '/patientGeneMutationPairs/'
snvPatients = np.load(mutDir + 'snvPatients.npy', allow_pickle=True, encoding='latin1').item()

svPatientsDel = np.load(mutDir + 'svPatientsDel.npy', allow_pickle=True, encoding='latin1').item()
svPatientsDup = np.load(mutDir + 'svPatientsDup.npy', allow_pickle=True, encoding='latin1').item()
svPatientsInv = np.load(mutDir + 'svPatientsInv.npy', allow_pickle=True, encoding='latin1').item()
svPatientsItx = np.load(mutDir + 'svPatientsItx.npy', allow_pickle=True, encoding='latin1').item()

cnvPatientsDel = np.load(mutDir + 'cnvPatientsDel.npy', allow_pickle=True, encoding='latin1').item()
cnvPatientsAmp = np.load(mutDir + 'cnvPatientsAmp.npy', allow_pickle=True, encoding='latin1').item()


disruptedTadExpression = []
nonDisruptedTadExpression = []

disruptedPairs = dict() #store for each patient a dict with genes, and in there, the expression of the gene. 
nonDisruptedPairs = dict()


#in case of tcga data, update the sample names that do not match the SVs.
if settings.general['source'] == 'TCGA':
	fixedSamples = []
	for sample in samples:

		if sample == '':
			fixedSamples.append('')
			continue
		if sample == 'Hybridization REF':
			continue

		splitPatientID = sample.split("-")
		shortPatientID = settings.general['cancerType'] + splitPatientID[2]
		
		if shortPatientID not in allPatientsWithDisruptions:
			continue #some patients never have svs, so no need to look at those. 
		
		fixedSamples.append(shortPatientID)
	samples = fixedSamples

#tcga and pcawg data is an identifier hell, so map things back here that apparently have no mutations...
#maybe I missed them earlier on but it is hard to tell when nothing maps to each other...
if settings.general['source'] == 'TCGA' or settings.general['source'] == 'PCAWG':
	for sample in samples:
		if sample not in snvPatients:
			snvPatients[sample] = []
		if sample not in svPatientsDel:
			svPatientsDel[sample] = []
		if sample not in svPatientsDup:
			svPatientsDup[sample] = []
		if sample not  in svPatientsInv:
			svPatientsInv[sample] = []
		if sample not in svPatientsItx:
			svPatientsItx[sample] = []
		if sample not in cnvPatientsAmp:
			cnvPatientsAmp[sample] = []
		if sample not in cnvPatientsDel:
			cnvPatientsDel[sample] = []

for patient in samples:
	
	if patient == '' or patient == 'Hyridization REF' or patient == 'feature':
		continue
	
	if patient not in disruptedPairs:
		disruptedPairs[patient] = dict()

for gene in allGenes:
	if gene[3].name not in nonDisruptedPairs:
		nonDisruptedPairs[gene[3].name] = dict()

#for each gene, remove the genes that overlap multiple TADs. we dno't know which TAD these belong to.
filteredGenes = []
for gene in allGenes:
	
	#check which tads is overlaps with
	tadChrSubset = tadData[tadData[:,0] == gene[0]]
	
	#gene starts before tad end, ends after start
	matches = tadChrSubset[(tadChrSubset[:,2] >= gene[1]) * (tadChrSubset[:,1] <= gene[2])]
	
	if len(matches) > 1:

		continue
	filteredGenes.append(gene)

filteredGenes = np.array(filteredGenes, dtype='object')

#for the non-disrupted, the patient does not matter, so keep all patient just in 1 array.

tadPositiveAndNegativeSet = []

tadInd = 0

#We define the positive/negative TADs based on ALL SVs. the negative set is only truly negative if there is also no other SV type disrupting it. 
for tad in tadDisruptions:

	patientCount = dict()
		
	for sv in tadDisruptions[tad]:
		
		patient = sv[0]
		if patient not in patientCount:
			patientCount[patient] = []
		patientCount[patient].append(sv[1][8].svType)
	
	#get the patient names that have disrupted TADs
	patientsWithDisruptions = []
	for tadDisr in tadDisruptions[tad]:

		patientsWithDisruptions.append(tadDisr[0])

	#find the genes that are in this TAD
	
	splitTad = tad.split("_")
	geneChrMatches = filteredGenes[filteredGenes[:,0] == splitTad[0]]
	
	#the gene is inside the tad if the start is inside the tad, or if the end is inside the tad
	
	#is the start of the gene inside the tad?
	geneMatches = (geneChrMatches[:,1] <= int(splitTad[2])) * (geneChrMatches[:,2] >= int(splitTad[1]))
	
	allMatches = geneChrMatches[geneMatches]
	#add genes only when these are not affected by any mutation

	#go through all the genes and all the patients and their expression values.
	positivePatients = []
	negativePatients = []
	svTypes = []
	for gene in allMatches:
		
		#extract the row for this gene
		if gene[3].name not in expressionData[:,0]:
			continue
		
		geneExpr = expressionData[expressionData[:,0] == gene[3].name][0]
		
		#for each patient, append the expression to either the disrupted or non-disrupted based on the tad patient list
		
		for sv in tadDisruptions[tad]:

			svType = sv[1][8].svType

			patient = sv[0]
			
			if patient not in samples:
				print(patient, 'not in samples')
				continue
			patientInd = samples.index(patient)
			
			if patient not in positivePatients:
				positivePatients.append(patient)
				svTypes.append(sv[1][8].svType)
			
			
			#we filter genes in the TAD based on the SV type.
			if svType == 'DEL':
				#in case of DEL, we filter genes with any mutation. Deleted genes are not relevant

				if gene[3].name in svPatientsDel[patient] or gene[3].name in svPatientsInv[patient] or \
				gene[3].name in svPatientsItx[patient] or gene[3].name in cnvPatientsDel[patient] or \
				gene[3].name in cnvPatientsAmp[patient] or gene[3].name in snvPatients[patient] or \
				gene[3].name in svPatientsDup[patient]:

					continue

			elif svType == 'DUP':
				
				
				#in case of a DUP, we keep genes that are disrupted by the TAD disrupting DUP,
				#because those are the ones that see the effect.
				#because the CNV amp may overlap with the dup, ignore that one too. 
				if gene[3].name in svPatientsDel[patient] or gene[3].name in svPatientsInv[patient] or \
				gene[3].name in svPatientsItx[patient] or gene[3].name in cnvPatientsDel[patient] or \
				gene[3].name in snvPatients[patient]:
				#gene[3].name in snvPatients[patient] or gene[3].name in cnvPatientsAmp[patient]:
					continue
				
			elif svType == 'INV':
				#only ignore genes that are in the INV.
				if gene[3].name in svPatientsDel[patient] or gene[3].name in svPatientsDup[patient] or \
				gene[3].name in svPatientsItx[patient] or gene[3].name in cnvPatientsDel[patient] or \
				gene[3].name in cnvPatientsAmp[patient] or gene[3].name in snvPatients[patient]:
					continue
				
			else:
				#if ITX, skip all mutations. 
				if gene[3].name in svPatientsDel[patient] or gene[3].name in svPatientsInv[patient] or \
				gene[3].name in svPatientsItx[patient] or gene[3].name in cnvPatientsDel[patient] or \
				gene[3].name in cnvPatientsAmp[patient] or gene[3].name in snvPatients[patient] or \
				gene[3].name in svPatientsDup[patient]:
					continue
			
			disruptedPairs[patient][gene[3].name] = float(geneExpr[patientInd])
			

		for patientInd in range(1, len(samples)):
			
			#for the negative set, we consider all TADs that are not disrupted by
			#ANY SV.
			
			
			patient = samples[patientInd]
			
			
			#in the negative set, filter all genes that have ANY mutation. 
			if gene[3].name in svPatientsDel[patient] or gene[3].name in svPatientsInv[patient] or \
				gene[3].name in svPatientsDup[patient] or gene[3].name in svPatientsItx[patient] or \
				gene[3].name in cnvPatientsDel[patient] or gene[3].name in cnvPatientsAmp[patient] or \
				gene[3].name in snvPatients[patient]:
					continue
			
			#only if there is at least 1 gene without a mutation in this patient,
			#do we add the TAD as part of the negative set. 
			if patient not in patientsWithDisruptions:
				if patient not in negativePatients:
					negativePatients.append(patient)
	

			if patient not in patientsWithDisruptions:
				
				nonDisruptedPairs[gene[3].name][patient] = float(geneExpr[patientInd])
				
	tadPositiveAndNegativeSet.append([tad, positivePatients, negativePatients, svTypes])
	tadInd += 1

tadPositiveAndNegativeSet = np.array(tadPositiveAndNegativeSet, dtype='object')
np.savetxt(specificOutDir + '/tadPositiveAndNegativeSet.txt', tadPositiveAndNegativeSet, fmt='%s', delimiter='\t')
print(specificOutDir + '/tadPositiveAndNegativeSet.txt')
#For each gene in the disrupted group, compute the z-score of the gene compared to the expression of all patients in the negative group

#store the z-scores in an aggregated way, gene/patient pairs. Then we can get the overall plots. 
zScores = []
pValues = []
zScoresPerGene = dict()
scorePairs = [] #keep a list of the pairs that have a z-score originally, and which were set to 0 for the ranks. 
for patient in disruptedPairs:

	for gene in disruptedPairs[patient]:
		
		if gene not in zScoresPerGene:
			zScoresPerGene[gene] = dict()
			

		expr = disruptedPairs[patient][gene]

		#get the negative set
		if gene not in nonDisruptedPairs:
			continue
		negExprPatients = nonDisruptedPairs[gene]
		
		
		##if we work with GTEx as normal controls, read the negative expression from that matrix instead.
		negExpr = []
		if settings.general['gtexControl'] == True:
			
			gtexExpressionForGene = gtexExpressionData[gtexExpressionData[:,0] == gene][0]
			negExpr = [float(i) for i in gtexExpressionForGene[1:]] #skip the gene name and convert to float

			for negPatient in negExprPatients:

				negExpr.append(negExprPatients[negPatient])
			
		else:	
			
			for negPatient in negExprPatients:
				
				negExpr.append(negExprPatients[negPatient])
			
		if np.std(negExpr) == 0:
			continue

		#compute the z-score
		z = (float(expr) - np.mean(negExpr)) / float(np.std(negExpr))
		
			
		posExpr = []
		for patientName in disruptedPairs:
			
			
			if gene in disruptedPairs[patientName]:
				posExpr.append(disruptedPairs[patientName][gene])
			
		zScoresPerGene[gene][patient] = z
		scorePairs.append(patient + '_' + gene)
		
#Go through the patients and compute the ranking of the z-scores inside the list.

for gene in zScoresPerGene:
	
	patients = list(zScoresPerGene[gene].keys())
	zScores = np.array(list(zScoresPerGene[gene].values()))

	#instead of ranks, we can also add a normalized/sigmoid value of the z-score. 

	for zScoreInd in range(0, len(zScores)):
		
		z = zScores[zScoreInd]
		patient = patients[zScoreInd]
		
		pValue = stats.norm.sf(abs(z))*2
		
		pValues.append([patient + '_' + gene, pValue, z])

pValues = np.array(pValues, dtype='object')

print(pValues.shape)

#Do MTC on the pValues
reject, pAdjusted, _, _ = multipletests(pValues[:,1], method='bonferroni') #fdr_bh or bonferroni
print(reject)
print(pAdjusted)
signPatients = []
for pValueInd in range(0, len(pValues[:,1])):

	signPatients.append([pValues[pValueInd][0], pValues[pValueInd][1], pAdjusted[pValueInd], reject[pValueInd], pValues[pValueInd][1], pValues[pValueInd][2]])

signPatients = np.array(signPatients, dtype='object')

print(signPatients.shape)
prefix = ''
if settings.general['gtexControl'] == True:
	suffix = '_gtex'

if randomize == 'True':
	np.savetxt(specificOutDir + '/zScores' + prefix + '_random.txt', signPatients, fmt='%s', delimiter='\t')
else:
	np.savetxt(specificOutDir + '/zScores' + prefix + '.txt', signPatients, fmt='%s', delimiter='\t')




