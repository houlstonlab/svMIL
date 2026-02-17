"""
	Generate files with dictionaries where each patient has a list of genes that are affected by a certain mutation type.
	We use this to filter out genes that are not affected just by the non-coding SV, but also by other stuff.

	Split by: SNVs, CNV amplifications (> 2.3), CNV deletions ( < 1.7), (coding) SV deletions, SV duplications, SV inversions, SV translocations.

	This works with HMF data, PCAWG data and TCGA data.

"""


import sys
import numpy as np
import os
from os import listdir
from os.path import isfile, join
import glob
import re
import gzip
path = sys.argv[1]
sys.path.insert(1, path)
sys.path.insert(1, 'linkSVsGenes/')

from inputParser import InputParser
import settings

outDir = sys.argv[2]


#Get the CNVs per gene for Nunes data
def getPatientsWithCNVGeneBased_nunes(cnvDir, cancerTypeIds):

	cnvPatientsDel = dict()
	cnvPatientsAmp = dict()

	for sampleId in cancerTypeIds:
		# For Nunes data, CNV files are named differently
		matchedFiles = glob.glob(cnvDir + '/*' + sampleId + '-T*_gene_copynumber.tsv')
		
		if len(matchedFiles) < 1:
			print('Skipping ', sampleId, ' which has no CNV data')
			continue
		
		matchedFile = matchedFiles[0]

		if sampleId not in cnvPatientsAmp:
			cnvPatientsAmp[sampleId] = []
		if sampleId not in cnvPatientsDel:
			cnvPatientsDel[sampleId] = []

		with open(matchedFile, 'r') as inF:

			lineCount = 0
			for line in inF:
				line = line.strip()

				if lineCount < 1: #skip header
					lineCount += 1
					continue

				splitLine = line.split("\t")

				gene = splitLine[3]

				if float(splitLine[5]) > 1.7 and float(splitLine[5]) < 2.3: #these are not CNVs
					continue

				if float(splitLine[5]) > 2.3:

					cnvPatientsAmp[sampleId].append(gene)
				elif float(splitLine[5]) < 1.7:

					cnvPatientsDel[sampleId].append(gene)

	return cnvPatientsAmp, cnvPatientsDel

#Get the CNVs per gene
def getPatientsWithCNVGeneBased_hmf(cnvDir, cancerTypeIds):

	cnvPatientsDel = dict()
	cnvPatientsAmp = dict()

	for sampleId in cancerTypeIds:
		patientId = cancerTypeIds[sampleId]


		matchedFile = glob.glob(cnvDir + '/*_' + patientId + '/' + sampleId + '.purple.cnv.gene.tsv')[0]

		if sampleId not in cnvPatientsAmp:
			cnvPatientsAmp[sampleId] = []
		if sampleId not in cnvPatientsDel:
			cnvPatientsDel[sampleId] = []

		with open(matchedFile, 'rb') as inF:

			lineCount = 0
			for line in inF:
				line = line.decode('ISO-8859-1')

				if lineCount < 1: #skip header
					lineCount += 1
					continue

				splitLine = line.split("\t")

				gene = splitLine[3]


				if float(splitLine[5]) > 1.7 and float(splitLine[5]) < 2.3: #these are not CNVs
					continue

				if float(splitLine[5]) > 2.3:

					cnvPatientsAmp[sampleId].append(gene)
				elif float(splitLine[5]) < 1.7:

					cnvPatientsDel[sampleId].append(gene)

	return cnvPatientsAmp, cnvPatientsDel

#Get the SNVs per gene for Nunes data
def getPatientsWithSNVs_nunes(snvDir, cancerTypeIds):

	patientsWithSNVs = dict()

	for sampleId in cancerTypeIds:
		# For Nunes data, look for SNV VCF files
		matchedFiles = glob.glob(snvDir + '/*_short_variants_genome_*' + sampleId + '-T*.vcf.gz')

		if len(matchedFiles) < 1:
			print('Skipping ', sampleId, ' which has no SNV data')
			continue

		matchedFile = matchedFiles[0]

		#open the .gz file
		with gzip.open(matchedFile, 'rb') as inF:

			for line in inF:
				line = line.strip().decode('ISO-8859-1')

				if re.search('^#', line): #skip header
					continue

				#skip the SNV if it did not pass.
				splitLine = line.split("\t")
				filterInfo = splitLine[6]
				if filterInfo not in ['PASS', '.']:
					continue

				#Check if this SNV has any affiliation with a gene. Look for gene annotation in INFO field.
				infoField = splitLine[7]

				# Look for ANN field in snpEff annotated VCF
				annSearch = re.search('ANN=([^;]+)', infoField)
				if annSearch:
					# ANN format: ALT|effect|impact|GENE_NAME|GENE_SYMBOL|...
					# Gene name is at position 3 when split by |
					annField = annSearch.group(1)
					# Split by comma to get first annotation (there can be multiple)
					firstAnnotation = annField.split(',')[0]
					# Split by pipe and get gene name (position 3)
					annParts = firstAnnotation.split('|')
					if len(annParts) > 3:
						geneName = annParts[3]
						if geneName and geneName != '':
							if sampleId not in patientsWithSNVs:
								patientsWithSNVs[sampleId] = []
							patientsWithSNVs[sampleId].append(geneName)	
						
	return patientsWithSNVs

#Get the SNVs per gene
def getPatientsWithSNVs_hmf(snvDir, cancerTypeIds):

	patientsWithSNVs = dict()

	for sampleId in cancerTypeIds:
		patientId = cancerTypeIds[sampleId]

		matchedFile = glob.glob(snvDir + '/*_' + patientId + '/' + sampleId + '.purple.somatic.vcf.gz')[0]

		#open the .gz file
		with gzip.open(matchedFile, 'rb') as inF:

			for line in inF:
				line = line.strip().decode('ISO-8859-1')

				if re.search('^#', line): #skip header
					continue

				#skip the SV if it did not pass.
				splitLine = line.split("\t")
				filterInfo = splitLine[6]
				if filterInfo not in ['PASS', '.']:
					continue

				#Check if this SNV has any affiliation with a gene. This means that in the info field, a gene is mentioned somewhere. That is, there is an ENSG identifier.
				infoField = splitLine[7]

				geneSearch = re.search('(ENSG\d+)', infoField)
				if geneSearch:
					#the gene name is always directly before the ENSG identifier
					geneMatch = re.search('.+[\|\=\(](.+)?\|ENSG\d+', infoField).group(1)


					if sampleId not in patientsWithSNVs:
						patientsWithSNVs[sampleId] = []
					patientsWithSNVs[sampleId].append(geneMatch)

	return patientsWithSNVs

#Get the SVs per gene for Nunes data
def getPatientsWithSVs_nunes(svDir, allGenes, cancerTypeIds):

	svPatientsDel = dict()
	svPatientsDup = dict()
	svPatientsInv = dict()
	svPatientsItx = dict()

	for sampleId in cancerTypeIds:

		if sampleId not in svPatientsDel:
			svPatientsDel[sampleId] = []
		if sampleId not in svPatientsDup:
			svPatientsDup[sampleId] = []
		if sampleId not in svPatientsInv:
			svPatientsInv[sampleId] = []
		if sampleId not in svPatientsItx:
			svPatientsItx[sampleId] = []

		# For Nunes data, look for SV VCF files and exclude .SV.vcf.gz files
		matchedFiles = glob.glob(svDir + '/*_sv_*' + sampleId + '-T*' + '.vcf.gz')
		matchedFiles = [f for f in matchedFiles if not f.endswith('.SV.vcf.gz')]

		if len(matchedFiles) < 1:
			print('Skipping ', sampleId, ' which has no SV data')
			continue

		matchedFile = matchedFiles[0]

		#open the .gz file
		with gzip.open(matchedFile, 'rb') as inF:

			for line in inF:
				line = line.strip().decode('ISO-8859-1')

				if re.search('^#', line): #skip header
					continue

				#skip the SV if it did not pass.
				splitLine = line.split("\t")
				filterInfo = splitLine[6]
				if filterInfo not in ['PASS', '.']:
					continue

				chr1 = splitLine[0]
				pos1 = int(splitLine[1])
				pos2Info = splitLine[4]
				infoField = splitLine[7]

				# Extract SVCLASS from INFO field
				svClassMatch = re.search('SVCLASS=([^;]+)', infoField)
				if svClassMatch:
					svClass = svClassMatch.group(1)
				else:
					svClass = None

				#match the end position and orientation. if there is no colon, this is an insertion, which we can skip.
				if not re.search(':', pos2Info):
					continue

				if re.match('[A-Z]*\[.*\:\d+\[$', pos2Info):
					o1 = '-'
					o2 = '-'
				elif re.match('[A-Z]*\].*\:\d+\]$', pos2Info):
					o1 = '+'
					o2 = '+'
				elif re.match('^\].*\:\d+\][A-Z]*', pos2Info):
					o1 = '+'
					o2 = '-'
				elif re.match('^\[.*\:\d+\[[A-Z]*', pos2Info):
					o1 = '-'
					o2 = '+'
				else:
					print('unmatched: ', pos2Info)
					print(line)
					continue

				#get the chr2 information
				chr2 = re.search('[\[\]]+(.*):(\d+).*', pos2Info).group(1)
				pos2 = int(re.search('.*\:(\d+).*', pos2Info).group(1))

				# Determine SV type using SVCLASS
				svType = ''
				if svClass:
					if svClass == 'deletion':
						svType = 'DEL'
					elif svClass == 'tandem-duplication':
						svType = 'DUP'
					elif svClass == 'inversion':
						svType = 'INV'
					elif svClass == 'translocation':
						svType = 'ITX'
					else:
						print('unknown sv type, ', line)
						continue
				else:
					print('unknown sv type, ', line)
					continue

				if svType not in ['DEL', 'DUP', 'INV', 'ITX']:
					continue

				s1 = pos1
				e1 = pos1
				s2 = pos2
				e2 = pos2
				orderedChr1 = chr1
				orderedChr2 = chr2

				#switch chromosomes if necessary
				if chr1 != chr2:
					chr1_clean = chr1.replace('chr', '')
					chr2_clean = chr2.replace('chr', '')

					if chr1_clean == 'Y' and chr2_clean == 'X':
						orderedChr1 = chr2
						orderedChr2 = chr1
					elif (chr1_clean in ['X', 'Y', 'MT']) and (chr2_clean not in ['X', 'Y', 'MT']):
						orderedChr1 = chr2
						orderedChr2 = chr1
					elif (chr1_clean not in ['X', 'Y', 'MT']) and (chr2_clean not in ['X', 'Y', 'MT']):
						chr1_num = int(chr1_clean)
						chr2_num = int(chr2_clean)
						if chr1_num > chr2_num:
							orderedChr1 = chr2
							orderedChr2 = chr1
					elif (chr1_clean in ['X', 'Y', 'MT']) and (chr2_clean in ['X', 'Y', 'MT']):
						if chr1_clean == 'Y' and chr2_clean == 'X':
							orderedChr1 = chr2
							orderedChr2 = chr1
						elif chr1_clean == 'MT' and chr2_clean in ['X', 'Y']:
							orderedChr1 = chr2
							orderedChr2 = chr1

					#always switch the coordinates as well if chromosomes are switched.
					if orderedChr1 == chr2:
						s1 = pos2
						e1 = pos2
						s2  = pos1
						e2 = pos1

				else: #if the chr are the same but the positions are reversed, change these as well.
					if pos2 < pos1:
						s1 = pos2
						e1 = pos2
						s2  = pos1
						e2 = pos1

				# Ensure chr notation is present
				chr1 = orderedChr1 if orderedChr1.startswith('chr') else 'chr' + orderedChr1
				chr2 = orderedChr2 if orderedChr2.startswith('chr') else 'chr' + orderedChr2

				#Check which genes are overlapped by this SV.
				#intrachromosomal SV
				if chr1 == chr2:

					geneChrSubset = allGenes[allGenes[:,0] == chr1]

					geneMatches = geneChrSubset[(geneChrSubset[:,1] <= e2) * (geneChrSubset[:,2] >= s1)]

					if svType == 'DEL':
						for match in geneMatches:
							svPatientsDel[sampleId].append(match[3].name)

					elif svType == 'DUP':
						for match in geneMatches:
							svPatientsDup[sampleId].append(match[3].name)
					elif svType == 'INV':
						for match in geneMatches:
							svPatientsInv[sampleId].append(match[3].name)

				else:

					#find breakpoints in the gene for each side of the SV
					geneChr1Subset = allGenes[allGenes[:,0] == chr1]
					geneChr2Subset = allGenes[allGenes[:,0] == chr2]

					#check if the bp start is within the gene.
					geneChr1Matches = geneChr1Subset[(s1 >= geneChr1Subset[:,1]) * (s1 <= geneChr1Subset[:,2])]
					geneChr2Matches = geneChr2Subset[(s2 >= geneChr2Subset[:,1]) * (s2 <= geneChr2Subset[:,2])]

					for match in geneChr1Matches:
						svPatientsItx[sampleId].append(match[3].name)

					for match in geneChr2Matches:
						svPatientsItx[sampleId].append(match[3].name)

	return svPatientsDel, svPatientsDup, svPatientsInv, svPatientsItx

# #Get the SV per gene
def getPatientsWithSVs_hmf(svDir, allGenes, cancerTypeIds):

	#Get all parsed and annotated SV type files from the main dir
	#use all genes because there is no gene in the file, so by overlap we determine which genes are affected by the SVs. 

	svPatientsDel = dict()
	svPatientsDup = dict()
	svPatientsInv = dict()
	svPatientsItx = dict()

	patientsWithSNVs = dict()

	for sampleId in cancerTypeIds:
		patientId = cancerTypeIds[sampleId]

		if sampleId not in svPatientsDel:
			svPatientsDel[sampleId] = []
		if sampleId not in svPatientsDup:
			svPatientsDup[sampleId] = []
		if sampleId not in svPatientsInv:
			svPatientsInv[sampleId] = []
		if sampleId not in svPatientsItx:
			svPatientsItx[sampleId] = []

		matchedFile = glob.glob(svDir + '/*_' + patientId + '/' + sampleId + '.purple.sv.ann.vcf.gz')[0]

		#open the .gz file
		with gzip.open(matchedFile, 'rb') as inF:

			for line in inF:
				line = line.strip().decode('ISO-8859-1')

				if re.search('^#', line): #skip header
					continue

				#skip the SV if it did not pass.
				splitLine = line.split("\t")
				filterInfo = splitLine[6]
				if filterInfo not in ['PASS', '.']:
					continue

				chr1 = splitLine[0]
				pos1 = int(splitLine[1])
				pos2Info = splitLine[4]

				#match the end position and orientation. if there is no orientation info, this is an insertion, which we can skip.
				if not re.search(':', pos2Info):
					continue

				if re.match('[A-Z]*\[.*\:\d+\[$', pos2Info):
					o1 = '+'
					o2 = '-'
				elif re.match('[A-Z]*\].*\:\d+\]$', pos2Info):
					o1 = '-'
					o2 = '+'
				elif re.match('^\].*\:\d+\][A-Z]*', pos2Info):
					o1 = '+'
					o2 = '+'
				elif re.match('^\[.*\:\d+\[[A-Z]*', pos2Info):
					o1 = '-'
					o2 = '-'
				else:
					print('unmatched: ', pos2Info)
					print(line)
					exit()

				#get the chr2 information
				chr2 = re.search('[\[\]]+(.*):(\d+).*', pos2Info).group(1)
				pos2 = int(re.search('.*\:(\d+).*', pos2Info).group(1))


				#get the SV type. This uses the rules as defined in the LUMPY paper.
				svType = ''
				if chr1 != chr2: #this is definitely a translocation.
					svType = 'ITX'
				else:
					#if the strands are the same, it is an inversion
					if o1 == o2:
						svType = 'INV'
					elif o1 == '+' and o2 == '-':
						svType = 'DEL'
					elif o1 == '-' and o2 == '+':
						svType = 'DUP'
					else:
						print('unknown sv type, ', line)
						exit(1)

				if svType not in ['DEL', 'DUP', 'INV', 'ITX']:
					continue

				chr1 = splitLine[0]
				pos1 = int(splitLine[1])
				pos2Info = splitLine[4]
				pos2 = int(re.search('.*\:(\d+).*', pos2Info).group(1))
				chr2 = re.search('[\[\]]+(.*):(\d+).*', pos2Info).group(1)

				s1 = pos1
				e1 = pos1
				s2 = pos2
				e2 = pos2
				orderedChr1 = chr1
				orderedChr2 = chr2

				#switch chromosomes if necessary
				if chr1 != chr2:
					if chr1 == 'Y' and chr2 == 'X':
						orderedChr1 = chr2
						orderedChr2 = chr1
					if (chr1 == 'X' or chr1 == 'Y' or chr1 == 'MT') and (chr2 != 'X' and chr2 != 'Y' and chr2 != 'MT'):
						orderedChr1 = chr2
						orderedChr2 = chr1
					if (chr1 != 'X' and chr1 != 'Y' and chr1 != 'MT') and (chr2 != 'X' and chr2 != 'Y' and chr2 != 'MT'):
						if int(chr1) > int(chr2):
							orderedChr1 = chr2
							orderedChr2 = chr1
					if (chr1 in ['X', 'Y', 'MT']) and (chr2 in ['X', 'Y', 'MT']): #order these as well
						if chr1 == 'Y' and chr2 == 'X':
							orderedChr1 = chr2
							orderedChr2 = chr1
						if chr1 == 'MT' and chr2 in ['X', 'Y']:
							orderedChr1 = chr2
							orderedChr2 = chr1


					#always switch the coordinates as well if chromosomes are switched.
					if orderedChr1 == chr2:
						s1 = pos2
						e1 = pos2
						s2  = pos1
						e2 = pos1

				else: #if the chr are the same but the positions are reversed, change these as well.
					if pos2 < pos1:
						s1 = pos2
						e1 = pos2
						s2  = pos1
						e2 = pos1

				chr1 = 'chr' + orderedChr1
				chr2 = 'chr' + orderedChr2

				#Check which genes are overlapped by this SV.
				#Keep track of the disrupted genes in the patient.

				#intrachromosomal SV
				if chr1 == chr2:

					geneChrSubset = allGenes[allGenes[:,0] == chr1]

					geneMatches = geneChrSubset[(geneChrSubset[:,1] <= e2) * (geneChrSubset[:,2] >= s1)]

					if svType == 'DEL':
						for match in geneMatches:
							svPatientsDel[sampleId].append(match[3].name)


					elif svType == 'DUP':
						for match in geneMatches:
							svPatientsDup[sampleId].append(match[3].name)
					elif svType == 'INV':
						for match in geneMatches:
							svPatientsInv[sampleId].append(match[3].name)

				else:

					#find breakpoints in the gene for each side of the SV
					geneChr1Subset = allGenes[allGenes[:,0] == chr1]
					geneChr2Subset = allGenes[allGenes[:,0] == chr2]

					#check if the bp start is within the gene.
					geneChr1Matches = geneChr1Subset[(s1 >= geneChr1Subset[:,1]) * (s1 <= geneChr1Subset[:,2])]
					geneChr2Matches = geneChr2Subset[(s2 >= geneChr2Subset[:,1]) * (s2 <= geneChr2Subset[:,2])]

					for match in geneChr1Matches:
						svPatientsItx[sampleId].append(match[3].name)



					for match in geneChr2Matches:
						svPatientsItx[sampleId].append(match[3].name)



	return svPatientsDel, svPatientsDup, svPatientsInv, svPatientsItx


causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.

#Combine the genes into one set.
allGenes = np.concatenate((causalGenes, nonCausalGenes), axis=0)

#get the right data based on the data input source.
if settings.general['source'] == 'Nunes':

	#read in the metadata file and get the right file identifiers
	metadataFile = settings.files['metadataNunes']

	#save the IDs of the patients with this cancer type
	cancerTypeIds = dict()
	with open(metadataFile, 'rt') as inF:

		for line in inF:
			line = line.strip()
			if not line or line.startswith('sample_id'):  # skip header and empty lines
				continue
			splitLine = line.split('\t')
			if len(splitLine) >= 2 and splitLine[1].strip() == settings.general['cancerType']:
				sampleId = splitLine[0].strip()

				#check if we have SV data for this sample
				matchedFiles = glob.glob(settings.files['svDir'] + '/*_sv_*' + sampleId + '-T*' + '.vcf.gz')
				matchedFiles = [f for f in matchedFiles if not f.endswith('.SV.vcf.gz')]

				#if we don't have SVs for this sample, skip it.
				if len(matchedFiles) < 1:
					print('Skipping ', sampleId, ' which has no SVs')
					continue

				#### check if we have expression for this sample.
				expressionDir = settings.files['expressionDir']
				matchedExpressionFiles = glob.glob(expressionDir + '/CPM_TMM.txt')

				if len(matchedExpressionFiles) < 1:
					print('Skipping ', sampleId, ' due to missing expression data')
					continue

				cancerTypeIds[sampleId] = sampleId

	svPatientsDel, svPatientsDup, svPatientsInv, svPatientsItx = getPatientsWithSVs_nunes(settings.files['svDir'], allGenes, cancerTypeIds)
	cnvPatientsAmp, cnvPatientsDel = getPatientsWithCNVGeneBased_nunes(settings.files['cnvDir'], cancerTypeIds)
	snvPatients = getPatientsWithSNVs_nunes(settings.files['snvDir'], cancerTypeIds)

elif settings.general['source'] == 'HMF':

	#from the metadata, get the right samples to use.
	#read in the metadata file and get the right file identifiers
	metadataFile = settings.files['metadataHMF']

	#save the IDs of the patients with this cancer type
	cancerTypeIds = dict()
	with open(metadataFile, 'rb') as inF:

		for line in inF:
			line = line.decode('ISO-8859-1')

			splitLine = line.split('\t')
			if splitLine[6] == settings.general['cancerType']:
				sampleId = splitLine[1]
				patientId = splitLine[0]

				#here already skip patients for which there is no SV/expression data.
				#we don't need to process these patients
				matchedFiles = glob.glob(settings.files['svDir'] + '/*_' + patientId + '/' + sampleId + '.purple.sv.ann.vcf.gz')

				#if we don't have SVs for this sample, skip it.
				if len(matchedFiles) < 1:
					print('Skipping ', sampleId, ' which has no SVs')
					continue

				#### check if we have expression for this sample.
				expressionDir = settings.files['expressionDir']
				matchedExpressionFiles = glob.glob(expressionDir + sampleId)

				if len(matchedExpressionFiles) < 1:
					print('Skipping ', sampleId, ' due to missing expression data')
					continue

				cancerTypeIds[sampleId] = patientId
	
	svPatientsDel, svPatientsDup, svPatientsInv, svPatientsItx = getPatientsWithSVs_hmf(settings.files['svDir'], allGenes, cancerTypeIds)
	cnvPatientsAmp, cnvPatientsDel = getPatientsWithCNVGeneBased_hmf(settings.files['cnvDir'], cancerTypeIds)
	snvPatients = getPatientsWithSNVs_hmf(settings.files['snvDir'], cancerTypeIds)


finalOutDir = outDir + '/patientGeneMutationPairs/'

if not os.path.exists(finalOutDir):
	os.makedirs(finalOutDir)

np.save(finalOutDir + 'svPatientsDel.npy', svPatientsDel)
np.save(finalOutDir + 'svPatientsDup.npy', svPatientsDup)
np.save(finalOutDir + 'svPatientsInv.npy', svPatientsInv)
np.save(finalOutDir + 'svPatientsItx.npy', svPatientsItx)
np.save(finalOutDir + 'cnvPatientsDel.npy', cnvPatientsDel)
np.save(finalOutDir + 'cnvPatientsAmp.npy', cnvPatientsAmp)
np.save(finalOutDir + 'snvPatients.npy', snvPatients)