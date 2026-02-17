"""
	Create settings files automatically.
"""

import os
import glob
import sys

mainOutFolder = 'settings/'

#First define the default settings that we will use everywhere
defaultFiles = dict(
	#to run on the HMF data, we need both the directory of the SVs, and also the metadata to link file identifiers to cancer types
	svDir = '/hpc/cuppen/shared_resources/HMF_data/DR-104/data/somatics/',
	snvDir = '/hpc/cuppen/shared_resources/HMF_data/DR-104/data/somatics/',
	cnvDir = '/hpc/cuppen/shared_resources/HMF_data/DR-104/data/somatics/',
	metadataHMF = '/hpc/cuppen/shared_resources/HMF_data/DR-104/metadata/metadata.tsv',
	expressionDir = '/hpc/cuppen/shared_resources/HMF_data/DR-104/data/isofix/',
	normalizedExpressionFile = '../data/expression/HMF_TMM.txt',

	causalGenesFile = '../data/genes/CCGC.tsv', #file with all causal genes from cosmic.
	nonCausalGenesFile = '../data/genes/hg19_proteinCodingGenes.bed', #file with all protein-coding genes.
	promoterFile = '../data/promoters/epdnew_hg38ToHg19_9vC8m.bed', #File with promoters in human, not cell-specific
	cpgFile = '../data/cpg/cpgIslandExt.txt', #All human CpG sites
	tfFile = '../data/tf/tf_experimentallyValidated.bed_clustered.bed', #All validated human TFs
	chromHmmFile = '../data/chromhmm/hmec/GSE57498_HMEC_ChromHMM.bed', #Always use HMEC since we only have that data type
	rankedGeneScoreDir = "linkedSVGenePairs", #File that the output scores will be written to. The output will be in a folder with the provided UUID under this main results folder
	hg19CoordinatesFile = "../data/chromosomes/hg19Coordinates.txt",
	geneNameConversionFile = '../data/genes/allGenesAndIdsHg19.txt', #used with HMF SNVs and expression, converting ENSG identifiers to gene names.
)

#Then select a cancer type t use (this will be in the general part)
cancerType = sys.argv[1]

if cancerType == 'Colorectal':
	cancerType = 'Colon/Rectum'
elif cancerType == 'NervousSystem':
	cancerType = 'Nervous system'
elif cancerType == 'UrinaryTract':
	cancerType = 'Urinary tract'

#Create the general part
general = dict(
	source = "'HMF'",
	cancerType = '"' + cancerType + '"', #is used to get the right cancer type from the data, use the cancer type name used in the metadata.
	shuffleTads = False, #Should TAD positions be shuffled
	crd = False,
	gtexControl = False, #Should we use GTEx expression as a normal control, or use the non-disrupted TADs in other patients as control?
	geneENSG = False, #True if in the expression data the gene names are ENSG IDs. Otherwise, use False if these are regular gene names.
	bagThreshold = 700 #Subsample the bags randomly if there are more than this amount. To save computational load.
)

#Then select a cancer type for which we get the regulatory information.
#this should match the folder names used for the regulatory data. so it may differ depending
#on the tissue used
# regulatoryTypes = ['kidney', 'hmec', 'coad', 'esophagus', 'gm12878', 'luad', 'nervousSystem',
				#    'ov', 'pancreas', 'prostate', 'skin', 'urinaryTract', 'uterus']
regulatoryTypes = ['coad']
#create a settings file for each combination
for regulatoryType in regulatoryTypes:
	#regulatoryType = 'kidney'

	#determine the sub-folder in the settings folder
	subOutFolder = mainOutFolder + 'settings_HMF_' + sys.argv[1] + '_' + regulatoryType
	if not os.path.exists(subOutFolder):
		os.makedirs(subOutFolder)


	outFile = 'settings.py' #always give this the same name for the tool to recognize it

	#Write the general part to a file
	with open(subOutFolder + '/' + outFile, 'w') as outF:
		outF.write('general = dict(\n')
		for key in general:
			outF.write('\t' + key + ' = ' + str(general[key]) + ',\n')
		outF.write(')\n')


	#append the cancer-type specific part to the default part
	#use clustered for eQTLs, histone marks
	#if not available, default to whole blood.
	regulatoryElements = dict()
	#name of the file in settings as key, folder on disk as value
	regulatoryElements['tadFile'] = 'tads'
	regulatoryElements['eQTLFile'] = 'eQTLs'
	regulatoryElements['enhancerFile'] = 'enhancers'
	regulatoryElements['h3k4me3'] = 'h3k4me3'
	regulatoryElements['h3k27me3'] = 'h3k27me3'
	regulatoryElements['h3k27ac'] = 'h3k27ac'
	regulatoryElements['h3k4me1'] = 'h3k4me1'
	regulatoryElements['dnaseIFile'] = 'dnase'
	regulatoryElements['rnaPolFile'] = 'rnapol'
	regulatoryElements['superEnhancerFile'] = 'superEnhancers'
	regulatoryElements['ctcfFile'] = 'ctcf'

	regulatoryDataPath = '../data/'
	for fileType in regulatoryElements:

		#get the correct file for the selected tissue type
		folderName = regulatoryElements[fileType]
		tissueFolder = regulatoryDataPath + folderName + '/' + regulatoryType

		if folderName in ['eQTLs', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1']:
			#use the clustered.bed file
			correctFile = glob.glob(tissueFolder + '/*_clustered.bed')
		else:
			#use the regular file
			correctFile = glob.glob(tissueFolder + '/*')

		if len(correctFile) < 1: #if we are missing the file, use whole blood instead
			tissueFolder = regulatoryDataPath + folderName + '/gm12878/'
	
			if folderName in ['eQTLs', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1']:
				#use the clustered.bed file
				correctFile = glob.glob(tissueFolder + '/*_clustered.bed')
			else:
				#use the regular file
				correctFile = glob.glob(tissueFolder + '/*')
		
		#append the files to the dictionary.
		defaultFiles[fileType] = correctFile[0]
		
	#save the default files to settings as well.
	with open(subOutFolder + '/' + outFile, 'a') as outF:
		outF.write('files = dict(\n')
		for key in defaultFiles:
			outF.write('\t' + key + ' = "' + str(defaultFiles[key]) + '",\n')
		outF.write(')\n')