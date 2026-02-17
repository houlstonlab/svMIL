from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import settings
import sys
import os
from six.moves import range
import re
import pickle as pkl

class GeneRanking:
	"""
		Class responsible for collecting all gains and losses of genes, and writing these to an output file.
		It used to do a ranking, but doesn't do that anymore.
	
	"""
	
	def __init__(self, genes, svData, runId, permutationRound):
		"""
			genes: (numpy array) array with the genes and their information. chr	start	end	Gene (object)
			svData: (numpy array) array with the SVs, as output by InputParser
			runId: (str) the uuid obtained by main.py
			permutationRound: (str) the permutation number, appended to the output file.

			
		"""
		self.scoreGenes(genes, svData, runId, permutationRound)
	
			
	def scoreGenes(self, genes, svData, runId, permutationRound):
		"""
			Doesn't score genes anymore, but instead assigns their gains and losses as 'scores'. 
				
			genes: (numpy array) array with the genes and their information. chr, start, end, geneObject
			svData: (numpy array) array with the SVs and their information. chr1, s1, e1, chr2, s2, e2, cancerType, sampleName, svObject
			runId: (str) the uuid obtained by main.py
			permutationRound: (str) the permutation number, appended to the output file.
		"""
		

		geneMap = dict() #not sure if this is still useful, used to get back the gene from a sorting.
		reverseGeneMap = dict() #also keep a map where we can search by index to later obtain back the gene

		#Make the gene maps
		geneIndex = 0
		for gene in genes:
			if gene not in geneMap:
				geneMap[gene] = geneIndex
				geneIndex += 1
		for gene in geneMap:
			index = geneMap[gene]
			reverseGeneMap[index] = gene


		print("collecting all gains and losses")
		# features = ['eQTL', 'enhancer', 'promoter', 'cpg', 'tf', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'dnaseI', 'rnaPol',
						# 'CTCF', 'CTCF+Enhancer', 'CTCF+Promoter', 'Enhancer', 'Heterochromatin', 'Poised_Promoter', 'Promoter', 'Repeat', 'Repressed', 'Transcribed', 'superEnhancer', 'ctcf']
		features = ['eQTL', 'enhancer', 'promoter', 'cpg', 'tf', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'dnaseI', 'rnaPol',
						'superEnhancer', 'ctcf',
						'1_TssA', '2_TssFlnk', '3_TssFlnkU', '4_TssFlnkD', '5_Tx', '6_TxWk', 
						'7_EnhG1', '8_EnhG2', '9_EnhA1', '10_EnhA2', '11_EnhWk', '12_ZNF/Rpts', 
						'13_Het', '14_TssBiv', '15_EnhBiv', '16_ReprPC', '17_ReprPCWk', '18_Quies']

		# features = ['eQTL', 'enhancer', 'promoter', 'cpg', 'tf', 'hic', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3', 'dnaseI', 'rnaPol',
		# 				'CTCF', 'CTCF+Enhancer', 'CTCF+Promoter', 'Enhancer', 'Heterochromatin', 'Poised_Promoter', 'Promoter', 'Repeat', 'Repressed', 'Transcribed', 'superEnhancer', 'ctcf']
		#



		svGeneMap = dict()
		svGeneIndices = []
		allLossScores = []
		allGainScores = []
		header = 'SV_metadata' #header for output file
		for feature in features:

			lossScores, svGeneMap, svGeneIndices = self.scoreByElementLossesSVs(genes, svGeneMap, svGeneIndices, feature)
			allLossScores.append(lossScores)

			gainScores, svGeneMap, svGeneIndices = self.scoreByElementGainsSVs(genes, svGeneMap, svGeneIndices, feature)
			allGainScores.append(gainScores)

		# Build header to match actual data layout: all losses first, then all gains
		for feature in features:
			header += '\t'
			header += feature + '_loss'
		for feature in features:
			header += '\t'
			header += feature + '_gain'

		#get the strength features
		strengthFeatures = ['enhancer', 'ctcf', 'rnaPol', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1']

		#strengthFeatures = ['enhancer', 'ctcf', 'rnaPol', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3']

		# DEBUG: Check if IGF2 enhancer gains are in allGainScores
		enhancer_idx = features.index('enhancer')
		if enhancer_idx < len(allGainScores):
			igf2_pairs_in_gains = [k for k in allGainScores[enhancer_idx].keys() if k.startswith('IGF2_')]
			with open('/data/scratch/DGE/DUDGE/MOPOPGEN/tyates/perturb/svMIL/debug_igf2_geneRanking.txt', 'a') as f:
				f.write(f"\n{'='*80}\n")
				f.write(f"DEBUG: After collecting all gains/losses\n")
				f.write(f"enhancer feature index: {enhancer_idx}\n")
				f.write(f"IGF2 pairs in allGainScores[{enhancer_idx}]: {len(igf2_pairs_in_gains)}\n")
				for pair in igf2_pairs_in_gains:
					f.write(f"  {pair}: {allGainScores[enhancer_idx][pair]}\n")
				f.write(f"Total svGeneIndices: {len(svGeneIndices)}\n")
				igf2_in_indices = [idx for idx, pair in enumerate(svGeneIndices) if pair.startswith('IGF2_')]
				f.write(f"IGF2 pairs in svGeneIndices: {len(igf2_in_indices)}\n")
				for idx in igf2_in_indices[:5]:  # Show first 5
					f.write(f"  Index {idx}: {svGeneIndices[idx]}\n")
				f.write(f"{'='*80}\n\n")

		allStrengthLossScores = []
		allStrengthGainScores = []
		for feature in strengthFeatures:
			lossScores = self.scoreByElementLossesStrengthsSVs(genes, feature)
			allStrengthLossScores.append(lossScores)

			gainScores = self.scoreByElementGainsStrengthsSVs(genes, feature)
			allStrengthGainScores.append(gainScores)

		# Add strength features to header
		for feature in strengthFeatures:
			header += '\t'
			header += feature + '_strength_loss'
		for feature in strengthFeatures:
			header += '\t'
			header += feature + '_strength_gain'


		delCount = 0
		dupCount = 0
		invCount = 0
		itxCount = 0
		pairScores = np.zeros([len(svGeneIndices), 80]) #46
		pairIds = []
		
		# DEBUG: Track IGF2 pairs during scoring
		debug_igf2_tracking = []
		
		for ind in range(0, len(svGeneIndices)):
			sv = svGeneIndices[ind]

			pairIds.append(sv)
			for featureInd in range(0, len(features)):
				if sv in allLossScores[featureInd]:
					pairScores[ind,featureInd] = allLossScores[featureInd][sv]
					# DEBUG: Track IGF2 losses
					if sv.startswith('IGF2_') and features[featureInd] == 'enhancer':
						debug_igf2_tracking.append(f"IGF2 Loss: ind={ind}, featureInd={featureInd}, col={featureInd}, value={allLossScores[featureInd][sv]}, sv={sv}")

				if sv in allGainScores[featureInd]:
					pairScores[ind,featureInd+len(features)] = allGainScores[featureInd][sv]
					# DEBUG: Track IGF2 gains
					if sv.startswith('IGF2_') and features[featureInd] == 'enhancer':
						debug_igf2_tracking.append(f"IGF2 Gain: ind={ind}, featureInd={featureInd}, col={featureInd+len(features)}, value={allGainScores[featureInd][sv]}, sv={sv}")

			for featureInd in range(0, len(strengthFeatures)):
				if sv in allStrengthLossScores[featureInd]:
					pairScores[ind,(featureInd+len(features)*2)] = allStrengthLossScores[featureInd][sv]

				if sv in allStrengthGainScores[featureInd]:
					pairScores[ind,(featureInd+len(features)*2)+len(strengthFeatures)] = allStrengthGainScores[featureInd][sv]

			# translocation = 0
			# deletion = 0
			# duplication = 0
			# inversion = 0
			# #Add some features for SV type
			# splitSV = sv.split("_")
			# svType = splitSV[12]
			# if re.search("ITX", svType, re.IGNORECASE):
			# 	translocation = 1
			# 	itxCount += 1
			# elif re.search("DEL", svType, re.IGNORECASE):
			# 	deletion = 1
			# 	delCount += 1
			# elif re.search("DUP", svType, re.IGNORECASE):
			# 	duplication = 1
			# 	dupCount += 1
			# elif re.search("INV", svType, re.IGNORECASE):
			# 	inversion = 1
			# 	invCount += 1
			# 
			# pairScores[ind,70] = deletion
			# pairScores[ind,71] = duplication
			# pairScores[ind,72] = inversion
			# pairScores[ind,73] = translocation

			#any other features beyond this point were removed, because they are not useful anymore.



		pairIds = np.array(pairIds)

		outDir = sys.argv[5] + '/' + settings.files['rankedGeneScoreDir']

		if not os.path.exists(outDir):
			os.makedirs(outDir)
		if not os.path.exists(outDir + "/" + runId):
			os.makedirs(outDir + "/" + runId)

		#keep at 80 to ensure that it still works with downstream processing although final features are now gone. Filter those out later where necessary.
		pairScoresWithPairIds = np.empty([len(svGeneIndices), 81], dtype="object") #47
		#pairScoresWithPairIds = np.empty([len(svGeneIndices), 1], dtype="object")
		pairScoresWithPairIds[:,0] = pairIds
		pairScoresWithPairIds[:,1:81] = pairScores #47

		# DEBUG: Check IGF2 rows before writing
		with open('/data/scratch/DGE/DUDGE/MOPOPGEN/tyates/perturb/svMIL/debug_igf2_geneRanking.txt', 'a') as f:
			f.write(f"\n{'='*80}\n")
			f.write(f"DEBUG: Before writing to file\n")
			f.write(f"Number of features: {len(features)}\n")
			f.write(f"Enhancer gain column should be: {features.index('enhancer') + len(features) + 1} (adding 1 for pairId column)\n")
			f.write(f"\nIGF2 tracking during loop:\n")
			for line in debug_igf2_tracking:
				f.write(f"  {line}\n")
			f.write(f"\nIGF2 rows in final array:\n")
			for idx, pair_id in enumerate(pairIds):
				if str(pair_id).startswith('IGF2_'):
					enhancer_loss_col = features.index('enhancer') + 1  # +1 for pairId column
					enhancer_gain_col = features.index('enhancer') + len(features) + 1  # +1 for pairId column
					f.write(f"  Row {idx}: {pair_id}\n")
					f.write(f"    enhancer_loss (col {enhancer_loss_col}): {pairScoresWithPairIds[idx, enhancer_loss_col]}\n")
					f.write(f"    enhancer_gain (col {enhancer_gain_col}): {pairScoresWithPairIds[idx, enhancer_gain_col]}\n")
					if idx < 3:  # Show first few full rows
						f.write(f"    Full row: {pairScoresWithPairIds[idx, :]}\n")
			f.write(f"{'='*80}\n\n")

		np.savetxt(outDir + '/' + runId + "/nonCoding_geneSVPairs.txt_" + str(permutationRound), pairScoresWithPairIds, delimiter='\t', fmt='%s', header = header)

		###3. Make the feature file for MIL for each sv-gene pair
		#Each SV-gene pair is a bag. A bag can contain a variable set of isntances, which represent the gained/lost elements
		#The feature vector was pre-defined in the gene class for each instance.

		bags = dict()
		for geneInd in range(0, genes.shape[0]):
			gene = reverseGeneMap[geneInd] #Get the gene back from the scoring matrix by index

			geneSVSamples = []
			for sv in gene.SVs:
				splitSV = sv.split("_")
				geneSVSamples.append(splitSV[len(splitSV)-1])

			for sv in gene.alteredElements:

				instances = []
				for element in gene.alteredElements[sv]:
					instances.append(gene.alteredElements[sv][element])
					
				#add number of instances as final feature
				#for instance in instances:
				#	instance.append(len(instances))
				
				if len(instances) > 0:
					bags[gene.name + "_" + sv] = instances

		print(len(bags))
	
		#output the bags to a file
		with open(outDir + '/' + runId + '/bags.pkl', 'wb') as handle:
			pkl.dump(bags, handle, protocol=pkl.HIGHEST_PROTOCOL)

	
	def scoreByElementLossesSVs(self, genes, svGeneMap, svGeneIndices, elementType):
		"""
			Determine for each SV-gene pair how many elements are lost.
			
			genes:  (numpy array) array with the genes and their information. chr, start, end, geneObject
			svGeneMap: (dictionary) each sample is a key, and the value is the index of where the pair is stored in the final output
			svGeneIndices: (list) list of each pair to determine the current order of pairs for the output.
			elementType: (string) type of the element that we should score the losses of. 
			
			return
			pairScores (dictionary): a 0 or 1 for each pair key to show if it lost ANY of this element type or not.
			svGeneMap: (dictionary) each sample is a key, and the value is the index of where the pair is stored in the final output
			svGeneIndices: (list) list of each pair to determine the current order of pairs for the output.
		"""
		
		#Idea: have sv-gene pairs on one axis, and a 1 or 0 for this particular feature in the first column
		pairScores = dict() #Keep a dictionary because we do not know how large the final matrix will be across all features
		for geneInd in range(0, len(genes)):
			gene = genes[geneInd]

			geneSVSamples = []
			for sv in gene.SVs:
				splitSV = sv.split("_")
				geneSVSamples.append(sv)
			
			
			if len(gene.lostElementsSVs) > 0:

				for sv in gene.lostElementsSVs:
					pairId = gene.name + "_" + sv

					geneSVSamples = []
					for geneSV in gene.SVs:
						splitSV = geneSV.split("_")
						geneSVSamples.append(splitSV[len(splitSV)-1])


					if pairId not in svGeneMap: #Check if we already used this index for a different feature
						svInd = len(svGeneIndices)
						svGeneMap[pairId] = svInd
						svGeneIndices.append(pairId)

				
					loss = False #Make sure that we only focus on lost elements of the provided type.
					for element in gene.lostElementsSVs[sv]:
						if element == elementType:
							loss = True
					if loss == True:

						pairScores[pairId] = 1 #assume that each SV can disrupt a gene only once

		return pairScores, svGeneMap, svGeneIndices

	def scoreByElementGainsSVs(self, genes, svGeneMap, svGeneIndices, elementType):
		"""
			Determine for each SV-gene pair how many elements are gained.

			genes:  (numpy array) array with the genes and their information. chr, start, end, geneObject
			svGeneMap: (dictionary) each sample is a key, and the value is the index of where the pair is stored in the final output
			svGeneIndices: (list) list of each pair to determine the current order of pairs for the output.
			elementType: (string) type of the element that we should score the gains of.

			return
			pairScores (dictionary): a 0 or 1 for each pair key to show if it gained ANY of this element type or not.
			svGeneMap: (dictionary) each sample is a key, and the value is the index of where the pair is stored in the final output
			svGeneIndices: (list) list of each pair to determine the current order of pairs for the output.
		"""
		
		#Idea: have sv-gene pairs on one axis, and a 1 or 0 for this particular feature in the first column
		pairScores = dict() #Keep a dictionary because we do not know how large the final matrix will be across all features
		for geneInd in range(0, len(genes)):
			gene = genes[geneInd]

			# DEBUG: Check if this is IGF2 when processing enhancers
			if gene.name == 'IGF2' and elementType == 'enhancer':
				debug_output = []
				debug_output.append("="*80)
				debug_output.append(f"DEBUG: scoreByElementGainsSVs for IGF2 and enhancers")
				debug_output.append(f"len(gene.gainedElementsSVs): {len(gene.gainedElementsSVs)}")
				
				if len(gene.gainedElementsSVs) > 0:
					debug_output.append(f"Number of SVs with gained elements: {len(gene.gainedElementsSVs)}")
					for sv_key in gene.gainedElementsSVs:
						debug_output.append(f"\nSV: {sv_key}")
						debug_output.append(f"  Element types in gainedElementsSVs[sv]:")
						for elem_type, count in gene.gainedElementsSVs[sv_key].items():
							debug_output.append(f"    {elem_type}: {count}")
						
						# Check if this SV will match
						pairId = gene.name + "_" + sv_key
						debug_output.append(f"  pairId would be: {pairId}")
						debug_output.append(f"  Will be added to pairScores? {'enhancer' in gene.gainedElementsSVs[sv_key]}")
				else:
					debug_output.append("WARNING: IGF2 has empty gainedElementsSVs dict!")
				
				debug_output.append("="*80)
				debug_output.append("")
				
				with open('/data/scratch/DGE/DUDGE/MOPOPGEN/tyates/perturb/svMIL/debug_igf2_geneRanking.txt', 'a') as f:
					f.write('\n'.join(debug_output) + '\n')

			if len(gene.gainedElementsSVs) > 0:

				for sv in gene.gainedElementsSVs:
					pairId = gene.name + "_" + sv

					geneSVSamples = []
					for geneSV in gene.SVs:
						
						splitSV = geneSV.split("_")
						geneSVSamples.append(geneSV)
						
					if pairId not in svGeneMap: #Check if we already used this index for a different feature
						svInd = len(svGeneIndices)
						svGeneMap[pairId] = svInd
						svGeneIndices.append(pairId)
					
					gain = False #Make sure that we only focus on lost elements of the provided type. 
					for element in gene.gainedElementsSVs[sv]:
						if element == elementType:
							gain = True
					if gain == True:
						pairScores[pairId] = 1 #assume that each SV can disrupt a gene only once
						
						# DEBUG: Log when IGF2 enhancer gain is added
						if gene.name == 'IGF2' and elementType == 'enhancer':
							with open('/data/scratch/DGE/DUDGE/MOPOPGEN/tyates/perturb/svMIL/debug_igf2_geneRanking.txt', 'a') as f:
								f.write(f"SUCCESS: Added {pairId} to pairScores with value 1\n\n")

		
		return pairScores, svGeneMap, svGeneIndices	

	def scoreByElementGainsStrengthsSVs(self, genes, elementType):
		"""
			For each SV-gene pair, determine the strengths of the elements that are gained. It just uses the strength of one element as assigned earlier on in gene.py when obtaining gains/losses due to
			TAD disruptions. 
			
			NOTE: these features are not used anymore later on!
			
			genes:  (numpy array) array with the genes and their information. chr, start, end, geneObject
			elementType: (string) type of the element that we should score the gains of. 
			
			return
			pairScores (dictionary): a strength score for each pair key for this element type that was gained. 
		"""
		
		pairScores = dict() #Keep a dictionary because we do not know how large the final matrix will be across all features
		for geneInd in range(0, len(genes)):
			gene = genes[geneInd]
			
			if len(gene.gainedElementsStrengthsSVs) > 0:
				
				for sv in gene.gainedElementsStrengthsSVs:
					pairId = gene.name + "_" + sv
					
					for element in gene.gainedElementsStrengthsSVs[sv]:
						if element == elementType:
							
							pairScores[pairId] = gene.gainedElementsStrengthsSVs[sv][element]
		return pairScores 
		
	def scoreByElementLossesStrengthsSVs(self, genes, elementType):
		"""
			For each SV-gene pair, determine the strengths of the elements that are lost. It just uses the strength of one element as assigned earlier on in gene.py when obtaining gains/losses due to
			TAD disruptions. 
			
			NOTE: these features are not used anymore later on!
			
			genes:  (numpy array) array with the genes and their information. chr, start, end, geneObject
			elementType: (string) type of the element that we should score the gains of. 
			
			return
			pairScores (dictionary): a strength score for each pair key for this element type that was lost. 
		"""
		
		pairScores = dict() #Keep a dictionary because we do not know how large the final matrix will be across all features
		for geneInd in range(0, len(genes)):
			gene = genes[geneInd]
			
			if len(gene.lostElementsStrengthsSVs) > 0:
				
				for sv in gene.lostElementsStrengthsSVs:
					pairId = gene.name + "_" + sv
					
					for element in gene.lostElementsStrengthsSVs[sv]:
						if element == elementType:
							
							pairScores[pairId] = gene.lostElementsStrengthsSVs[sv][element]
		return pairScores 	
