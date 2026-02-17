"""
	The goal of this class is to take a set of SVs as input, and then for each of these determine what the derivate TADs are on the part of the genome that they disrupt.
	In the derivative TADs, we assign gains and losses of genomic elements to the genes that are affected when the derivative TADs are made.
	
	For deletions, duplications, inversions, and (intra & interchromosomal) translocations.
	
"""

from __future__ import absolute_import
from __future__ import print_function
import re
import numpy as np
from tad import TAD
from sv import SV
from gene import Gene
from six.moves import range

class DerivativeTADMaker:
	
	
	def __init__(self, svData, genes, tadData):
		
		self.linkSVEffectsToGenes(svData, genes, tadData)
		
	
	def linkSVEffectsToGenes(self, svData, genes, tadData):
		
		"""
			For every SV, determine the type.
			For the collected SV or SVs, we first make derivative TADs, specific for each SV type.
			Then after these derivative TADs have been made, we go through the genes that are present within these new TADs, and add the newly gained or remove newly lost genomic elements for these genes.
			
			genes: (numpy array) array with the genes and their information. chr, start, end, geneObject
			svData: (numpy array) array with the SVs and their information. chr1, s1, e1, chr2, s2, e2, cancerType, sampleName, svObject.
			tadData: (numpy array) array with the TADs and their information. chr, start, end, tadObject
			
		"""
		print("Linking SV effects to genes")
		invCount = 0
		dupCount = 0
		delCount = 0

		import time
		startTime = time.time() #Keep run time of the program

		#Inversions
		print('Checking inversions')
		for sv in svData:
			
			typeMatch = re.search("inv", sv[8].svType, re.IGNORECASE)
			
			if typeMatch is not None:


				self.determineDerivativeTADs(sv, tadData, "inv") #make derivative TADs specific for this SV type
				invCount += 1
				#print("inversion count: ", invCount) #could use to check progress
		
		endTime = time.time()
		print("The program took ", endTime-startTime, " for inv")
		
		
		#Duplications
		print('Checking duplications')
		for sv in svData:

			typeMatch = re.search("dup", sv[8].svType, re.IGNORECASE)
			if typeMatch is not None:
				
				self.determineDerivativeTADs(sv, tadData, "dup")
				dupCount += 1
				#print("duplication count: ", dupCount)
		
		#Deletions
		print('Checking deletions')
		for sv in svData:
			typeMatch = re.search("del", sv[8].svType, re.IGNORECASE)
			if typeMatch is not None:
				
				self.determineDerivativeTADs(sv, tadData, "del")
				delCount += 1
				#print("deletion count: ", delCount)		
			
		#For the translocations separately	
		# 1. For each SV, determine which TAD the SVs are in
		print("matching TADs with translocations")
		tadsPerSV = self.matchTADsWithTranslocations(svData, tadData)
		
		#2. Groups are not used anymore, used to group TADs as 'chains' to determine their collective effect, but now we only look at each breakpoint individually, to model the effect we would see in
		#the clinic, in case 1 translocation is found. 
		print("making SV groups")
		svGroups = self.defineGroupsOfTranslocations(tadsPerSV)
		
		print("determine derivative TADs")
		#3. Call the derivative TAD maker on this group of SVs and let it assign the gains/losses to the genes
		import time
		startTime = time.time()
		self.determineDerivativeTADs([svGroups, tadsPerSV], tadData, "trans")
		endTime = time.time()
		print("Took ", endTime - startTime, " to determine the derivative TADs")
		
		print("done making derivative TADs")
	
	def defineGroupsOfTranslocations(self, tadsPerSV):
		"""
			This was originally meant to group translocations together, based on their 'collective' effect.
			Now it just gets the SV, and attaches the TADs that it disrupts to it.

			tadsPerSV: (dictionary) the keys of this dictionary are an SV object, the values are lists with the TADs that are affected by this SV. Output from matchTADsWithTranslocations()
			
			return
			svGroups: (numpy array) the samples and for each of those which translocations are part of that sample.
		"""

		svsPerTad = dict()
		for sv in tadsPerSV:
			if tadsPerSV[sv][0][0][3] not in svsPerTad:
				svsPerTad[tadsPerSV[sv][0][0][3]] = []
			svsPerTad[tadsPerSV[sv][0][0][3]].append(sv)
			if tadsPerSV[sv][1][0][3] not in svsPerTad:
				svsPerTad[tadsPerSV[sv][1][0][3]] = []
			svsPerTad[tadsPerSV[sv][1][0][3]].append(sv)
		
		filteredSVs = dict()
		for tad in svsPerTad:
			sv = svsPerTad[tad][0]
			filteredSVs[sv] = tadsPerSV[sv]
			
		sampleGroups = dict()
		svGroups = dict()
		samples = []
		for sv in filteredSVs:
			if sv.sampleName not in sampleGroups:
				sampleGroups[sv.sampleName] = []
			sampleGroups[sv.sampleName].append([sv.s1, sv])
		
		svGroups = []
		for sampleGroup in sampleGroups:
			currentGroup = np.array(sampleGroups[sampleGroup], dtype="object")
			currentGroupSorted = currentGroup[currentGroup[:,0].argsort()]
			samples.append(sampleGroup)
			groupList = []
			for sv in currentGroupSorted:
				groupList.append(sv[1])
			
			svGroups.append(groupList)
				
		#sort by sample name
		samples = np.array(samples)
		svGroups = np.array(svGroups)
		sortedInd = np.argsort(samples)

		return svGroups[sortedInd]	


	def matchTADsWithTranslocations(self, svData, tadData):
		"""
			For every SV, find the TAD that it ends in on the left and on the right. This can be the same TAD.
			
			svData: (numpy array) array with the SVs and their information. chr1, s1, e1, chr2, s2, e2, cancerType, sampleName, svObject.
			tadData: (numpy array) array with the TADs and their information. chr, start, end, tadObject
			
			return:
			tadsPerSV: (dictionary) dictionary with the SV object as the key, and the left TAD (in np array format, like tadData) as first array element,
									right TAD as the second. 
			
		"""
		tadsPerSV = dict()
		
		for sv in svData:
			
			#Focus on translocations only
			transTypeMatch = re.search("trans", sv[8].svType, re.IGNORECASE)
			itxTypeMatch = re.search("itx", sv[8].svType, re.IGNORECASE)
			if itxTypeMatch is None and transTypeMatch is None:
				continue
			
			#1. Check which TADs are affected by the breakpoint (there should be only 1 on each side of the SV)
			
			#Check if the SV is intrachromosomal or interchromosomal
			if sv[0] == sv[3]:
				tadChrSubset = tadData[tadData[:,0] == sv[0]]
				
				startMatches = (sv[1] > tadChrSubset[:,1]) * (sv[1] < tadChrSubset[:,2])
				matchingTadStart = tadChrSubset[startMatches]
				
				if len(matchingTadStart) < 1:
					continue #if there is no TAD affected, we skip it for now
			
				endMatches = (sv[5] > tadChrSubset[:,1]) * (sv[5] < tadChrSubset[:,2])
				matchingTadEnd = tadChrSubset[endMatches]
			
				if len(matchingTadEnd) < 1:
					continue
			
				tadsPerSV[sv[8]] = [matchingTadStart, matchingTadEnd]
			else: #interchromosomal SV
				tadChr1Subset = tadData[tadData[:,0] == sv[0]]
				
				startMatches = (sv[1] > tadChr1Subset[:,1]) * (sv[1] < tadChr1Subset[:,2])
				matchingTadStart = tadChr1Subset[startMatches]
				
				if len(matchingTadStart) < 1:
					continue #if there is no TAD affected, we skip it for now
				
				tadChr2Subset = tadData[tadData[:,0] == sv[3]]
			
				endMatches = (sv[5] > tadChr2Subset[:,1]) * (sv[5] < tadChr2Subset[:,2])
				matchingTadEnd = tadChr2Subset[endMatches]
			
				if len(matchingTadEnd) < 1:
					continue
			
				tadsPerSV[sv[8]] = [matchingTadStart, matchingTadEnd]
			
		return tadsPerSV
		
		
	def determineDerivativeTADs(self, svData, tadData, svType):
		
	
		"""
			Given an SV or a set of SVs, depending on the type of SVs, we compute how the affected region of the genome will look after the SV.
			We then make a set of new TAD objects that are located next to each other, and update all elements that are now inside these new/affected TADs. 
		
			svData: (numpy array) array with the SVs and their information. chr1, s1, e1, chr2, s2, e2, cancerType, sampleName, svObject. 
			tadData: (numpy array) array with the TADs and their information. chr, start, end, tadObject
			svType: (string) type of SV that we should determine the derivative TAD for. Either del, inv, dup or trans.
			
		
		"""

		### TRANSLOCATIONS ###
		if svType == "trans":
			
			svGroups = svData[0]
			tadsPerSV = svData[1]
			
			#1. For each translocation, determine by orientation how the TADs are disrupted
			#2. Collect the affected genes and elements
			#3. Set the gains and losses correctly for the genes. 
			
			updatedTadPos = dict() #keep the TADs and the new start/ends after 
			for group in svGroups: #This is not a group but just a sample now
				
				
				gains = dict()
				losses = dict()
				
				for sv in group: #for each SV of the sample
					
					leftTad = tadsPerSV[sv][0][0][3]
					rightTad = tadsPerSV[sv][1][0][3]
						
					
					#These are the scenarios that we need to incorporate for translocations:
					#For +-, genes in the left side of chr1 gain elements from the right side of chr2, and vice versa. 
					#For -+, genes in the left side of chr2 gain elements from the right side of chr1, and vice versa. 
					#For ++, genes in the left side of chr1 gain elements from the inverted left side of chr2, so starting from the breakpoint until the TAD start. Vice versa for the genes in the other part. 
					#For --, genes in the right side of chr1 gain elements from the inverted right side of chr2, so starting from the breakpoint until the TAD end. Vice versa for the genes in the other part.

					if sv.o1 == "+" and sv.o2 == "-": 
						
						#2. Get the elements on the left of the breakpoint and the right of the breakpoint
						leftSideElements = leftTad.getElementsByRange(leftTad.start, sv.s1)
						leftSideGenes = leftTad.getGenesByRange(leftTad.start, sv.s1)
						
						rightSideElements = rightTad.getElementsByRange(sv.e2, rightTad.end)
						rightSideGenes = rightTad.getGenesByRange(sv.e2, rightTad.end)
						
						# #Also get the elements in the remaining part
						remainingElementsLeft = leftTad.getElementsByRange(sv.s1, leftTad.end)
						remainingGenesLeft = leftTad.getGenesByRange(sv.s1, leftTad.end)
						
						remainingElementsRight = rightTad.getElementsByRange(rightTad.start, sv.e2)
						remainingGenesRight = rightTad.getGenesByRange(rightTad.start, sv.e2)
						
					elif sv.o1 == "-" and sv.o2 == "+":
				
						#2. Get the elements on the left of the breakpoint and the right of the breakpoint
						#Left side is the first chromosome, right side the second chromosome. 
						leftSideElements = leftTad.getElementsByRange(sv.s1, leftTad.end)
						leftSideGenes = leftTad.getGenesByRange(sv.s1, leftTad.end)
						
						rightSideElements = rightTad.getElementsByRange(rightTad.start, sv.e2)
						rightSideGenes = rightTad.getGenesByRange(rightTad.start, sv.e2)

						# #Also get the elements in the remaining part
						remainingElementsLeft = leftTad.getElementsByRange(leftTad.start, sv.s1)
						remainingGenesLeft = leftTad.getGenesByRange(leftTad.start, sv.s1)
						
						remainingElementsRight = rightTad.getElementsByRange(sv.e2, rightTad.end)
						remainingGenesRight = rightTad.getGenesByRange(sv.e2, rightTad.end)
						
						
					elif sv.o1 == "+" and sv.o2 == "+":
					
						leftSideElements = leftTad.getElementsByRange(leftTad.start, sv.s1)
						leftSideGenes = leftTad.getGenesByRange(leftTad.start, sv.s1)
						
						#This part is inverted, so we start from the SV until the TAD start
						rightSideElements = rightTad.getElementsByRange(rightTad.start, sv.e2)
						rightSideGenes = rightTad.getGenesByRange(rightTad.start, sv.e2)
						
						# #Also get the elements in the remaining part
						#Left side is the first chromosome, right side the second chromosome.
						remainingElementsLeft = leftTad.getElementsByRange(sv.s1, leftTad.end)
						remainingGenesLeft = leftTad.getGenesByRange(sv.s1, leftTad.end)
						
						remainingElementsRight = rightTad.getElementsByRange(rightTad.start, sv.e2)
						remainingGenesRight = rightTad.getGenesByRange(rightTad.start, sv.e2)
	
					elif sv.o1 == "-" and sv.o2 == "-":
						

						#This is the left part of chr1
						leftSideElements = leftTad.getElementsByRange(sv.s1, leftTad.end)
						leftSideGenes = leftTad.getGenesByRange(sv.s1, leftTad.end)
						
						#This part is inverted, so we start SV until the TAD end
						rightSideElements = rightTad.getElementsByRange(sv.e2, rightTad.end)
						rightSideGenes = rightTad.getGenesByRange(sv.e2, rightTad.end)
						
						
						# #Also get the elements in the remaining part
						#Left side is the first chromosome, right side the second chromosome.
						remainingElementsLeft = leftTad.getElementsByRange(leftTad.start, sv.s1)
						remainingGenesLeft = leftTad.getGenesByRange(leftTad.start, sv.s1)
						
						remainingElementsRight = rightTad.getElementsByRange(rightTad.start, sv.e2)
						remainingGenesRight = rightTad.getGenesByRange(rightTad.start, sv.e2)
					
					else:
						
						print('sv without correct orientation: ')
						print(sv, sv.o1, sv.o2)
					
					svStr = sv.chr1 + "_" + str(sv.s1) + "_" + str(sv.e1) + "_" + sv.chr2 + "_" + str(sv.s2) + "_" + str(sv.e2) + "_" + sv.sampleName + "_" + str(leftTad.startStrength) + '_' + str(leftTad.endStrength) + '_' + str(rightTad.startStrength) + '_' + str(rightTad.endStrength) + '_' + sv.svType + "_" + str(leftTad.startStrengthSignal) + '_' + str(leftTad.endStrengthSignal) + '_' + str(rightTad.startStrengthSignal) + '_' + str(rightTad.endStrengthSignal)
			
						
					for gene in leftSideGenes:		
						gene.addLostElements(remainingElementsLeft, sv.sampleName)
						gene.addLostElementsSVs(remainingElementsLeft, svStr)
						
						gene.addGainedElements(rightSideElements, sv.sampleName)
						gene.addGainedElementsSVs(rightSideElements, svStr)
						
					for gene in rightSideGenes:		
						gene.addLostElements(remainingElementsRight, sv.sampleName)
						gene.addLostElementsSVs(remainingElementsRight, svStr)
					
						gene.addGainedElements(leftSideElements, sv.sampleName)
						gene.addGainedElementsSVs(leftSideElements, svStr)

			
		### DELETIONS ###
		if svType == "del":
			
			#1. For every SV, determine the TADs that it occurs in.
			#2. Collect all elements within the SV region of these TADs
			#3. Assign the elements as lost for these genes (ignore genes that are deleted themselves). If the elements were linked to genes and are lost, or if these are TAD-wide, is determined
			#by the gene object itself. 
			
			#Determine all overlapping TADs.
			
			tadChrSubsetInd = svData[0] == tadData[:,0]
			tadChrSubset = tadData[tadChrSubsetInd]
			
			#If the SV start is before the end of the TAD, and the SV end after the start of the TAD, the TAD is overlapped.
			startMatches = svData[1] <= tadChrSubset[:,2]
			endMatches = svData[5] >= tadChrSubset[:,1]
			
			tadMatches = tadChrSubset[startMatches * endMatches]

			if tadMatches.shape[0] < 2: #no matches, or overlapping just 1 TAD. 
				return

			#Get the leftmost and rightmost TADs
			farLeftTad = tadMatches[0] #This list is sorted
			farRightTad = tadMatches[tadMatches.shape[0]-1]

			#The genes in the far left TAD, only in the part that is not overlapped by the deletion, gain the elements that are not overlapped by the deletion
			#in the far right tad.
			remainingLeftGenes = farLeftTad[3].getGenesByRange(farLeftTad[1], svData[1])
			remainingLeftElements = farLeftTad[3].getElementsByRange(farLeftTad[1], svData[1])
			
			remainingRightGenes = farRightTad[3].getGenesByRange(svData[5], farRightTad[2])
			remainingRightElements = farRightTad[3].getElementsByRange(svData[5], farRightTad[2])
			
			deletedLeftGenes = farLeftTad[3].getGenesByRange(svData[1], farLeftTad[2])
			deletedLeftElements = farLeftTad[3].getElementsByRange(svData[1], farLeftTad[2])
			
			deletedRightGenes = farRightTad[3].getGenesByRange(farRightTad[1], svData[5])
			deletedRightElements = farRightTad[3].getElementsByRange(farRightTad[1], svData[5])
			
			#for the losses, the remaining genes in the left lose the lost left elements
			#for the gains, the remaining genes in the left gain the remaining elements on the right.
			svStr = svData[0] + "_" + str(svData[1]) + "_" + str(svData[2]) + "_" + svData[3] + "_" + str(svData[4]) + "_" + str(svData[5]) + "_" + svData[8].sampleName + "_" + str(farLeftTad[3].startStrength) + '_' + str(farLeftTad[3].endStrength) + '_' + str(farRightTad[3].startStrength) + '_' + str(farRightTad[3].endStrength) + '_' + svData[8].svType + '_'  + str(farLeftTad[3].startStrengthSignal) + '_' + str(farLeftTad[3].endStrengthSignal) + '_' + str(farRightTad[3].startStrengthSignal) + '_' + str(farRightTad[3].endStrengthSignal)
			
			for gene in remainingLeftGenes: #add gains from the right
				
				if len(remainingRightElements) > 0:
					gene.addGainedElements(remainingRightElements, svData[8].sampleName)
					gene.addGainedElementsSVs(remainingRightElements, svStr)

				####Have removed losses, because these are not as a result of the TAD disruption!!! Simply because they are removed by the deletion. 
				#gene.addLostElements(deletedLeftElements, svData[8].sampleName)
				#gene.addLostElementsSVs(deletedLeftElements, svData[0] + "_" + str(svData[1]) + "_" + str(svData[2]) + "_" + svData[3] + "_" + str(svData[4]) + "_" + str(svData[5]) + "_" + svData[8].sampleName + "_" + svData[8].svType)
				
			for gene in remainingRightGenes: #add gains from the left
				
				if len(remainingRightElements) > 0:
					gene.addGainedElements(remainingLeftElements, svData[8].sampleName)
					gene.addGainedElementsSVs(remainingLeftElements, svStr)
			
				####Have removed losses, because these are not as a result of the TAD disruption!!! Simply because they are removed by the deletion. 
				#gene.addLostElements(deletedRightElements, svData[8].sampleName)
				#gene.addLostElementsSVs(deletedRightElements, svData[0] + "_" + str(svData[1]) + "_" + str(svData[2]) + "_" + svData[3] + "_" + str(svData[4]) + "_" + str(svData[5]) + "_" + svData[8].sampleName + "_" + svData[8].svType)
			
		
		### INVERSION ###
		if svType == "inv":
			
			#1. Get the two TADs at the start and end of the inversions
			
			tadChrSubsetInd = svData[0] == tadData[:,0]
			tadChrSubset = tadData[tadChrSubsetInd]
			
			svStr = svData[0] + '_' + str(svData[1]) + '_' + str(svData[2]) + '_' + svData[3] + '_' + str(svData[4]) + '_' + str(svData[5]) + '_' + str(svData[7]) + '_' + svData[8].svType
			
			#Get all TADs overlapped by the inversion
			#First get all TADs overlapped by the start of the inversion
			startMatches = svData[1] <= tadChrSubset[:,2]
			endMatches = svData[1] >= tadChrSubset[:,1]
			
			invStartMatches = startMatches * endMatches #either the start or end needs to match
			
			#Then get the TADs overlapped by the end of the inversion
			startMatches = svData[5] >= tadChrSubset[:,1]
			endMatches = svData[5] <= tadChrSubset[:,2]
			
			invEndMatches = startMatches * endMatches
			
			leftMostTad = tadChrSubset[invStartMatches]
			rightMostTad = tadChrSubset[invEndMatches]
			

			
			
			if len(leftMostTad) < 1 or len(rightMostTad) < 1:
				return #for now skip all inversions that do not end in a TAD on either side. 

			#These are only the cases where the inversion ends in a TAD on both sides. 
			if len(leftMostTad) > 0 and len(rightMostTad) > 0:

				if leftMostTad[0][1] == rightMostTad[0][1] and leftMostTad[0][2] == rightMostTad[0][2]: #skip if the SV is within a TAD entirely
					
					return
			
				
				leftMostTad = leftMostTad[0]
				rightMostTad = rightMostTad[0]
				
				
				
				#2. Collect all elements until the right TAD boundary inside the inversion.
				
				leftSideElements = leftMostTad[3].getElementsByRange(svData[1], leftMostTad[2]) #From the start of the inversion until the end of the left most TAD
				unaffectedElementsLeft = leftMostTad[3].getElementsByRange(leftMostTad[1], svData[1])
	
				#Also get the genes
				leftSideGenes = leftMostTad[3].getGenesByRange(svData[1], leftMostTad[2]) #From the start of the inversion until the end of the left most TAD
				unaffectedGenesLeft = leftMostTad[3].getGenesByRange(leftMostTad[1], svData[1])
				
				#3. Collect all elements from the left TAD boundary until the end of the inversion.
				
				rightSideElements = rightMostTad[3].getElementsByRange(rightMostTad[1], svData[5]) #from the start of the rightmost TAD until the end of the inversion
				
				unaffectedElementsRight = rightMostTad[3].getElementsByRange(svData[5], rightMostTad[2])
				
				rightSideGenes = rightMostTad[3].getGenesByRange(rightMostTad[1], svData[5]) #from the start of the rightmost TAD until the end of the inversion
				unaffectedGenesRight = rightMostTad[3].getGenesByRange(svData[5], rightMostTad[2])
				
			#Assigning the gains and losses to the genes is independent of the type of inversion
			#All genes that were originally in the left TAD (outisde of the inversion) will gain elements of the right side of the inversion
			#All unaffected genes on the left will lose the eQTLs that are in the left side of the inversion
			svStr = svData[0] + "_" + str(svData[1]) + "_" + str(svData[2]) + "_" + svData[3] + "_" + str(svData[4]) + "_" + str(svData[5]) + "_" + svData[8].sampleName + "_" + str(leftMostTad[3].startStrength) + '_' + str(leftMostTad[3].endStrength) + '_' + str(rightMostTad[3].startStrength) + '_' + str(rightMostTad[3].endStrength) + '_' + svData[8].svType  + "_" + str(leftMostTad[3].startStrengthSignal) + '_' + str(leftMostTad[3].endStrengthSignal) + '_' + str(rightMostTad[3].startStrengthSignal) + '_' + str(rightMostTad[3].endStrengthSignal)
			
			for gene in unaffectedGenesLeft:
				
				gene.addGainedElements(rightSideElements, svData[7])
				gene.addGainedElementsSVs(rightSideElements, svStr)
				
				gene.addLostElements(leftSideElements, svData[7])
				gene.addLostElementsSVs(leftSideElements, svStr)
				
			#All genes in the right side of the inversion will gain elements from the original left TAD.
			#All genes in the right side will lose interactions with eQTLs in the unaffected right TAD.
			for gene in rightSideGenes:
				
				gene.addGainedElements(unaffectedElementsLeft, svData[7])
				gene.addGainedElementsSVs(unaffectedElementsLeft, svStr)
				#print "Number of unaffected elements right: ", len(unaffectedElementsRight), " for genes ", len(rightSideGenes)
				gene.addLostElements(unaffectedElementsRight, svData[7])
				gene.addLostElementsSVs(unaffectedElementsRight, svStr)
			
			#vice versa but then for the right TAD and right side of the inversion.
			#The lost eQTLs are the ones that are in the right side of the inversion
			for gene in unaffectedGenesRight:
				
				gene.addGainedElements(leftSideElements, svData[7])
				gene.addGainedElementsSVs(leftSideElements, svStr)
				
				gene.addLostElements(rightSideElements, svData[7])
				gene.addLostElementsSVs(rightSideElements, svStr)
			
			#The lost eQTLs are the ones that are in the unaffected original left TAD
			for gene in leftSideGenes:
				
				gene.addGainedElements(unaffectedElementsRight, svData[7])
				gene.addGainedElementsSVs(unaffectedElementsRight, svStr)
				
				gene.addLostElements(unaffectedElementsLeft, svData[7])
				gene.addLostElementsSVs(unaffectedElementsLeft, svStr)

			return
		
		### DUPLICATION ###
		if svType == "dup":

			#1. Determine which TADs are involved in the duplication (only the outmost 2 are affected, the rest can be kept in tact)
			tadChrSubsetInd = svData[0] == tadData[:,0]
			tadChrSubset = tadData[tadChrSubsetInd]
			tadChrSubset = tadChrSubset[tadChrSubset[:,1].argsort()]
			
			startMatches = svData[1] < tadChrSubset[:,2]
			endMatches = svData[5] > tadChrSubset[:,1]  
			
			matches = startMatches * endMatches
	
			matchingTads = tadChrSubset[matches]
			
			
			#Remove all matches where the SV is exclusively within a TAD
			filteredTads = []
			for tad in matchingTads:
				if svData[1] > tad[1] and svData[5] < tad[2]:
					continue
				filteredTads.append(tad)

			if len(filteredTads) < 1:
				return
			
			#2. Make the derivative positions for the TAD boundaries
			
			#The original TAD boundaries of the left TAD remains where it is.
			#Then get the insert point in the right TAD (dup end position)
			#Derive the position of the duplicated boundary (piece of left TAD until TAD boundary)
			
			#The start position of the new TAD is the end of the leftmost TAD in which the duplication starts.
			#If the duplication overlaps with the start of the first TAD, so there is no TAD on the left, this should be the start of rightmost tad in which the duplication ends. 
	
			
			#In case only 1 boundary is overlapped and there is no TAD to the right, the first new TAD is also the last TAD.
			
			#Possible cases:
			#- The duplication ends in a TAD on both sides.
			#- We skip cases where the duplication does not end in a TAD on either side. 
			
			if len(filteredTads) > 1: #The SV spans multiple TADs
				
				#First get all TADs overlapped by the start of the inversion
				startMatches = svData[1] <= tadChrSubset[:,2]
				endMatches = svData[1] >= tadChrSubset[:,1]
				
				dupStartMatches = startMatches * endMatches #either the start or end needs to match
				
				#Then get the TADs overlapped by the end of the inversion
				startMatches = svData[5] >= tadChrSubset[:,1]
				endMatches = svData[5] <= tadChrSubset[:,2]
				
				dupEndMatches = startMatches * endMatches
				
				leftMostTad = tadChrSubset[dupStartMatches]
				rightMostTad = tadChrSubset[dupEndMatches]
				
				if len(leftMostTad) < 1 or len(rightMostTad) < 1:
					return #For now here only use cases where the duplication ends in 2 TADs. 

				#Assign the gained elements to the TAD.

				#Assign the elements to the new TADs in the right order.
				#The first TAD gets the eQTLs within the SV of the last TAD.
				#The last TAD gets the eQTLs within the SV of the last TAD.
				
				svInteractionsFirstTad = leftMostTad[0][3].getElementsByRange(svData[1], leftMostTad[0][2])
				svInteractionsLastTad = rightMostTad[0][3].getElementsByRange(rightMostTad[0][1], svData[5])
				
				#Determine the gains for every gene. Also for the copied TADs, there are now multiple of these genes. 
				
				#For the new TADs, this is the same principle as for the eQTLs.
				#For the duplicated TADs, we can do * 2 of the elements
				
				#For TAD 1, the first part of C can interact with the second half of A.
				svGenesFirstTad = leftMostTad[0][3].getGenesByRange(svData[1], leftMostTad[0][2])
				svGenesLastTad = rightMostTad[0][3].getGenesByRange(rightMostTad[0][1], svData[5])
				
				svStr = svData[0] + "_" + str(svData[1]) + "_" + str(svData[2]) + "_" + svData[3] + "_" + str(svData[4]) + "_" + str(svData[5]) + "_" + svData[8].sampleName + "_" + str(leftMostTad[0][3].startStrength) + '_' + str(leftMostTad[0][3].endStrength) + '_' + str(rightMostTad[0][3].startStrength) + '_' + str(rightMostTad[0][3].endStrength) + '_' + svData[8].svType  + "_" + str(leftMostTad[0][3].startStrengthSignal) + '_' + str(leftMostTad[0][3].endStrengthSignal) + '_' + str(rightMostTad[0][3].startStrengthSignal) + '_' + str(rightMostTad[0][3].endStrengthSignal)
				
				# DEBUG: Check if this is an IGF2 duplication to save debug info
				igf2_genes = [g for g in svGenesFirstTad if g.name == 'IGF2']
				if len(igf2_genes) > 0 and len(svInteractionsLastTad) > 0:
					debug_output = []
					debug_output.append("="*80)
					debug_output.append(f"DEBUG: IGF2 DUPLICATION DETECTED")
					debug_output.append(f"Patient: {svData[8].sampleName}")
					debug_output.append(f"SV: {svData[0]}:{svData[1]}-{svData[5]}")
					debug_output.append(f"Left TAD: {leftMostTad[0][0]}:{leftMostTad[0][1]}-{leftMostTad[0][2]}")
					debug_output.append(f"Right TAD: {rightMostTad[0][0]}:{rightMostTad[0][1]}-{rightMostTad[0][2]}")
					debug_output.append("")
					debug_output.append(f"svGenesFirstTad (genes in left TAD): {len(svGenesFirstTad)} genes")
					for g in svGenesFirstTad:
						debug_output.append(f"  - {g.name} ({g.chromosome}:{g.start}-{g.end})")
					debug_output.append("")
					debug_output.append(f"svInteractionsLastTad (elements from right TAD): {len(svInteractionsLastTad)} elements")
					
					# Count elements by type
					element_counts = {}
					enhancer_examples = []
					for elem in svInteractionsLastTad:
						elem_type = elem[3]
						if elem_type not in element_counts:
							element_counts[elem_type] = 0
						element_counts[elem_type] += 1
						if elem_type == 'enhancer' and len(enhancer_examples) < 10:
							# elem format: [chr, start, end, type, gene_name, ...]
							enhancer_examples.append(f"    {elem[0]}:{elem[1]}-{elem[2]} -> {elem[4] if len(elem) > 4 else 'NA'}")
					
					debug_output.append(f"Element counts by type:")
					for elem_type, count in sorted(element_counts.items()):
						debug_output.append(f"  {elem_type}: {count}")
					
					if len(enhancer_examples) > 0:
						debug_output.append("")
						debug_output.append("Example enhancers (first 10):")
						debug_output.extend(enhancer_examples)
					
					# Save to file before adding to genes
					with open('/data/scratch/DGE/DUDGE/MOPOPGEN/tyates/perturb/svMIL/debug_igf2_dup_BEFORE.txt', 'a') as f:
						f.write('\n'.join(debug_output) + '\n')
				
				for gene in svGenesFirstTad:
					
					gene.addGainedElements(svInteractionsLastTad, svData[7])
					gene.addGainedElementsSVs(svInteractionsLastTad, svStr)
					
					# DEBUG: Save IGF2's gained elements after adding
					if gene.name == 'IGF2' and len(svInteractionsLastTad) > 0:
						debug_output = []
						debug_output.append("")
						debug_output.append(f"AFTER addGainedElementsSVs for IGF2:")
						debug_output.append(f"svStr: {svStr}")
						
						if svStr in gene.gainedElementsSVs:
							debug_output.append(f"gene.gainedElementsSVs[svStr] keys: {list(gene.gainedElementsSVs[svStr].keys())}")
							debug_output.append(f"Counts:")
							for elem_type, count in gene.gainedElementsSVs[svStr].items():
								debug_output.append(f"  {elem_type}: {count}")
						else:
							debug_output.append(f"WARNING: svStr not found in gene.gainedElementsSVs!")
							debug_output.append(f"Available keys: {list(gene.gainedElementsSVs.keys())}")
						
						debug_output.append("")
						debug_output.append(f"gene.gainedElements[{svData[7]}]:")
						if svData[7] in gene.gainedElements:
							debug_output.append(f"Keys: {list(gene.gainedElements[svData[7]].keys())}")
							for elem_type, count in gene.gainedElements[svData[7]].items():
								debug_output.append(f"  {elem_type}: {count}")
						else:
							debug_output.append(f"WARNING: Sample not found in gene.gainedElements!")
						
						debug_output.append("="*80)
						debug_output.append("")
						
						with open('/data/scratch/DGE/DUDGE/MOPOPGEN/tyates/perturb/svMIL/debug_igf2_dup_AFTER.txt', 'a') as f:
							f.write('\n'.join(debug_output) + '\n')
				
				for gene in svGenesLastTad:
				
					gene.addGainedElements(svInteractionsFirstTad, svData[7])
					gene.addGainedElementsSVs(svInteractionsFirstTad, svStr)
				
				