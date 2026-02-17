"""
	The goal of this script is to generate the similarity matrices for MIL
	The similarity matrices are pre-generated for the CV type that these will be used for.
	There is:
	- leave-one-chromosome-out CV
	- leave-one-patient-out CV
	- leave-bags-out CV
	- the similarity matrix on the whole dataset (used for feature importance)
	If we do feature elimination, we have similarity matrices for each feature
	output these as well.
"""


import sys
import os
import numpy as np
import pickle as pkl
import random
random.seed(785)
np.random.seed(785)

import matplotlib
matplotlib.use('Agg')

path = sys.argv[7]
sys.path.insert(1, path)

import settings

featureElimination = sys.argv[2]
leaveOnePatientOut = sys.argv[3] #make the similarity matrices for each left out patient
leaveOneChromosomeOut = sys.argv[4] #1 chromosome at a time in the test set
leaveBagsOut = sys.argv[5] #random bags in each CV fold
fullDataset = sys.argv[6] #generate sim matrix for the whole dataset.

svTypes = ['DEL', 'DUP', 'INV', 'ITX']
#svTypes = ['ALL']

outDir = sys.argv[1]
finalOutDir = outDir + '/multipleInstanceLearning/similarityMatrices/'

if not os.path.exists(finalOutDir):
	os.makedirs(finalOutDir)

if featureElimination == "True":
	featureEliminationOutDir = finalOutDir + '/featureSelection'
	if not os.path.exists(featureEliminationOutDir):
		os.makedirs(featureEliminationOutDir)

if leaveOnePatientOut == 'True':
	leaveOnePatientOutDir = finalOutDir + '/leaveOnePatientOut'
	if not os.path.exists(leaveOnePatientOutDir):
		os.makedirs(leaveOnePatientOutDir)

if leaveOneChromosomeOut == 'True':
	leaveOneChromosomeOutDir = finalOutDir + '/leaveOneChromosomeOut'
	if not os.path.exists(leaveOneChromosomeOutDir):
		os.makedirs(leaveOneChromosomeOutDir)

if leaveBagsOut == 'True':
	leaveBagsOutDir = finalOutDir + '/leaveBagsOut'
	if not os.path.exists(leaveBagsOutDir):
		os.makedirs(leaveBagsOutDir)

#input the normalized bags
with open(outDir + '/linkedSVGenePairs/normalizedBags.pkl', 'rb') as handle:
	bagDict = pkl.load(handle)

#get the information for the bag labels
degPairs = np.loadtxt(outDir + '/tadDisruptionsZScores/zScores.txt', dtype='object') #labels

print(degPairs)

print("initial number of bags: ", len(bagDict))
print('initial deg pairs: ', degPairs.shape[0])

mutDir = outDir + '/patientGeneMutationPairs/'
cnvPatientsAmp = np.load(mutDir + 'cnvPatientsAmp.npy', allow_pickle=True, encoding='latin1').item()
svPatientsDup = np.load(mutDir + 'svPatientsDup.npy', allow_pickle=True, encoding='latin1').item()
svGenePairs = np.loadtxt(outDir + '/linkedSVGenePairs/nonCoding_geneSVPairs.txt_', dtype='object')
splitSVGenePairs = []
for pair in svGenePairs:

	splitPair = pair[0].split('_')

	splitSVGenePairs.append(splitPair[0] + '_' + splitPair[7] + '_' + splitPair[12])



def getSimilarityMatrix(bags, instances, reverseBagMap):
	"""
		function to get the similarity matrix. This is mainly used to make the sim matrix for the training set.
		To make the test set sim matrix, use the function below.
		
		bags (numpy array): all bags that we use for this matrix
		instances (numpy array): all instances in the bags
		reverseBagMap (dictionary): bag index as key, instance indices as values. Used to find out which instances are in which bag. 
	"""

	bagIndices = np.arange(bags.shape[0])
	similarityMatrix = np.zeros([bags.shape[0], instances.shape[0]])

	for bagInd in range(0, bags.shape[0]):

		#Get the indices of the instances that are in this bag
		instanceIndices = reverseBagMap[bagInd]
		instanceSubset = instances[instanceIndices,:]

		#get the average of all instances in this bag
		instanceAvg = np.mean(instanceSubset, axis=0)

		#compute distance to all other instances from this bag average
		distance = np.abs(instanceAvg - instances)

		#sum the distances to get 1 similarity score
		summedDistance = np.sum(distance,axis=1)
		similarityMatrix[bagInd,:] = summedDistance

	return similarityMatrix

def getSimilarityMatrixTest(testBags, trainInstances, labels):
	"""
		function to get the similarity matrix specific for the test case.
		The instances that we map the distance to are the provided train instances.
		testBags (numpy array): all test bags that we use for this matrix
		trainInstances (numpy array): all instances in the bags of the training data, we compute distance to these instances from the test bags
		labels (list): obsolete.
	"""

	similarityMatrix = np.zeros([testBags.shape[0], trainInstances.shape[0]])

	#print(similarityMatrix.shape)

	for bagInd in range(0, testBags.shape[0]):
		#print(labels[bagInd])
		#get the average of all instances in this test patient bag
		testInstances = testBags[bagInd]

		instanceAvg = np.mean(testInstances, axis=0)

		#compute distance to all other instances from this bag average
		distance = np.abs(instanceAvg - trainInstances)

		#sum the distances to get 1 similarity score
		summedDistance = np.sum(distance,axis=1)
		#print(summedDistance)
		similarityMatrix[bagInd,:] = summedDistance

	return similarityMatrix


#Generate the similarity matrices for the SV types
for svType in svTypes:
	
	bagLabels = []
	positiveBagPairNames = []
	negativeBagPairNames = []
	positiveInstanceLabels = []
	positiveBags = []
	negativeBags = []

	#for each SV-gene pair, get the instances
	for pair in bagDict:

		#check if the SV type matches our selection
		splitPair = pair.split("_")
		shortPair = splitPair[7] + '_' + splitPair[0]

		if svType != '' and svType != 'ALL':
			if splitPair[12] != svType:
				continue

		#get the label of the bag by checking if it exists in degPairs, some pairs do not have a z-score because the gene is excluded due to mutations.
		if shortPair in degPairs[:,0]:

			#get the z-score of the pair.
			degPairInfo = degPairs[degPairs[:,0] == shortPair][0]

			#if the z-score matches this criterion, the SV-gene pair is positive
			if float(degPairInfo[5]) > 1.5 or float(degPairInfo[5]) < -1.5:
				#go through the instances of this SV-gene pair, and include only those that have gains and losses, and more than 1 instance. This should in principle not happen, but good to keep a check.
				instances = []
				for instance in bagDict[pair]:

					if instance[0] == 0 and instance[1] == 0:
						continue


					instances.append(instance)

				if len(instances) < 1:
					continue

				###Here do an extra check:
				#to make fig2A, we only look at TADs with SVs across the boundary, so those z-scores are in the set.
				#BUT some of these genes are not actually affected by the SV, since this doesn't lead to
				#regulatory elements gained/lost. SO, we need to remove those here to get the actual pairs.
				#This only goes wrong for duplications, because we also keep CNV amps that could be the same event,
				#but then the duplication does not lead to gains/losses, while the CNV amp does because it is slightly
				#longer. So if there is evidence of a cnv AMP, but no non-coding duplication linked, we can remove
				#this as a positive pair.
				if splitPair[7] not in cnvPatientsAmp:
					positiveBagPairNames.append(pair)
					positiveBags.append(instances)
				else:

					dupMatch = splitPair[0] + '_' + splitPair[7] + '_DUP'

					if splitPair[0] in cnvPatientsAmp[splitPair[7]] and dupMatch not in splitSVGenePairs:
						negativeBags.append(instances)
						negativeBagPairNames.append(pair)

					else:
						positiveBagPairNames.append(pair)
						positiveBags.append(instances)
			else: #if the z-score is anything else, this bag will be labeled negative.

				#get the right number of features per instance
				instances = []
				for instance in bagDict[pair]:

					if instance[0] == 0 and instance[1] == 0:
						continue


					instances.append(instance)
					
				if len(instances) < 1:
					continue

				negativeBags.append(instances)
				negativeBagPairNames.append(pair)


	positiveBags = np.array(positiveBags)
	negativeBags = np.array(negativeBags)
	positiveBagPairNames = np.array(positiveBagPairNames)
	negativeBagPairNames = np.array(negativeBagPairNames)

	#fail-safe in case there are not enough SVs of this type
	if positiveBags.shape[0] < 2 or negativeBags.shape[0] < 2:
		continue

	
	#add the number of instances per bag as feature to the instances
	for bag in positiveBags:
		instCount = len(bag)

		for instance in bag:
			instance.append(instCount / positiveBags.shape[0])


	for bag in negativeBags:
		instCount = len(bag)

		for instance in bag:
			instance.append(instCount / negativeBags.shape[0])

	#remove instances with no variance
	posInstances = np.vstack(positiveBags)
	negInstances = np.vstack(negativeBags)

	allInstances = np.concatenate((posInstances, negInstances))

	#remove instances with 0 variance across all instances. These are not useful for the classifier. 
	from sklearn.feature_selection import VarianceThreshold
	t = 0
	vt = VarianceThreshold(threshold=t)
	vt.fit(allInstances)
	idx = np.where(vt.variances_ > t)[0]
	badIdx = np.where(vt.variances_ <= t)[0]
	np.savetxt(finalOutDir + '/lowVarianceIdx_' + svType + '.txt', badIdx)

	newPositiveBags = []
	newNegativeBags = []
	for bag in positiveBags:
		instances = []
		for instance in bag:
			filteredInstance = []
			featureInd = 0
			for feature in instance:
				if featureInd in idx:
					filteredInstance.append(feature)
				featureInd += 1
			instances.append(filteredInstance)

		newPositiveBags.append(instances)

	for bag in negativeBags:
		instances = []
		for instance in bag:
			filteredInstance = []
			featureInd = 0
			for feature in instance:
				if featureInd in idx:
					filteredInstance.append(feature)
				featureInd += 1
			instances.append(filteredInstance)

		newNegativeBags.append(instances)

	# Filter out empty bags after variance filtering
	filteredPositiveBags = []
	filteredPositiveBagPairNames = []
	for i, bag in enumerate(newPositiveBags):
		if len(bag) > 0:  # Only keep bags with instances
			filteredPositiveBags.append(bag)
			filteredPositiveBagPairNames.append(positiveBagPairNames[i])
	
	filteredNegativeBags = []
	filteredNegativeBagPairNames = []
	for i, bag in enumerate(newNegativeBags):
		if len(bag) > 0:  # Only keep bags with instances
			filteredNegativeBags.append(bag)
			filteredNegativeBagPairNames.append(negativeBagPairNames[i])
	
	positiveBags = np.array(filteredPositiveBags)
	negativeBags = np.array(filteredNegativeBags)
	positiveBagPairNames = np.array(filteredPositiveBagPairNames)
	negativeBagPairNames = np.array(filteredNegativeBagPairNames)
	
	print(positiveBags.shape)
	
	#subsample the positive bags if there are too many
	threshold = int(settings.general['bagThreshold'])
	if positiveBags.shape[0] > threshold:
		
		random.seed(785)
			
		#subsample the positive set to the threshold
		positiveBagsSubsampled = np.random.choice(positiveBags, threshold)
	
		positiveBagsSubsampleInd = np.random.choice(np.arange(positiveBags.shape[0]), threshold)
		positiveBagsSubsampled = positiveBags[positiveBagsSubsampleInd]

		positiveBags = positiveBagsSubsampled
	
	

	print('Number of positive bags: ', positiveBags.shape)
	print('Number of negative bags: ', negativeBags.shape)

	print('Number of positive instances: ', len(positiveInstanceLabels))

	if positiveBags.shape[0] == 0 or negativeBags.shape[0] == 0:
		continue

	#subsample negative to the same number of positives. 

	random.seed(785) #somehow the global setting doesn't work in the second loop? so set it here.

	#subsample the negative set to the same number of positives.
	negativeBagsSubsampled = np.random.choice(negativeBags, positiveBags.shape[0])

	negativeBagsSubsampleInd = np.random.choice(np.arange(negativeBags.shape[0]), positiveBags.shape[0])
	negativeBagsSubsampled = negativeBags[negativeBagsSubsampleInd]

	negativeBagPairNamesSubsampled = negativeBagPairNames[negativeBagsSubsampleInd]

	#posInstances = np.vstack(positiveBags)
	#negInstances = np.vstack(negativeBagsSubsampled)

	bagPairLabelsSubsampled = np.concatenate((positiveBagPairNames, negativeBagPairNamesSubsampled))

	#save the bag pair labels for later
	np.save(finalOutDir + '/bagPairLabelsSubsampled_' + svType + '.npy', bagPairLabelsSubsampled)

	#merge the bags so that we can easily get to 1 similarity matrix and do all-to-all computations
	bagsSubsampled = np.concatenate((positiveBags, negativeBagsSubsampled))

	#assign bag labels
	bagLabelsSubsampled = np.array([1]*positiveBags.shape[0] + [0]*negativeBagsSubsampled.shape[0])

	np.save(finalOutDir + '/bagLabelsSubsampled_' + svType + '.npy', bagLabelsSubsampled)

	#stack the instances in the bags so that we can easily compute bag-instance distances
	instancesSubsampled = np.vstack(bagsSubsampled)

	#also output the instances for later
	np.save(finalOutDir + '/instancesSubsampled_' + svType + '.npy', instancesSubsampled)

	#and save the bags.
	np.save(finalOutDir + '/bagsSubsampled_' + svType + '.npy', bagsSubsampled)


	#in case of leave-one-patient out, we subsample later on
	bagPairLabels = np.concatenate((positiveBagPairNames, negativeBagPairNames))

	#save the bag pair labels for later
	np.save(finalOutDir + '/bagPairLabelsNotSubsampled_' + svType + '.npy', bagPairLabels)

	#merge the bags so that we can easily get to 1 similarity matrix and do all-to-all computations
	bags = np.concatenate((positiveBags, negativeBags))
	#assign bag labels
	bagLabels = np.array([1]*positiveBags.shape[0] + [0]*negativeBags.shape[0])

	np.save(finalOutDir + '/bagLabelsNotSubsampled_' + svType + '.npy', bagLabels)

	#stack the instances in the bags so that we can easily compute bag-instance distances
	instances = np.vstack(bags)

	#also output the instances for later
	np.save(finalOutDir + '/instancesNotSubsampled_' + svType + '.npy', instances)

	#and save the bags.
	np.save(finalOutDir + '/bagsNotSubsampled_' + svType + '.npy', bags)

	print('size comparison: ')
	print('bags not subsampled: ', bags.shape)
	print('instances not subsampled: ', instances.shape)
	print('bag labels not subsampled: ', bagLabels.shape)
	print('bag pair labels not subsampled: ', bagPairLabels.shape)

	print('bags subsampled: ', bagsSubsampled.shape)
	print('instances subsampled: ', instancesSubsampled.shape)
	print('bag labels subsampled: ', bagLabelsSubsampled.shape)
	print('bag pair labels subsampled: ', bagPairLabelsSubsampled.shape)



	# if fullDataset == 'False':
	# 	bags = bagsSubsampled
	# 	instances = instancesSubsampled
	# 	bagPairLabels = bagPairLabelsSubsampled
	# 	bagLabels = bagLabelsSubsampled
	bags = bagsSubsampled
	instances = instancesSubsampled
	bagPairLabels = bagPairLabelsSubsampled
	bagLabels = bagLabelsSubsampled



	#Make an index where we can lookup at which position the instances are in the concatenated bag array.
	reverseBagMap = dict() #lookup instance by bag index
	bagMap = dict() #lookup bag by instance index
	instanceInd = 0
	for bagInd in range(0, bags.shape[0]):
		reverseBagMap[bagInd] = []
		for instance in bags[bagInd]:
			reverseBagMap[bagInd].append(instanceInd)
			bagMap[instanceInd] = bagInd

			instanceInd += 1

	#save bagmap for later
	np.save(finalOutDir + '/bagMap_' + svType + '.npy', bagMap)

	#if we do feature elimination, randomize the features here
	featureCount = instances.shape[1]
	featureStart = featureCount-1
	if featureElimination == "True":
		featureStart = 0 #set this to featureCount to run with all features.

	#if featureStart is not updated, this will run once
	#otherwise it will randomize a new feature each time
	for featureInd in range(featureStart, featureCount):
		print('current feature: ', featureInd+1)

		if featureElimination == "True": #in case of feature elimination, we use the leave-one-chromosome-out CV setting

			#per chromosome, shuffle the features in the training set.
			#then output the original test set
			#so we have per SV type, per chromosome CV, X files for the number of features shuffled
			chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
						   'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
						   'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
						   'chr20', 'chr21', 'chr22']

			positiveBagsPerChromosome = dict()
			negativeBagsPerChromosome = dict()
			for labelInd in range(0, len(positiveBagPairNames)):
				label = positiveBagPairNames[labelInd]
				splitLabel = label.split('_')

				chromosome = splitLabel[1]
				if chromosome not in positiveBagsPerChromosome:
					positiveBagsPerChromosome[chromosome] = []
				positiveBagsPerChromosome[chromosome].append(positiveBags[labelInd])

			for labelInd in range(0, len(negativeBagPairNames)):
				label = negativeBagPairNames[labelInd]
				splitLabel = label.split('_')

				chromosome = splitLabel[1]
				if chromosome not in negativeBagsPerChromosome:
					negativeBagsPerChromosome[chromosome] = []
				negativeBagsPerChromosome[chromosome].append(negativeBags[labelInd])

			trainBags = dict()
			testBags = dict()
			trainLabels = dict()
			testLabels = dict()
			for chromosome in chromosomes:

				if chromosome not in positiveBagsPerChromosome:
					continue
				if chromosome not in negativeBagsPerChromosome:
					continue
				

				#make stratified
				testPositiveBags = positiveBagsPerChromosome[chromosome]
				testNegativeBags = negativeBagsPerChromosome[chromosome]

				testPositiveBags = np.array(testPositiveBags)
				testNegativeBags = np.array(testNegativeBags)

				
				randInd = random.sample(range(0, testNegativeBags.shape[0]), testPositiveBags.shape[0])
				testSubsetNegativeBags = testNegativeBags[randInd]

				allTestBags = []
				for bag in testPositiveBags:
					allTestBags.append(bag)
				for bag in testSubsetNegativeBags:
				#for bag in testNegativeBags:
					allTestBags.append(bag)

				allTestBags = np.array(allTestBags)

				testBags[chromosome] = allTestBags
				testLabels[chromosome] = [1]*testPositiveBags.shape[0] + [0]*testSubsetNegativeBags.shape[0]
				#testLabels[chromosome] = [1]*testPositiveBags.shape[0] + [0]*testNegativeBags.shape[0]

				testPositiveInstances = np.vstack(testPositiveBags)
				testNegativeInstances = np.vstack(testSubsetNegativeBags)
				testPositiveLabels = [1]*testPositiveInstances.shape[0]
				testNegativeLabels = [0]*testNegativeInstances.shape[0]

				#make training set from the rest
				trainingSet = []
				trainingLabels = []
				allTrainInstances = []
				allTrainLabels = []
				for chromosome2 in chromosomes:

					if chromosome == chromosome2:
						continue
					
					if chromosome2 not in positiveBagsPerChromosome:
						continue
					if chromosome2 not in negativeBagsPerChromosome:
						continue

					#make stratified
					chrPositiveBags = positiveBagsPerChromosome[chromosome2]
					chrNegativeBags = negativeBagsPerChromosome[chromosome2]

					chrPositiveBags = np.array(chrPositiveBags)
					chrNegativeBags = np.array(chrNegativeBags)

					random.seed(785)
					randInd = random.sample(range(0, chrNegativeBags.shape[0]), chrPositiveBags.shape[0])
					subsetNegativeBags = chrNegativeBags[randInd]

					for bag in chrPositiveBags:
						trainingSet.append(bag)

					for bag in subsetNegativeBags:
					#for bag in chrNegativeBags:
						trainingSet.append(bag)

					trainingLabels += [1]*chrPositiveBags.shape[0]
					trainingLabels += [0]*subsetNegativeBags.shape[0]
					#trainingLabels += [0]*chrNegativeBags.shape[0]

					trainPositiveInstances = np.vstack(chrPositiveBags)
					trainNegativeInstances = np.vstack(subsetNegativeBags)
					trainPositiveLabels = [1]*trainPositiveInstances.shape[0]
					trainNegativeLabels = [0]*trainNegativeInstances.shape[0]


				trainBags[chromosome] = np.array(trainingSet)
				trainLabels[chromosome] = trainingLabels

				trainInstances = np.vstack(trainBags[chromosome])
				
				#shuffle the training instances
				shuffledInstanceValues = trainInstances[:,featureInd]
				randomInd = np.arange(0, shuffledInstanceValues.shape[0])
				np.random.shuffle(randomInd)
				
				#we compute the similarity matrix based on the instances
				#but the instance values need to be reset every iteration
				shuffledInstances = np.zeros(trainInstances.shape)
				for col in range(0, trainInstances.shape[1]):
					if col != featureInd:
						shuffledInstances[:,col] = trainInstances[:,col]
					else:
						shuffledInstances[:,col] = trainInstances[randomInd,col]
				

				reverseBagMapOtherPatients = dict() #lookup instance by bag index
				instanceInd = 0
				for bagInd in range(0, trainBags[chromosome].shape[0]):
					reverseBagMapOtherPatients[bagInd] = []
					for instance in trainBags[chromosome][bagInd]:
						reverseBagMapOtherPatients[bagInd].append(instanceInd)
						instanceInd += 1

				similarityMatrixTrain = getSimilarityMatrix(trainBags[chromosome], shuffledInstances, reverseBagMapOtherPatients)

				#now the curent patient bags need to be to the instances of the training set
				similarityMatrixTest = getSimilarityMatrixTest(testBags[chromosome], shuffledInstances, testLabels)

				#output these to a file
				#write these data to disk so that we can access it later on
				np.save(featureEliminationOutDir + '/' + 'similarityMatrixTrain_' + svType + '_' + chromosome + '_' + str(featureInd) + '.npy', similarityMatrixTrain)
				np.save(featureEliminationOutDir + '/' + 'similarityMatrixTest_' + svType + '_' + chromosome + '_' + str(featureInd) + '.npy', similarityMatrixTest)

				#also save the labels
				np.save(featureEliminationOutDir + '/' + 'bagLabelsTrain_' + svType + '_' + chromosome + '_' + str(featureInd) + '.npy', trainLabels[chromosome])
				np.save(featureEliminationOutDir + '/' + 'bagLabelsTest_' + svType + '_' + chromosome + '_' + str(featureInd) + '.npy', testLabels[chromosome])


		elif featureElimination == 'False' and leaveOneChromosomeOut == 'True': ### leave-one-chromosome-out CV setting


			chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
						   'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
						   'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
						   'chr20', 'chr21', 'chr22']


			positiveBagsPerChromosome = dict()
			negativeBagsPerChromosome = dict()
			for labelInd in range(0, len(positiveBagPairNames)):
				label = positiveBagPairNames[labelInd]
				splitLabel = label.split('_')

				chromosome = splitLabel[1]
				if chromosome not in positiveBagsPerChromosome:
					positiveBagsPerChromosome[chromosome] = []
				positiveBagsPerChromosome[chromosome].append(positiveBags[labelInd])

			for labelInd in range(0, len(negativeBagPairNames)):
				label = negativeBagPairNames[labelInd]
				splitLabel = label.split('_')

				chromosome = splitLabel[1]
				if chromosome not in negativeBagsPerChromosome:
					negativeBagsPerChromosome[chromosome] = []
				negativeBagsPerChromosome[chromosome].append(negativeBags[labelInd])


			trainBags = dict()
			testBags = dict()
			trainLabels = dict()
			testLabels = dict()
			for chromosome in chromosomes:
				print(chromosome)

				#if chromosome != 'chr16':
				#	continue
				
				if chromosome not in positiveBagsPerChromosome:
					continue
				if chromosome not in negativeBagsPerChromosome:
					continue


				#make stratified
				testPositiveBags = positiveBagsPerChromosome[chromosome]
				testNegativeBags = negativeBagsPerChromosome[chromosome]

				testPositiveBags = np.array(testPositiveBags)
				testNegativeBags = np.array(testNegativeBags)

				random.seed(785)
				randInd = random.sample(range(0, testNegativeBags.shape[0]), testPositiveBags.shape[0])
				testSubsetNegativeBags = testNegativeBags[randInd]

				allTestBags = []
				for bag in testPositiveBags:
					allTestBags.append(bag)
				for bag in testSubsetNegativeBags:
				#for bag in testNegativeBags:
					allTestBags.append(bag)

				allTestBags = np.array(allTestBags)


				print(allTestBags.shape)


				testBags[chromosome] = allTestBags
				testLabels[chromosome] = [1]*testPositiveBags.shape[0] + [0]*testSubsetNegativeBags.shape[0]


				#make training set from the rest
				trainingSet = []
				trainingLabels = []
				allTrainInstances = []
				allTrainLabels = []
				for chromosome2 in chromosomes:

					if chromosome == chromosome2:
						continue
					
					if chromosome2 not in positiveBagsPerChromosome:
						continue
					if chromosome2 not in negativeBagsPerChromosome:
						continue

					#make stratified
					chrPositiveBags = positiveBagsPerChromosome[chromosome2]
					chrNegativeBags = negativeBagsPerChromosome[chromosome2]

					chrPositiveBags = np.array(chrPositiveBags)
					chrNegativeBags = np.array(chrNegativeBags)

					random.seed(785)
					randInd = random.sample(range(0, chrNegativeBags.shape[0]), chrPositiveBags.shape[0])
					subsetNegativeBags = chrNegativeBags[randInd]

					for bag in chrPositiveBags:
						trainingSet.append(bag)

					for bag in subsetNegativeBags:
					#for bag in chrNegativeBags:
						trainingSet.append(bag)

					trainingLabels += [1]*chrPositiveBags.shape[0]
					trainingLabels += [0]*subsetNegativeBags.shape[0]
				
				trainBags[chromosome] = np.array(trainingSet)
				trainLabels[chromosome] = trainingLabels
				

				trainInstances = np.vstack(trainBags[chromosome])

				reverseBagMapOtherPatients = dict() #lookup instance by bag index
				instanceInd = 0
				for bagInd in range(0, trainBags[chromosome].shape[0]):
					reverseBagMapOtherPatients[bagInd] = []
					for instance in trainBags[chromosome][bagInd]:
						reverseBagMapOtherPatients[bagInd].append(instanceInd)
						instanceInd += 1

				similarityMatrixTrain = getSimilarityMatrix(trainBags[chromosome], trainInstances, reverseBagMapOtherPatients)
				print(similarityMatrixTrain.shape)
				#now the curent patient bags need to be to the instances of the training set
				similarityMatrixTest = getSimilarityMatrixTest(testBags[chromosome], trainInstances, testLabels)

				#output these to a file
				#write these data to disk so that we can access it later on
				np.save(leaveOneChromosomeOutDir + '/' + 'similarityMatrixTrain_' + chromosome + '_' + svType + '.npy', similarityMatrixTrain)
				np.save(leaveOneChromosomeOutDir + '/' + 'similarityMatrixTest_' + chromosome + '_' + svType + '.npy', similarityMatrixTest)

				#also save the labels
				np.save(leaveOneChromosomeOutDir + '/' + 'bagLabelsTrain_' + chromosome + '_' + svType + '.npy', trainLabels[chromosome])
				np.save(leaveOneChromosomeOutDir + '/' + 'bagLabelsTest_' + chromosome + '_' + svType + '.npy', testLabels[chromosome])

		elif featureElimination == 'False' and leaveBagsOut == 'True': ### leave-bags-out CV setting
			
			#divide into X bags, regardless of patients
			foldSize = 10
			import math
			bagsPerFold = math.ceil((bags.shape[0] / foldSize) / 2)

			#in each fold, randomly sample positive bags and negative bags of same size
			trainBags = dict()
			testBags = dict()

			trainLabels = dict()
			testLabels = dict()

			testPairLabels = dict()
			#set random bags to use for each fold
			random.seed(785)
			randInd = random.sample(range(0, positiveBags.shape[0]), positiveBags.shape[0])
			randIndNegative = random.sample(range(0, negativeBags.shape[0]), positiveBags.shape[0])
			currentInd = 0
			currentUntil = currentInd + bagsPerFold

			for foldInd in range(0, foldSize):

				#randomly sample x positive and negative bags
				randomPositive = positiveBags[randInd[currentInd:currentUntil]]
				randomNegative = negativeBags[randInd[currentInd:currentUntil]]

				#and labels
				positiveLabels = [1]*randomPositive.shape[0]
				negativeLabels = [0]*randomNegative.shape[0]

				testBags[foldInd] = np.concatenate((randomPositive, randomNegative))
				testLabels[foldInd] = positiveLabels + negativeLabels

				#also get the pair labels
				testPairLabels[foldInd] = positiveBagPairNames[randInd[currentInd:currentUntil]]

				#then the training set will be all other bags

				otherPosInd = []
				for ind in randInd:
					if ind not in randInd[currentInd:currentUntil]:
						otherPosInd.append(ind)
				otherNegInd = []
				for ind in randInd:
					if ind not in randInd[currentInd:currentUntil]:
						otherNegInd.append(ind)

				positiveTrain = positiveBags[otherPosInd]
				negativeTrain = negativeBags[otherPosInd]

				trainPairLabels = positiveBagPairNames[otherPosInd]

				splitTrainPairs = dict()
				for pair in trainPairLabels:
					splitLabel = pair.split('_')
					splitTrainPairs[splitLabel[7] + '_' + splitLabel[0]] = pair
					
				for pair in testPairLabels[foldInd]:
					
					splitPair = pair.split('_')
					
				trainBags[foldInd] = np.concatenate((positiveTrain, negativeTrain))
				trainLabels[foldInd] = [1]*len(otherPosInd) + [0]*len(otherNegInd)

				currentInd += bagsPerFold
				if currentUntil + bagsPerFold > positiveBags.shape[0]:
					currentUntil = positiveBags.shape[0]
				else:
					currentUntil += bagsPerFold
			
			for foldInd in range(0, foldSize):
				print(foldInd)

				#get instances
				trainInstances = np.vstack(trainBags[foldInd])

				#this needs a bag map, which is changed each time we make subsets.
				reverseBagMapOtherPatients = dict() #lookup instance by bag index
				instanceInd = 0
				for bagInd in range(0, trainBags[foldInd].shape[0]):
					reverseBagMapOtherPatients[bagInd] = []
					for instance in trainBags[foldInd][bagInd]:
						reverseBagMapOtherPatients[bagInd].append(instanceInd)
						instanceInd += 1

				#collect all this information as total bags/labels

				similarityMatrixTrain = getSimilarityMatrix(trainBags[foldInd], trainInstances, reverseBagMapOtherPatients)
				print(similarityMatrixTrain.shape)
				#now the curent patient bags need to be to the instances of the training set
				similarityMatrixTest = getSimilarityMatrixTest(testBags[foldInd], trainInstances, testLabels)
				print(similarityMatrixTest.shape)
				
				np.save(leaveBagsOutDir + '/similarityMatrixTrain_' + svType + '_' + str(foldInd) + '.npy', similarityMatrixTrain)
				np.save(leaveBagsOutDir + '/similarityMatrixTest_' + svType + '_' + str(foldInd) + '.npy', similarityMatrixTest)
				np.save(leaveBagsOutDir + '/bagLabelsTrain_' + svType + '_' + str(foldInd) + '.npy', trainLabels[foldInd])
				np.save(leaveBagsOutDir + '/bagLabelsTest_' + svType + '_' + str(foldInd) + '.npy', testLabels[foldInd])

		elif featureElimination == 'False' and leaveOnePatientOut == 'True': ### leave-one-patient-out CV setting

			random.seed(785)
			#first, get the bags and labels per patient
			perPatientPositiveBags = dict()
			for bagInd in range(0, positiveBags.shape[0]):

				#get the label of this bag
				bagPairLabel = positiveBagPairNames[bagInd]
				splitLabel = bagPairLabel.split('_')

				shortPair = splitLabel[7] + '_' + splitLabel[0]

				patientId = splitLabel[7]
				if patientId not in perPatientPositiveBags:
					perPatientPositiveBags[patientId] = dict()
					perPatientPositiveBags[patientId]['bags'] = []
					perPatientPositiveBags[patientId]['pairLabels'] = []

				perPatientPositiveBags[patientId]['bags'].append(positiveBags[bagInd])
				perPatientPositiveBags[patientId]['pairLabels'].append(bagPairLabel)

			perPatientNegativeBags = dict()
			for bagInd in range(0, negativeBags.shape[0]):

				#get the label of this bag
				bagPairLabel = negativeBagPairNames[bagInd]
				splitLabel = bagPairLabel.split('_')


				patientId = splitLabel[7]
				if patientId not in perPatientNegativeBags:
					perPatientNegativeBags[patientId] = dict()
					perPatientNegativeBags[patientId]['bags'] = []
					perPatientNegativeBags[patientId]['pairLabels'] = []

				perPatientNegativeBags[patientId]['bags'].append(negativeBags[bagInd])
				perPatientNegativeBags[patientId]['pairLabels'].append(bagPairLabel)
	
			print(len(perPatientPositiveBags))
			print(len(perPatientNegativeBags))

			#for each patient, randomly subsample as many negative bags as there are positives
			perPatientBags = dict()
			skippedPatients = 0
			for patient in perPatientPositiveBags:

				if patient not in perPatientNegativeBags:
					skippedPatients += 1
					continue

				if patient not in perPatientBags:
					perPatientBags[patient] = dict()
					perPatientBags[patient]['bags'] = []
					perPatientBags[patient]['labels'] = []
					perPatientBags[patient]['pairLabels'] = []

				patientNegativeBags = perPatientNegativeBags[patient]['bags']
				patientNegativeBags = np.array(patientNegativeBags)
				
				#also get the labels
				patientNegativeBagLabels = perPatientNegativeBags[patient]['pairLabels']
				patientNegativeBagLabels = np.array(patientNegativeBagLabels)
				
				#add the same number of positives/negatives

				if len(perPatientPositiveBags[patient]['bags']) > patientNegativeBags.shape[0]:
					sampleCount = patientNegativeBags.shape[0]

					randomInd = random.sample(range(0, patientNegativeBags.shape[0]), sampleCount)
					print(patient)
					print(randomInd)
					randomNegativeBags = patientNegativeBags[randomInd]
					randomNegativeBagPairLabels = patientNegativeBagLabels[randomInd]

					for bag in randomNegativeBags:
						perPatientBags[patient]['bags'].append(bag)
						perPatientBags[patient]['labels'].append(0)
					for label in randomNegativeBagPairLabels:
						perPatientBags[patient]['pairLabels'].append(label)

					for bag in perPatientPositiveBags[patient]['bags']:
						perPatientBags[patient]['bags'].append(bag)
						perPatientBags[patient]['labels'].append(1)
					for label in perPatientPositiveBags[patient]['pairLabels']:
						perPatientBags[patient]['pairLabels'].append(label)

				else:
					sampleCount = len(perPatientPositiveBags[patient]['bags'])

					randomInd = random.sample(range(0, patientNegativeBags.shape[0]), sampleCount)
					print(patient)
					print(randomInd)
					randomNegativeBags = patientNegativeBags[randomInd]
					randomNegativeBagPairLabels = patientNegativeBagLabels[randomInd]

					for bag in randomNegativeBags:
						perPatientBags[patient]['bags'].append(bag)
						perPatientBags[patient]['labels'].append(0)
					for label in randomNegativeBagPairLabels:
						perPatientBags[patient]['pairLabels'].append(label)

					for bag in perPatientPositiveBags[patient]['bags']:
						perPatientBags[patient]['bags'].append(bag)
						perPatientBags[patient]['labels'].append(1)
					for label in perPatientPositiveBags[patient]['pairLabels']:
						perPatientBags[patient]['pairLabels'].append(label)
				
			print(skippedPatients)

			#go through each patient, and divide into train/test
			#the training set will be a merge of the bags of all other patients.
			#then for each patient, get the train/test combination
		
			testPatients = dict()
			trainPatients = dict()
			ind = 1
			foldInd = 0
			foldSize = 1
			for patient in perPatientBags:
				if foldInd not in testPatients:
					testPatients[foldInd] = []
					trainPatients[foldInd] = []

				testPatients[foldInd].append(patient)
				if ind % foldSize == 0:

					for patient2 in perPatientBags:

						if patient2 not in testPatients[foldInd]:
							trainPatients[foldInd].append(patient2)
					foldInd += 1

				ind += 1

			#add remaining patients
			if foldInd in testPatients and len(trainPatients[foldInd]) < 1:
				for patient2 in perPatientBags:
					if patient2 not in testPatients[foldInd]:
						trainPatients[foldInd].append(patient2)


			for fold in testPatients:

				testBags = []
				trainBags = []
				testLabels = []
				trainLabels = []
				testPairLabels = []
				for patient in perPatientBags:

					patientBags = perPatientBags[patient]['bags']
					patientLabels = perPatientBags[patient]['labels']
					patientPairLabels = perPatientBags[patient]['pairLabels']

					if patient in testPatients[fold]:
						testBags += patientBags
						testLabels += patientLabels
						testPairLabels += patientPairLabels
					else:
						trainBags += patientBags
						trainLabels += patientLabels

				testBags = np.array(testBags)
				trainBags = np.array(trainBags)

				#get instances
				trainInstances = np.vstack(trainBags)

				#this needs a bag map, which is changed each time we make subsets.
				reverseBagMapOtherPatients = dict() #lookup instance by bag index
				instanceInd = 0
				for bagInd in range(0, trainBags.shape[0]):
					reverseBagMapOtherPatients[bagInd] = []
					for instance in trainBags[bagInd]:
						reverseBagMapOtherPatients[bagInd].append(instanceInd)
						instanceInd += 1

				#collect all this information as total bags/labels
				similarityMatrixTrain = getSimilarityMatrix(trainBags, trainInstances, reverseBagMapOtherPatients)
				#now the curent patient bags need to be to the instances of the training set
				similarityMatrixTest = getSimilarityMatrixTest(testBags, trainInstances, testLabels)


				#write these data to disk so that we can access it later on
				np.save(leaveOnePatientOutDir + '/' + 'similarityMatrixTrain_' + str(fold) + '_' + svType + '.npy', similarityMatrixTrain)
				np.save(leaveOnePatientOutDir + '/' + 'similarityMatrixTest_' + str(fold) + '_' + svType + '.npy', similarityMatrixTest)

				#also save the labels
				np.save(leaveOnePatientOutDir + '/' + 'bagLabelsTrain_' + str(fold) + '_' + svType + '.npy', trainLabels)
				np.save(leaveOnePatientOutDir + '/' + 'bagLabelsTest_' + str(fold) + '_' + svType + '.npy', testLabels)

				#and the test pair labels
				np.save(leaveOnePatientOutDir + '/' + 'bagPairLabelsTest_' + str(fold) + '_' + svType + '.npy', testPairLabels)


		else: #output the whole similarity matrix

			#similarityMatrix = getSimilarityMatrix(bags, instances, reverseBagMap)
			similarityMatrix = getSimilarityMatrix(bags, instances, reverseBagMap)
			
			print('similarity matrix sizes: ', similarityMatrix.shape)
			print('bag size: ', bags.shape)

			#output these to a file
			#write these data to disk so that we can access it later on
			np.save(finalOutDir + '/' + 'similarityMatrix_' + svType + '.npy', similarityMatrix)