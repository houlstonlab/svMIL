"""
	Cluster ENCODE signal files using sliding windows

"""

import sys
import numpy as np
import glob

#read file

#provide source and automatically cluster the bed files
folder = sys.argv[1]

# exclude already clustered files
# dataFiles = glob.glob(folder + '*[!_clustered].bed')
# specifically choose lifted over hg38 files
dataFiles = glob.glob(folder + '*_hg38.bed')
for dataFile in dataFiles:


#
# dataFile = '../../data/histones/hmec/ENCFF154XFN_H3K27ac.bed'
# dataFile = '../../data/histones/hmec/ENCFF065TIH_H3K4me3.bed'
# dataFile = '../../data/histones/hmec/ENCFF291WFP_H3K27me3.bed'
# dataFile = '../../data/histones/hmec/ENCFF336DDM_H3K4me1.bed'
# #dataFile = 'tf_experimentallyValidated.bed'
	data = np.loadtxt(dataFile, dtype='object')
	
	# Detect number of columns
	if data.ndim == 1:
		# Single row case
		num_cols = len(data)
	else:
		num_cols = data.shape[1] if data.ndim > 1 else 3
	
	has_extra_cols = num_cols > 3

	#sort the file for this to work

	#map chr to numbers first
	chrMap = {'chr1': 1, 'chr2' : 2, 'chr3' : 3, 'chr4' : 4, 'chr5' : 5, 'chr6' : 6, 'chr7' : 7,
			  'chr8' : 8, 'chr9' : 9, 'chr10' : 10, 'chr11' : 11, 'chr12' : 12, 'chr13' : 13,
			  'chr14' : 14, 'chr15' : 15, 'chr16' : 16, 'chr17' : 17, 'chr18' : 18, 'chr19' : 19,
			  'chr20' : 20, 'chr21' : 21, 'chr22' : 22, 'chrX' : 23, 'chrY' : 24}

	mappedData = []
	for mark in data:
		if mark[0] not in chrMap:
			continue
		mappedChr = chrMap[mark[0]]
		
		if has_extra_cols:
			mappedData.append([mappedChr, int(mark[1]), int(mark[2]), mark[3], mark[4]])
		else:
			mappedData.append([mappedChr, int(mark[1]), int(mark[2])])

	mappedData = np.array(mappedData)

	# fixed to correctly interpret the data as ints not strings
	sortedData = mappedData[np.lexsort((mappedData[:,1].astype(int), mappedData[:,0].astype(int)))]

	import pandas as pd
	# changed to use sortedData instead of mappedData
	if has_extra_cols:
		df = pd.DataFrame(sortedData, columns=['chr', 'start', 'end', '.', 'signal'])
	else:
		df = pd.DataFrame(sortedData, columns=['chr', 'start', 'end'])

	# this line was not sorting correctly
	# df = df.sort_values(['chr', 'start', 'end'])

	print(df)
	data = df.to_numpy()

	#add back the chr notation
	chrData = []
	for mark in data:
		chrNotation = ''
		if mark[0] == 23:
			chrNotation = 'chrX'
		elif mark[0] == 24:
			chrNotation = 'chrY'
		else:
			chrNotation = 'chr' + str(mark[0])

		if has_extra_cols:
			chrData.append([chrNotation, mark[1], mark[2], mark[3], mark[4]])
		else:
			chrData.append([chrNotation, mark[1], mark[2]])

	chrData = np.array(chrData, dtype='object')

	data = chrData
	windowSize = 1000
	clusters = []
	currentCluster = []
	for markInd in range(0, data.shape[0]):
		mark = data[markInd,:]

		#if this is the first eQTL, it is always the start of a cluster.
		if markInd == 0:
			if has_extra_cols:
				currentCluster.append([mark[0], mark[1], mark[2], mark[3], mark[4]])
			else:
				currentCluster.append([mark[0], mark[1], mark[2]])
		else:
			#otherwise, check if the distance to the previous is within the window size.
			previousMark = data[markInd-1,:]

			# abs() avoids unordered marks (giving negative distances) from being clustered together.
			if abs(float(mark[1]) - float(previousMark[1])) <= windowSize and previousMark[0] == mark[0]:
				if has_extra_cols:
					currentCluster.append([mark[0], mark[1], mark[2], mark[3], mark[4]])
				else:
					currentCluster.append([mark[0], mark[1], mark[2]])
			else:
				clusters.append(currentCluster)
				currentCluster = []
				if has_extra_cols:
					currentCluster.append([mark[0], mark[1], mark[2], mark[3], mark[4]])
				else:
					currentCluster.append([mark[0], mark[1], mark[2]])


	#remaining eQTLs
	clusters.append(currentCluster)
	clusters = np.array(clusters, dtype='object')

	print(clusters)

	#write clusters to a new file
	clusterFile = dataFile + '_clustered.bed'
	with open(clusterFile, 'w') as outF:

		for cluster in clusters:

			clusterStart = 0
			clusterEnd = 0
			markInd = 0
			avgSignal = 0
			for mark in cluster:
				if markInd == 0:
					clusterStart = int(mark[1])
					clusterEnd = int(mark[2])
					markInd += 1
				clusterEnd = np.max([clusterEnd, int(mark[2])]) # now uses maximum end position of any mark in cluster
				if has_extra_cols:
					avgSignal += int(mark[4])

			chrom = cluster[0][0]
			
			if has_extra_cols:
				avgSignal = int(avgSignal / len(cluster))
				clusterLine = chrom + '\t' + str(clusterStart) + '\t' + str(clusterEnd) + '\t.\t' + str(avgSignal) + '\n'
			else:
				clusterLine = chrom + '\t' + str(clusterStart) + '\t' + str(clusterEnd) + '\n'
			
			outF.write(clusterLine)

