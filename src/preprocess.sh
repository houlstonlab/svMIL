### PREPROCESSING OF REGULATORY INFORMATION AND INPUT DATA###

# Some regulatory data files and other files need to be parsed to work properly as input to the tool.

#set 'false' to 'true' to run a specific block.

### PARSE CPG ISLANDS ###
run=false

if $run; then
	runFolder='./DataProcessing'
	inFile='../data/cpg/cpgIslandExt.txt'
	outFile='../data/cpg/cpg2.bed'

	python "$runFolder/cpgToBed.py" "$inFile" "$outFile"

fi

### CANCER TYPE-SPECIFIC DATA PROCESSING ###
#These are the cancer types we want to run the tool for, so preprocess their data.

# types=('hmec' 'ov' 'gm12878' 'coad' 'luad' 'urinaryTract' 'prostate' 'esophagus' 'skin' 'pancreas' 'uterus' 'nervousSystem' 'kidney')
types=('coad')

### PARSE EQTLS ###
#This requires download of GTEx_Analysis_v7.metasoft.txt.gz. See eQTLs/readme.txt.
#FOr kidney, GTEx_Analysis_v8.metasoft.txt needs to be downloaded.
run=false

if $run; then
	runFolder='./DataProcessing'
	inFile='../data/eQTLs/GTEx_Analysis_v8.metasoft.txt'
	geneLookupFile='../data/genes/ensemblGenesHg19'

	#loop over the cancer types
	for type in ${types[@]}; do
		outFolder="../data/eQTLs/$type/"

		# if [ "$type" == "kidney" ]; then
		# 	inFile='../data/eQTLs/GTEx_Analysis_v8.metasoft.txt'
		# fi

		python "$runFolder/filterEQTLs.py" "$inFile" "$geneLookupFile" "$outFolder" "$type"


	done

fi

run=true
if $run; then
	runFolder='./DataProcessing/'

	#loop over the different types
	for type in ${types[@]}; do

		#Cluster histones
		python "$runFolder/genericClustering.py" "../data/h3k27ac/$type/"

		#cluster eQTLs
		python "$runFolder/eQTLClustering.py" "../data/eQTLs/$type/"
	done

fi

#convert the SEdb files to actual bed files
run=false
types=('uterus' 'kidney' 'nervousSystem', 'prostate')
if $run; then
	runFolder='./DataProcessing'
	dataFolder='../data/superEnhancers/'

	#loop over the different types
	for type in ${types[@]}; do

		python "$runFolder/convertSeDbFiles.py" "$dataFolder/$type/" "$dataFolder/$type/$type.bed"

	done
fi



### TMM NORMALIZATION ###
run=false

if $run; then
	runFolder='./DataProcessing'
	settingsFile='./settings/settings_HMF_Breast_hmec' #which file does not matter, all we need
	#are the paths to the svDir, expressionDir and metadata file.

	#first merge the HMF expression files
	python "$runFolder/mergeHMFExpression.py" "$settingsFile"
	
	Rscript "$runFolder/normalizeEdgeR.R"

fi

