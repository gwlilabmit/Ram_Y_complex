'''
Gets all unique kmers for a given value of k for all 21 sites and calculates their frequency
This is meant to be used as a helper function--I'm going to iterate through k = 3->10 and identify the most enriched kmers 
'''

def enrichment(k): 

	#################
	#input files
	#################

	cleavagesites = open("21_cleavage_final.txt", 'r') #windows for all 21 sites in subtilis

	#background datasets--200nt windows from 5' ends of all subtilis cds's, 3' ends or an equivalent number of 200nt sequences randomly found through the genome

	fiveprime = open("200nt_5prime_windows.txt", 'r') 
	threeprime = open("200nt_3prime_windows.txt", 'r')
	random = open("random_sequences.txt", 'r') 

	###############################
	#global vars and dictionaries
	###############################

	kmers = set() #a set of all unique kmers, produced by get_kmers() so I don't redo it on every iter of enrichment()--that's why it's a set, not a list

	subtiseqs = [] #an array with all subtilis sequence windows

	#arrays for all sequences in the background datasets
	fiveprimeseqs = []
	threeprimeseqs = []
	randomseqs = [] 


	subtiliscounts = {} #key: kmer, value: number of times it occurs in all subtilis sequences

	#dictionaries of the form of subtiliscounts--one for each background set

	fiveprimecounts = {}
	threeprimecounts = {}
	randomcounts = {}

	'''
	I define kmer enrichment as frequency(21 sites)/frequency(background set) for a given background set. This is based on https://www.ncbi.nlm.nih.gov/pubmed/29883606
	The enrichment dictionaries are of the form key:kmer, value:enrichment for that background set 
	'''

	fiveprimeenrichment = {}
	threeprimeenrichment = {}
	randomenrichment = {}

	##################################################
	#get all sequences for background sets
	##################################################

	#subtilis 

	for line in cleavagesites: 
		if ">" not in line: #i.e. if the line isn't a fasta header
			subtiseqs.append(line[:-1]) #the [:-1] gets rid of the \n at the end

	#5'

	for line in fiveprime: 
		if ">" not in line: #i.e. if the line isn't a fasta header
			fiveprimeseqs.append(line[:-1]) #the [:-1] gets rid of the \n at the end

	#3'

	for line in threeprime: 
		if ">" not in line: #i.e. if the line isn't a fasta header
			threeprimeseqs.append(line[:-1]) #the [:-1] gets rid of the \n at the end

	#random

	for line in random: 
		if ">" not in line: #i.e. if the line isn't a fasta header
			randomseqs.append(line[:-1]) #the [:-1] gets rid of the \n at the end

	##########################################
	#identify all unique kmers in the 21 sites
	##########################################

	#take a rolling window of size k through each sequence in subtiseqs--if it's unique add it to kmers

	for sequence in subtiseqs: 
		for i in range(len(sequence)-k+1): 
			kmer = sequence[i:i+k]
			kmers.add(kmer)

	#we don't need kmers to be a set anymore, convert to a list so we can manipulate it later
	kmers = list(kmers)

	#################################################################
	#calculate frequency of each kmer in subtilis+background windows
	#################################################################

	for kmer in kmers: 
		subtiliscounts[kmer] = 0 #add it as a key to subtiliscounts
		
		fiveprimecounts[kmer] = 0
		threeprimecounts[kmer] = 0
		randomcounts[kmer] = 0

		for sequence in subtiseqs: 
			frequency = sequence.count(kmer)
			subtiliscounts[kmer]+= frequency

		#################
		#background sets
		#################

		for sequence in fiveprimeseqs: 
			frequency = sequence.count(kmer)
			fiveprimecounts[kmer]+= frequency

		for sequence in threeprimeseqs: 
			frequency = sequence.count(kmer)
			threeprimecounts[kmer]+= frequency

		for sequence in randomseqs: 
			frequency = sequence.count(kmer)
			randomcounts[kmer]+= frequency


	###################################
	#calculate enrichment of each kmer
	###################################

	for key in subtiliscounts.keys(): #each key is a unique kmer
		
		#get the frequency for each key for 21 sites+ each background set 
		subtilis = subtiliscounts[key]
		five = fiveprimecounts[key]
		three = threeprimecounts[key]
		rand = randomcounts[key]

		#if a kmer is present in the 21 sites but not in the background sets, set its enrichment to 0--there's not enough information to determine
		#if this is the result of actual enrichment or simply being unlucky when picking sites in the background sets
	
		if float(five)>0.0: 	
			fiveprimeenrichment[key] = float(subtilis)/float(five)
		else: 
			fiveprimeenrichment[key] = 0
		
		if float(three)>0.0: 
			threeprimeenrichment[key] = float(subtilis)/float(three)
		else: 
			threeprimeenrichment[key] = 0

		if float(rand)>0.0:
			randomenrichment[key] = float(subtilis)/float(rand)
		else: 
			randomenrichment[key] = 0
		


	random.close()
	threeprime.close()
	fiveprime.close()

	cleavagesites.close()

	return fiveprimeenrichment, threeprimeenrichment, randomenrichment

'''
The output is kinda ugly but it gets the job done. 
Basically, it's of the form ({'GCG':0.002},{'GCG':0.001},{'GCG':0.003}). This is in the order 5':3':random. 

'''
	
allkmersfive = {} #enrichment relative to 5' for all values of k--formatted the same as fiveprimeenrichment [key=kmer, value=enrichment]
allkmersthree = {}  #enrichment relative to 3' for all values of k--formatted the same as threeprimeenrichment [key=kmer, value=enrichment]
allkmersrandom = {} #enrichment relative to random for all values of k--formatted the same as randomenrichment [key=kmer, value=enrichment]

for k in range(3, 11): #so we capture all kmers of size 3->10
	kenriched = enrichment(k)
	
	#extract the different dictionaries from the output of enrichment
	fiveprime = kenriched[0]
	threeprime = kenriched[1]
	random = kenriched[2]
	
	#add the values of fiveprime/threeprime/random to the corresponding "allkmers" dictionary--see https://thispointer.com/how-to-merge-two-or-more-dictionaries-in-python/
	allkmersfive.update(fiveprime)
	allkmersthree.update(threeprime)
	allkmersrandom.update(random)



#from https://pybit.es/dict-ordering.html -- orders dictionaries such that the key with largest value [so, the most frequent kmer] appears first
######print the top 10 enriched kmers

print "5' enrichment: "
print sorted(allkmersfive.items(), key=lambda x: x[1], reverse=True)[:10]
#print sorted(allkmersfive.items(), key=lambda x: x[1], reverse=True)


print "3' enrichment: "
print sorted(allkmersthree.items(), key=lambda x: x[1], reverse=True)[:10]
#print sorted(allkmersthree.items(), key=lambda x: x[1], reverse=True)


print "random enrichment: "
print sorted(allkmersrandom.items(), key=lambda x: x[1], reverse=True)[:10]  
#print sorted(allkmersrandom.items(), key=lambda x: x[1], reverse=True)  



