'''
Gets all kmers of a specified length k
'''

from itertools import product

kmersize = 4

#Generate kmers as per this approach: https://stackoverflow.com/questions/48677692/generating-all-possible-k-mers-string-combinations-from-a-given-list
kmers = list(product('ATCG', repeat=kmersize))

#This gives you a list of tuples. To convert to strings, use https://stackoverflow.com/questions/19641579/python-convert-tuple-to-string
for i in range(len(kmers)): 
	kmers[i] = "".join(kmers[i])

print kmers 


