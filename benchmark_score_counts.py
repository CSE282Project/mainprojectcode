import numpy as np
import time
import matplotlib.pyplot as plt

from GenerateRandomSeqWithMotifs import *
import greedy_ctcf_matchings as grdy
import ctcf_matchings as match
import ctcf_matchings_complete as cmplt

recursiveCount=0
count = 0

def individual_test(f,genome,motifs,k):
	if f=='Greedy':
		from greedy_ctcf_matchings import count
		grdy.init_counter()  # initialize counters
		matching = grdy.greedy_matching(genome, motifs, k) # run matching
		weight = matching.get_weight()[0]
		recursiveCount = 0
		count = 0
	elif f=='Dynamic':
		from ctcf_matchings import recursiveCount
		from ctcf_matchings import count
		match.initCounter()  # initialize counters
		matchings = match.maximal_matching(genome, motifs, k) # run matching
		weightsL = [matching.get_weight()[0] for matching in matchings]
		weight   = sum(weightsL)/len(weightsL)
	elif f=='Complete':
		from ctcf_matchings_complete import recursiveCount
		from ctcf_matchings_complete import count
		cmplt.initCounter()  # initialize counters
		matchings = cmplt.maximal_matching(genome, motifs, k) # run matching
		weightsL = [matching.get_weight()[0] for matching in matchings]
		weight   = sum(weightsL)/len(weightsL)

	print '-'*20
	print recursiveCount
	print count
	print weight

	return recursiveCount,count,weight


if __name__ == '__main__':
#def main():
	# Which algorithms
	algorithms = ['Greedy','Dynamic']#,'Complete']
	# Genome lengths
	#L   = [round(10**i) for i in range(3,4,.5)]
	L = [150,200,300,400,500,600]
	# n replicates of m motifs of length k
	m = [3,4,5,6,10,12,14,16]#,500,1000,5000]
	n = m
	k = 5

	# Initialize data tensors
	rec_count = np.zeros( ( len(algorithms) , len(L) , len(m) ) )
	T_n       = np.zeros( ( len(algorithms) , len(L) , len(m) ) )
	weight    = np.zeros( ( len(algorithms) , len(L) , len(m) ) )
	elapsed   = np.zeros( ( len(algorithms) , len(L) , len(m) ) )
	# Which files
	#files=['greedy_ctcf_matchings','ctcf_matchings','ctcf_matchings.complete']

	# Iterate over procedures and parameters
	for f in range(len(algorithms)):
		for i in range(len(L)):
			for j in range(len(m)): 
				print j
				before = time.time()
				genome,motifs = random_genome_wrapper(n[i],L[i],m[i],k)
				after  = time.time()
				rec_count[f,i,j], T_n[f,i,j], weight[f,i,j] = individual_test(algorithms[f],genome,motifs,k)
				elapsed[f,i,j] = after-before


print rec_count

fixL = 3
fixm = 7
plt.figure(1)
# recursive counts vs L
mat = rec_count
# red dashes, blue dashes and green dashes
plt.subplot(421)
plt.plot(L,mat[0,:,fixm], 'r--', L, mat[1,:,fixm],'b--')
plt.xlabel('Length of Genome')
plt.ylabel('Recursive Counts')
# recursive counts vs m
# red dashes, blue dashes and green dashes
plt.subplot(422)
plt.plot(m,mat[0,fixL,:], 'r--', m, mat[1,fixL,:],'b--')
plt.xlabel('Number of Forward CTCF Binding Motifs')
plt.ylabel('Recursive Counts')

# T_n vs L
mat = T_n
# red dashes, blue dashes and green dashes
plt.subplot(423)
plt.plot(L,mat[0,:,fixm], 'r--', L, mat[1,:,fixm],'b--')
plt.xlabel('Length of Genome')
plt.ylabel('T(n)')
# recursive counts vs m
# red dashes, blue dashes and green dashes
plt.subplot(424)
plt.plot(m,mat[0,fixL,:], 'r--', m, mat[1,fixL,:],'b--')
plt.xlabel('Number of Forward CTCF Binding Motifs')
plt.ylabel('T(n)')

# weights vs L
mat = weight
# red dashes, blue dashes and green dashes
plt.subplot(425)
plt.plot(L,mat[0,:,fixm], 'r--', L, mat[1,:,fixm],'b--')
plt.xlabel('Length of Genome')
plt.ylabel('Average Score')
# recursive counts vs m
# red dashes, blue dashes and green dashes
plt.subplot(426)
plt.plot(m,mat[0,fixL,:], 'r--', m, mat[1,fixL,:],'b--')
plt.xlabel('Number of Forward CTCF Binding Motifs')
plt.ylabel('Average Score')

# elapsed vs L
mat = elapsed
# red dashes, blue dashes and green dashes
plt.subplot(427)
plt.plot(L,mat[0,:,fixm], 'r--', L, mat[1,:,fixm],'b--')
plt.xlabel('Length of Genome')
plt.ylabel('Time Elapsed (sec)')
# recursive counts vs m
# red dashes, blue dashes and green dashes
plt.subplot(428)
plt.plot(m,mat[0,fixL,:], 'r--', m, mat[1,fixL,:],'b--')
plt.xlabel('Number of Forward CTCF Binding Motifs')
plt.ylabel('Time Elapsed (sec)')


plt.show()

