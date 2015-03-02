 from scipy.stats import rv_discrete

'''
randomSequence will return a random DNA sequence of length L.
L, int is the length of the output seequence
prob, optional probabilities list of lenght 4 corresponding to the probabilities of choosing A, T, G, or C.
'''
def randomSequence(L,prob=[.25,.25,.25,.25]):
	nuc = ['A','T','G','C']
	distrib = rv_discrete(values=(nuc, prob))
	sequence = ''.join(distrib.rvs(size=L))
	return sequence

def revComp(g):
	return [{'A':'T','T':'A','C':'G','G':'C'}[g[i]] for i in range(len(g)-1,-1,-1)]

'''
randomSequence will return a random DNA sequence of length L.
L, int is the length of the output sequence
k, int is the length of the motifs in the sequnece
n, int is the number of k-mer motifs to be included in the sequence
prob, optional probabilities list of lenght 4 corresponding to the probabilities of choosing A, T, G, or C.
'''
def randomSequence_motifs(L,k,n,prob=[.25,.25,.25,.25]):
	nuc = ['A','T','G','C']
	distrib = rv_discrete(values=(nuc, prob))
	# Initialize base sequence
	sequence = ''.join(distrib.rvs(size=L))

	# add motifs
	prob_motif = [1./L]*L
	for i in range(n):
		# get start position
		distrib_motif = (values=(range(L), prob_motif))
		start = ''.join(distrib_motif.rvs(size=1))
		# get motif
		motif_i = ''.join(distrib.rvs(size=k))
		# insert motif
		sequence = sequence[0:start] + motif_i + sequence[(start+k+1)::]
		# update probability distribution
		prob_motif[start:start+k] = [0]*k


random.choice
rand.arrange

# scaling of greedy and dp as a function 


