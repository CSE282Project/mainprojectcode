import numpy as np

def individual_test(f,genome,motifs,k):
	# load specific algorithm
	execfile(f)
	# initialize counters
    initCounter()
    # run matching
    matchings = maximal_matching(genome, motifs, k)

    print "Recursive Count: " + str(recursiveCount)
    print "T(n) = " + str(count)

    weightsL = [matching.get_weight() for matching in matchings]
    weight   = sum(weightsL)/len(weightsL)

    return recursiveCount,count,weight


if __name__ == '__main__':
	algorithms = ['Greedy','Dynamic','Complete']
	data = np.array(len(algo))
	files=['greedy_ctcf_matchings.py','ctcf_matchings.py','ctcf_matchings.complete.py']
	for i in files:
		for g in genomes
		RC_count, T_n, weights = individual_test(i,genome,motifs,k)