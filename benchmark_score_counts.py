import numpy as np
import time
import matplotlib.pyplot as plt

from GenerateRandomSeqWithMotifs import *
import greedy_ctcf_matchings as grdy
import greedy_matching2 as grdy2
import ctcf_matchings as match
import ctcf_matchings_complete as cmplt

#recursiveCount=0
#count = 0

def individual_test(f,genome,motifs,k):
	if f=='Greedy':
		grdy.init_counter()  # initialize counters
		matching = grdy.greedy_matching(genome, motifs, k) # run matching
		weight = matching.get_weight()[0]
		from ctcf_matchings import recursiveCount
		from greedy_ctcf_matchings import count
	elif f=='Dynamic_Programming':
		match.initCounter()  # initialize counters
		matchings = match.maximal_matching(genome, motifs, k) # run matching
		weightsL = [matching.get_weight()[0] for matching in matchings]
		weight   = sum(weightsL)/len(weightsL)
		from ctcf_matchings import recursiveCount
		from ctcf_matchings import count
	elif f=='Complete':
		cmplt.initCounter()  # initialize counters
		matchings = cmplt.maximal_matching(genome, motifs, k) # run matching
		weightsL = [matching.get_weight()[0] for matching in matchings]
		weight   = sum(weightsL)/len(weightsL)
		from ctcf_matchings_complete import recursiveCount
		from ctcf_matchings_complete import count
	if f=='Greedy2':
		grdy.init_counter()  # initialize counters
		matching = grdy.greedy_matching(genome, motifs, k) # run matching
		weight = matching.get_weight()[0]
		from ctcf_matchings import recursiveCount
		from greedy_ctcf_matchings import count

	print '-'*20
	print recursiveCount
	print count
	print weight

	return recursiveCount,count,weight

def plotin3d(fig,sub,x,y,Z,xlab,ylab,title):
	from mpl_toolkits.mplot3d import Axes3D
	from matplotlib import cm
	from matplotlib.ticker import LinearLocator, FormatStrFormatter
	import matplotlib.pyplot as plt
	import numpy as np

	#fig = plt.figure()
	fig.add_subplot(sub)
	ax = fig.gca(projection='3d')
	#X = np.arange(-5, 5, 0.25)
	#Y = np.arange(-5, 5, 0.25)
	X = np.array(x)
	Y = np.array(y)
	X, Y = np.meshgrid(X, Y)
	#R = np.sqrt(X**2 + Y**2)
	#Z = np.sin(R)
	surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
	        linewidth=0, antialiased=False)
	ax.set_zlim(np.min(Z)-abs(np.min(Z)), np.max(Z)+abs(np.max(Z)))

	ax.zaxis.set_major_locator(LinearLocator(10))
	ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

	fig.colorbar(surf, shrink=0.5, aspect=5)

	ax.set_xlabel(xlab)
	ax.set_ylabel(ylab)
	ax.set_title(title)

	#plt.show()


if __name__ == '__main__':
	test = False
	d2 = True

	# Which algorithms
	algorithms = ['Greedy','Dynamic_Programming', 'Greedy2']#,'Complete']
	# Genome lengths
	#L   = [round(10**i) for i in range(3,4,.5)]
	if test:
		L = [150,200]
		m = [3,4]
	else:
		L = [500,600]
		# n replicates of m motifs of length k
		m = [3  ,4  ,5  ,6  ,10 ,12 ,14]#,500,1000,5000]
	n = m
	k = 5

	# Initialize data tensors
	rec_count = np.zeros( ( len(algorithms) , len(L) , len(m) ) )
	T_n       = np.zeros( ( len(algorithms) , len(L) , len(m) ) )
	weight    = np.zeros( ( len(algorithms) , len(L) , len(m) ) )
	elapsed   = np.zeros( ( len(algorithms) , len(L) , len(m) ) )
	out = np.zeros( (len(algorithms)*len(L)*len(m) , 7) )
	# Which files
	#files=['greedy_ctcf_matchings','ctcf_matchings','ctcf_matchings.complete']

	outcount=0
	# Iterate over procedures and parameters
	for f in range(len(algorithms)):
		for i in range(len(L)):
			for j in range(len(m)): 
				print j
				genome,motifs = random_genome_wrapper(n[i],L[i],m[i],k)
				before = time.time()
				rec_count[f,i,j], T_n[f,i,j], weight[f,i,j] = individual_test(algorithms[f],genome,motifs,k)
				after  = time.time()
				elapsed[f,i,j] = after-before
				out[outcount,:] = [ f , L[i] , m[j] , rec_count[f,i,j], T_n[f,i,j], weight[f,i,j] , elapsed[f,i,j] ]
				outcount+=1
 # Save 
import pickle
bench = [rec_count,T_n,weight,elapsed]
pickle.dump( bench, open( "benchmark.p", "wb" ) )
pickle.dump( out, open("benchmark_out.p","wb"))
np.savetxt("benchmark.csv", out, delimiter=",")


if(d2):
	fixL = 5
	fixm = 6
	plt.figure(1)
	plt.title('Count, Score and Runtime Analysis')
	# recursive counts vs L
	mat = rec_count
	# red dashes, blue dashes and green dashes
	plt.subplot(421)
	plt.plot(L,mat[0,:,fixm], 'r--', L, mat[1,:,fixm],'b--')
	plt.ylabel('Recursive Counts')
	#plt.ylim(-100,100000)
	# recursive counts vs m
	# red dashes, blue dashes and green dashes
	plt.subplot(422)
	plt.plot(m,mat[0,fixL,:], 'r--', m, mat[1,fixL,:],'b--')
	plt.ylabel('Recursive Counts')
	# T_n vs L
	mat = T_n
	# red dashes, blue dashes and green dashes
	plt.subplot(423)
	plt.plot(L,mat[0,:,fixm], 'r--', L, mat[1,:,fixm],'b--')
	plt.ylabel('T(n)')
	#plt.ylim(-100,6000000)
	# recursive counts vs m
	# red dashes, blue dashes and green dashes
	plt.subplot(424)
	plt.plot(m,mat[0,fixL,:], 'r--', m, mat[1,fixL,:],'b--')
	plt.ylabel('T(n)')
	# weights vs L
	mat = weight
	# red dashes, blue dashes and green dashes
	plt.subplot(425)
	plt.plot(L,mat[0,:,fixm], 'r--', L, mat[1,:,fixm],'b--')
	plt.ylabel('Average Score')
	# recursive counts vs m
	# red dashes, blue dashes and green dashes
	plt.subplot(426)
	plt.plot(m,mat[0,fixL,:], 'r--', m, mat[1,fixL,:],'b--')
	plt.ylabel('Average Score')
	# elapsed vs L
	mat = elapsed
	# red dashes, blue dashes and green dashes
	plt.subplot(427)
	plt.plot(L,mat[0,:,fixm], 'r--', L, mat[1,:,fixm],'b--')
	plt.xlabel('Length of Genome (L) '+'(fixed m='+str(m[fixm])+')')
	plt.ylabel('Time Elapsed (sec)')
	plt.ylim(-5,np.max(mat[1,fixL,:])+np.max(mat[1,fixL,:])*.1)
	# recursive counts vs m
	# red dashes, blue dashes and green dashes
	plt.subplot(428)
	plt.plot(m,mat[0,fixL,:], 'r--', m, mat[1,fixL,:],'b--')
	plt.xlabel('Number of Forward CTCF Binding Motifs (m)'+'(fixed L='+str(L[fixL])+')')
	plt.ylabel('Time Elapsed (sec)')
	plt.ylim(-5,np.max(mat[1,fixL,:])+np.max(mat[1,fixL,:])*.1)
	plt.show()
else:
	x = L
	y = m
	fig = plt.figure(1)
	plt.title('Count, Score and Runtime Analysis')
	# recursive counts vs L vs m, dp
	mat = rec_count
	Z = mat[1,:,:]
	sub  = '321'
	xlab = 'Length of Genome'
	ylab = 'Number of Forward CTCF Binding Motifs'
	zlab = 'Recursive Counts'
	title= 'Dynamic Programming Approach: Recursive Counts'
	plotin3d(fig,sub,x,y,Z,xlab,ylab,title)
	#  counts vs L vs m, dp
	mat = T_n
	Z = mat[1,:,:]
	sub  = '322'
	xlab = 'Length of Genome'
	ylab = 'Number of Forward CTCF Binding Motifs'
	zlab = 'Counts, T(n)'
	title= 'Dynamic Programming Approach: Counts, T(n)'
	plotin3d(fig,sub,x,y,Z,xlab,ylab,title)
	#  score vs L vs m, dp
	mat = weight
	Z = mat[1,:,:]
	sub  = '323'
	xlab = 'Length of Genome'
	ylab = 'Number of Forward CTCF Binding Motifs'
	zlab = 'Score'
	title= 'Dynamic Programming Approach: Score'
	plotin3d(fig,sub,x,y,Z,xlab,ylab,title)
	#  score vs L vs m, greedy
	Z = mat[0,:,:]
	sub  = '324'
	xlab = 'Length of Genome'
	ylab = 'Number of Forward CTCF Binding Motifs'
	zlab = 'Score'
	title= 'Greedy Approach: Score'
	plotin3d(fig,sub,x,y,Z,xlab,ylab,title)
	#  time vs L vs m, dp
	mat = elapsed
	Z = mat[1,:,:]
	sub  = '325'
	xlab = 'Length of Genome'
	ylab = 'Number of Forward CTCF Binding Motifs'
	zlab = 'Time Elapsed (sec)'
	title= 'Dynamic Programming Approach: Time Elapsed'
	plotin3d(fig,sub,x,y,Z,xlab,ylab,title)
	#  time vs L vs m, greedy
	Z = mat[0,:,:]
	sub  = '326'
	xlab = 'Length of Genome'
	ylab = 'Number of Forward CTCF Binding Motifs'
	zlab = 'Time Elapsed (sec)'
	title= 'Greedy Approach: Time Elapsed'
	plotin3d(fig,sub,x,y,Z,xlab,ylab,title)


	plt.show()


