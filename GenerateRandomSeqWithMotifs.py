__author__ = 'sjroth'

import random

#Create a function to generate a random sequence of length L.
#Input: L
#Output: randomly generated data of length L
def generate_genome(L):
    genome = ''.join([random.choice(['A','T','C','G']) for i in range(L)])
    return genome
'''
Create a function that generates the reverse complement of a given sequence.
Input: sequence
Output: reverse complement
'''
def rev_comp(sequence):
    out = ''
    comp_dict = {'A':'T','T':'A','C':'G','G':'C'}
    for char in sequence:
        out += comp_dict[char]
    return out[::-1]

#Create a function that checks if the coordinates for two overlaps.
#Input: i, j
#Output: boolean indicating whether the coordinates overlap.
def is_overlap(i,j):
    overlap = False
    if i[0] >= j[0] and i[0] <= j[1]: overlap = True
    if i[1] >= j[0] and i[1] <= j[1]: overlap = True
    return overlap

'''
Create a function that finds a new location for inserting a motif.
WARNING: Be sure that you are inserting a motif where there is room. Otherwise,
this function will loop forever.
Input: motif_locations, genome_length, motif_length
'''
def find_insertion_location(motif_locations,genome_length,motif_length):

    #Keep iterating until you find a new location.
    j = 0
    not_found_location = True
    while not_found_location:

        #Pick an initial random location.
        rand_location = random.randrange(genome_length-motif_length)
        motif_location = (rand_location,rand_location+motif_length)

        #Iterate through the locations that have been used.
        used = False
        i = 0
        while i < len(motif_locations) and not used:
            used_location = motif_locations[i]
            if is_overlap(used_location,motif_location): used = True
            i += 1
        if not used:
            not_found_location = False

    return motif_location

#Create a function to insert n motifs into a random sequence of length L.
#Input: n, L, motifs
#Output: randomly generated sequence with n motifs
def generate_genome_with_motifs(n,L,motifs):

    #Generate random genome.
    genome = generate_genome(L)

    #Motif length.
    motif_len = len(motifs[0])

    #Insert n motifs. Make sure to insert the motifs in distinct locations.
    #Keep track of which motifs are in the genome.
    motif_locations = []
    print genome
    unused_motifs = set(motifs)
    for i in range(n):

        #If there is an unused motif, pick it. Otherwise randomly choose a motif.
        if unused_motifs:
            motif = unused_motifs.pop()
        else:
            motif = random.choice(motif)
        rev_comp_motif = rev_comp(motif)

        #Find locations for the insertions.
        motif_location = find_insertion_location(motif_locations,L,motif_len)
        motif_locations.append(motif_location)
        rev_comp_location = find_insertion_location(motif_locations,L,motif_len)
        motif_locations.append(rev_comp_location)
        #Insert the motif into the genome.
        genome = genome[:motif_location[0]]+motif+genome[motif_location[1]:]
        genome = genome[:rev_comp_location[0]]+rev_comp_motif+\
                 genome[rev_comp_location[1]:]

    return genome

'''
Create a function that generates m motifs of length k.
Input: m, k
Output: list of motifs
'''
def generate_motifs(m,k):
    motifs = []
    for i in range(m):
        motif = generate_genome(k)
        while motif in motifs:
            motif = generate_genome(k)
        motifs.append(motif)
    return motifs

'''
Create a function that generates a genome of length L with n motifs of length k
inserted. There are m motifs that are randomly generated.
Input: L, n, m, k
Output: tuple containing the genome and motifs
'''
def random_genome_wrapper(n,L,m,k):

    #Generate the random motifs first.
    motifs = generate_motifs(m,k)
    genome = generate_genome_with_motifs(n,L,motifs)

    return genome,motifs