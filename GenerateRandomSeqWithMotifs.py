__author__ = 'sjroth'

import random

#Create a function to generate a random sequence of length L.
#Input: L
#Output: randomly generated data of length L
def generate_genome(L):
    i = 0
    genome = ''
    while i < L:
        genome += random.choice(['A','T','C','G'])
        i += 1
    return genome

#Create a function that checks if the coordinates for two overlaps.
#Input: i, j
#Output: boolean indicating whether the coordinates overlap.
def is_overlap(i,j):
    overlap = False
    if i[0] >= j[0] and i[0] <= j[1]: overlap = True
    if i[1] >= j[0] and i[1] <= j[1]: overlap = True
    return overlap

#Create a function to insert n motifs into a random sequence of length L.
#Input: n, L, motifs
#Output: randomly generated sequence with n motifs
def generate_genome_with_motifs(n,L,motifs):

    #Generate random genome.
    genome = generate_genome(L)

    #Insert n motifs. Make sure to insert the motifs in distinct locations.
    motif_locations = []
    for i in range(n):
        motif = random.choice(motifs)
        rand_location = random.randrange(L-4)
        motif_location = (rand_location,rand_location+4)
        not_used = True
        j = 0
        while j < len(motif_locations) and not_used:
            used_location = motif_locations[j]
            if is_overlap(used_location,motif_location): not_used = False
            j += 1
        motif_locations.append(motif_location)
        if not_used:
            genome = genome[:motif_location[0]]+motif+genome[motif_location[1]:]
    print motif_locations
    return genome
