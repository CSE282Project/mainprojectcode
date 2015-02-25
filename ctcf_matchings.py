from sys import argv

class Matching:
    '''
    A matching object is a collection of vertex-disjoint edges
    '''
    def __init__(self, edges = None):
        if edges == None:
            edges = []
        self.vertices = set()
        self.edges = []
        self.add_edges(edges)
    
    def get_edges(self):
        return self.edges
    
    def add_edge(self, edge):
        u, v = edge
        assert u > 0
        assert v < 0
        
        assert u not in self.vertices
        assert v not in self.vertices
        
        self.vertices.update(edge)
        self.edges.append(edge)
    
    def add_edges(self, edges):
        for edge in edges:
            self.add_edge(edge)
        
    def add_matching(self, matching):
        '''
        merge with other matching
        '''
        if matching != None:
            self.add_edges(matching.get_edges())
    
    def get_weight(self):
        '''
        weight a matching as simply the average distance between vertices
        '''
        if len(self.edges) == 0:
            '''
            if matching is empty just return 0, 0
            this makes it easy when combining the scores of multiple matchings
            admittedly this may not be 100% robust but for now i think it's fine
            '''
            return 0, 0
        total = 0.0
        n = len(self.edges)
        for u, v in self.edges:
            total += abs(u - abs(v))
        return total / n, n
        
    def __add__(self, other):
        m = Matching()
        m.add_matching(self)
        m.add_matching(other)
        return m
        
    def __iter__(self):
        for edge in self.edges:
            yield edge
            
    def __repr__(self):
        return str(self.edges)

complements = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G'}

def reverse_complement(dna):
    revc = []
    for i in range(len(dna) - 1, -1, -1):
        base = dna[i]
        revc.append(complements[base])
    return ''.join(revc)
    
def get_vertices(genome, motifs, k):
    reverse_complements = map(reverse_complement, motifs)
    motif_nodes = []
    revc_nodes = []
    
    for i in range(len(genome) - k + 1):
        '''
        use 1-based indexing, so that positive/negative values delineate motifs from
        reverse complements
        '''
        pos = i + 1
        kmer = genome[i:i + k]
        if kmer in motifs:
            motif_nodes.append(pos)
        elif kmer in reverse_complements:
            revc_nodes.append(-pos)
    
    '''
    we know the graph is bipartite and semi-complete
    as such there is no need to create an adjacency list and matrix, we can infer the
    edge relation by simply bipartitioning the vertex set
    '''
    return motif_nodes, revc_nodes

def matching_helper(motif_nodes, revc_nodes, k, start, end, sub_matchings):
    '''
    helper function to find the minimum-weight maximal on a certain subsection of the 
    genome
    
    Args:
        motif_nodes: the locations of CTCF binding motifs
        
        revc_nodes: the locations of reverse complements to CTCF binding motifs
        
        k: length of the motifs
        
        start, end: the endpoints of the region on which we should find our matching
        ensures that whatever matching the callee returns does not cross with the caller's
        matching
        
        sub_matchings: dictionary of minimally-weighted matchings on previously computed
        intervals, to avoid repeated computation
        
    Returns:
        The non-crossing maximal matching that minimizes the average distance between 
        nodes in the matching among the space of all non-crossing maximal matchings on 
        the interval.  Returns an empty matching if there are no edges within this
        interval
    '''
    
    if (start, end) in sub_matchings:
        return sub_matchings[(start, end)]
    
    
    inbounds = lambda x : start <= x and x <= end
    motifs = filter(inbounds, motif_nodes)
    reverses = filter(lambda x : inbounds(-x), revc_nodes)
    
    '''
    instead of initializing best_matching to None we initialize it to an empty matching
    if there are no edges within the interval, then an empty matching will be returned
    which is what we want
    '''
#    best_matching = Matching()
    best_matching = []
    best_score = float("inf")
    
    for motif in motifs:
        assert motif > 0
        for revc in reverses:
            assert revc < 0
            
            dist = abs(motif - abs(revc))
            if dist < k:
                continue
            '''
            adding edge (motif, revc) to the matching, then finding the best 
            non-crossing matching that includes this edge
            
            by (greedily) adding whatever non-crossing edges possible we assure that our
            matching is maximally non-crossing.  We take the best-weighted matching 
            among maximal non-crossing matchings.  We use memoization to avoid repated
            computation.
            '''
            
            first, last = sorted([motif, abs(revc)])
            
            # get best matching to the left of the edge, but within the interval
            left = matching_helper(motif_nodes, revc_nodes, k, start, first - k, sub_matchings)
            
            # get best matching to the right of the edge, 
            right = matching_helper(motif_nodes, revc_nodes, k, last + k, end, sub_matchings)
            
            # get best matching within the loop formed by this edge
            mid = matching_helper(motif_nodes, revc_nodes, k, first + k, last - k, sub_matchings)
            
            left_mean, left_n = left.get_weight()
            right_mean, right_n = left.get_weight()
            mid_mean, mid_n = left.get_weight()
            
            # get the weighted sum of the optimal sub-matchings that do not cross 
            # (motif, revc)
            total = (left_mean * left_n) + (right_mean * right_n) + (mid_mean * mid_n)
            
            total += dist
            n = left_n + right_n + mid_n + 1
            score = float(total) / n
            
            # if this matching is optimal on this interval, then combine all the optimal
            # sub-matchings and add this edge
            if score == best_score:
            #if score > best_score - espsilon and score < best_score + epsilon
#                best_matching = left + right + mid
#                best_matching.add_edge((motif, revc))
                best_matching_i = Matching()
                best_matching_i = left + right + mid
                best_matching_i.add_edge((motif, revc))
                best_matching.append(best_matching_i)
#                best_score = score
            elif score < best_score:
                best_matching = []
                best_matching_i = Matching()
                best_matching_i = left + right + mid
                best_matching_i.add_edge((motif, revc))
                best_matching.append(best_matching_i)
                best_score = score

    sub_matchings[(start, end)] = best_matching
    return best_matching
    
def maximal_matching(genome, motifs, k):
    # generate the graph
    motif_nodes, revc_nodes = get_vertices(genome, motifs, k)
    '''
    call the helper function, starting with an empty dictionary and setting the 
    endpoints to be the entire genome
    1-based indexing of nodes, so initial endpoints are (1, |G| - k)
    have not computed anything, so sub_matchings is empty, but since dictionaries are
    pointers we can pass by reference, and when it is updated by the callee during
    recursive calls the changes will be reflected in the caller's stack frame
    '''
    return matching_helper(motif_nodes, revc_nodes, k, 1, len(genome) - k, {})
    
if __name__ == '__main__':
    genome = 'TTGAACTGAGCGAAGAAAGATCGCAGACTGTCAA'
    motifs = ['TTGA', 'GCGA']  # revc = ['TCAA', 'TCGC']
    k = 4
    matching = maximal_matching(genome, motifs, k)
    print matching
    print matching.get_weight()
