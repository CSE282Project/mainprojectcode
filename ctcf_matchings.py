from sys import argv

class Matching:
    def __init__(self, edges = None):
        if edges == None:
            self.edges = []
        self.vertices = set()
        self.add_edges(edge)
    
    def get_edges(self):
        return self.edges
    
    def add_edge(self, edge):
        u, v = edge
        
        assert u not in self.vertices
        assert v not in self.vertices
        
        self.vertices.update(edge)
        self.edges.append(edge)
    
    def add_edges(self, edges):
        for edge in edges:
            self.add_edge(edge)
        
    def add_matching(self, matching):
        self.add_edges(matching.get_edges)
    
    def get_score(self):
        total = 0.0
        n = len(self.edges)
        for u, v in self.edges:
            total += abs(abs(u) - abs(v))
        return total / n, n
        
    def __add__(self, other):
        m = Matching()
        m.add_matching(self)
        m.add_matching(other)
        return m

complements = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G'}

def reverse_complements(dna):
    revc = []
    for base in range(len(dna) - 1, -1, -1):
        revc.append(complements[base])
    return ''.join(revc)
    
def get_vertices(genome, motifs, k):
    reverse_complements = map(reverse_complement, motifs)
    
    motif_nodes = []
    revc_nodes = []
    
    for i in range(len(genome) - k + 1):
        kmer = genome[i:i + k]
        if kmer in motifs:
            motif_nodes.append(i)
        elif kmer in revc_nodes:
            revc_nodes.append(-1)
    
    '''
    we know the graph is bipartite and semi-complete
    as such there is no need to create an adjacency list and matrix, we can infer the
    edge relation by simply bipartitioning the vertex set
    '''
    return motif_nodes, revc_nodes

def matching_helper(motif_nodes, revc_nodes, k, start, end, sub_matchings):
    if (start, end) in sub_matchings:
        return sub_matchings[(start, end)]
    
    inbounds = lambda x : start <= x and x <= end
    motifs = filter(inbounds, motif_nodes)
    reverses = filter(inbounds, revc_nodes)
    
    best_matching = None
    best_score = float("inf")
    
    for motif in motifs:
        for revc in revc_nodes:
            left = matching_helper(motif_nodes, revc_nodes, k, 0, start - 1, sub_matchings)
            right = matching_helper(motif_nodes, revc_nodes, k, end + 1, len(genome), sub_matchings)
            mid = matching_helper(motif_nodes, revc_nodes, k, start + k, end - 1, sub_matchings)
            
            left_mean, left_n = left.get_score()
            right_mean, right_n = left.get_score()
            mid_mean, mid_n = left.get_score()
            
            total = (left_mean * left_n) + (right_mean * right_n) + (mid_mean * mid_n)
            dist = abs(motif - abs(revc))
            total += dist
            n = left_n + right_n + mid_n + 1
            score = float(total) / n
            if score < best_score:
                best_matching = left + right + mid
                best_matching.add_edge((motif, revc))
                
    sub_matchings[(start, end)] = best_matching
    return best_matching
    
def maximal_matching(genome, motifs, k):
    motif_nodes, revc_nodes = get_vertices(genome, motifs, k)
    return matching_helper(motif_nodes, revc_nodes, k, 0, len(genome) - 1, {})