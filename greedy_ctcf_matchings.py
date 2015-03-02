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
            return None, 0
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

    def __contains__(self,item):
        u,v = item
        if v in self.vertices and u in self.vertices:  # avg: O(1), worst: O(n)
            for e in self.edges:
                if u == e[0] and v == e[1]:
                    return True
        return False


count = 0
def init_counter():
    '''
    Initializes the recursion counter
    '''
    global count
    count = 0

complements = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G'}

def reverse_complement(dna):
    global count
    revc = []
    for i in range(len(dna) - 1, -1, -1):
        count += 1
        base = dna[i]
        revc.append(complements[base])
    return ''.join(revc)
    
def get_vertices(genome, motifs, k):
    global count
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
        count += len(motifs) # O(n) to check if kmer in motifs
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

def merge_graph(motifs, revc):
    global count
    graph = []
    i = 0
    j = 0
    while i < len(motifs) and j < len(revc):
        count += 1
        elem1 = motifs[i]
        elem2 = revc[j]
        if elem1 < abs(elem2):
            graph.append(elem1)
            i += 1
        else:
            graph.append(elem2)
            j += 2
    graph += motifs[i:] + revc[j:]
    return graph
    
def greedy_matching(genome, motifs, k):
    global count
    motifs, revc = get_vertices(genome, motifs, k)
    graph = merge_graph(motifs, revc)
    print graph
    matching = Matching()
    up = False
    down = False
    prev = None
    for vertex in graph:
        count += 1
        if vertex > 0:
            if down:
                assert not up
                assert prev < 0
                matching.add_edge((vertex, prev))
                prev = None
                down = False
            elif not up:
                assert prev == None
                up = True
                prev = vertex
        elif vertex < 0:
            if up:
                assert not down
                assert prev > 0
                matching.add_edge((prev, vertex))
                prev = None
                up = False
            elif not down:
                assert prev == None
                down = True
                prev = vertex
                
    return matching

def parse_file(fname):
    f = open(fname)
    input = f.read().splitlines()
    genome = input[0]
    k = int(input[1])
    motifs = input[2:]
    return genome, k, motifs
    
if __name__ == '__main__':
    fname = argv[1]
    genome, k, motifs = parse_file(fname)
    matching = greedy_matching(genome, motifs, k)
    print matching
    print matching.get_weight()