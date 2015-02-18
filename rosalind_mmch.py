from sys import argv
from collections import defaultdict
from rosalind_lia import ncr, factorial
from rosalind_utils import parse_file_fasta

def matching_count(min, max):
    return int(ncr(max, min)) * int(factorial(min))

def maximum_matchings(dna_str):
    base_counts = defaultdict(int)
    for base in dna_str:
        base_counts[base] += 1
    
    au_counts = [base_counts['A'], base_counts['U']]
    gc_counts = [base_counts['C'], base_counts['G']]
    
    max_au = max(au_counts)
    min_au = min(au_counts)
    
    max_gc = max(gc_counts)
    min_gc = min(gc_counts)
    
    return matching_count(min_au, max_au) * matching_count(min_gc, max_gc)
    
if __name__ == '__main__':
    fasta = parse_file_fasta(argv[1])
    fasta = fasta[0]
    dna_str = fasta.get_sequence()
    print maximum_matchings(dna_str)