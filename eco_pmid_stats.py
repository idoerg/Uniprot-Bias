
from collections import defaultdict
def parse_eco_pmid(infile):
    eco_id_count = defaultdict(int)
    eco_name_count = defaultdict(int)
    eco_id_meaning = defaultdict(str)
    pmid_eco = defaultdict(list)
    for inline in infile:
        inrec = inline.strip().split('\t')
        pmid = inrec[0]
        print inrec
        for eco_rec in inrec[1:]:
    
            eco_name, eco_id = eco_rec.split('@')
            eco_id = eco_id.strip()
            pmid_eco[pmid].append((eco_id, eco_name))
            eco_id_count[eco_id] += 1
            eco_name_count[eco_name] += 1
            eco_id_meaning[eco_id] = eco_name
        
    return eco_id_count, eco_name_count, eco_id_meaning, pmid_eco

from itertools import combinations
def amat_pmid_eco(pmid_eco):
    adj_mat = defaultdict(int)
    for eco_list in pmid_eco.values():
        eco_id_list = [i[0] for i in eco_list]
        for eco_pair in combinations(eco_id_list,2):
            adj_mat[eco_pair] += 1

    return adj_mat

import networkx as nx
import matplotlib.pyplot as plt
def graph_pmid_eco(pmid_eco,eco_name_count):
    eco_name_numbering = {}
    eco_numbering_name = {}
    for i,eco_name in enumerate(eco_name_count):
        eco_name_numbering[eco_name] = i
        eco_numbering_name[i] = eco_name
        
    # Draw a graph (should I weight edges?)
    G = nx.Graph()
    for eco_list in pmid_eco.values():
        eco_name_list = [i[1] for i in eco_list]
        for eco_pair in combinations(eco_name_list,2):
            eco1 = eco_name_numbering[eco_pair[0]]
            eco2 = eco_name_numbering[eco_pair[1]]
            G.add_edge(eco1, eco2) 
    for i in range(len(eco_numbering_name)):
        print i, eco_numbering_name[i]
    nx.draw(G)
    plt.show()
    return G,eco_numbering_name, eco_name_numbering

def n_neighbors(G, eco_numbering_name):
    # Number of neighbors per node. Also weighted by the number of papers a node (assertion
    # code) appears in.

    foo = open("foo","w")
    bar = open("bar","w")
    for n in G.nodes():
        print >> foo, "%d\t%s\t%d" % (n, eco_numbering_name[n], len(G.neighbors(n)))
    foo.close() 

def term_count_in_papers(pmid_eco):
    papers_per_term = defaultdict(int)
    for pmid in pmid_eco:
        for term in pmid_eco[pmid]:
            term_id, term_name = term
            papers_per_term[term_name] += 1
    return papers_per_term
         
        

if __name__ == '__main__':
    eco_id_count, eco_name_count, eco_id_meaning, pmid_eco = \
        parse_eco_pmid(open(sys.argv[1]))

    
    




            
            
        
    
    
