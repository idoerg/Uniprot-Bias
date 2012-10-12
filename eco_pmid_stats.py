#!/usr/bin/env python
import sys
from collections import defaultdict
import networkx as nx
import matplotlib.pyplot as plt
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

def graph_pmid_eco(pmid_eco,eco_id_count):
        
    # Draw a graph (should I weight edges?)
    G = nx.Graph()
    for eco_list in pmid_eco.values():
        eco_id_list = [i[0] for i in eco_list]
        for eco_pair in combinations(eco_id_list,2):
            G.add_edge(*eco_pair)
    return G

def get_n_neighbors(G, pmid_eco, eco_id_meaning):
    # Number of neighbors per node. Also weighted by the number of papers a node (assertion
    # code) appears in.

    papers_per_term = term_count_in_papers(pmid_eco)
    foo = open("papers_per_term_2-OCT-2012.csv","w")
    for eco_id in G.nodes():
        cur_name = eco_id_meaning[eco_id]
        n_papers = float(papers_per_term[eco_id])
        num_neighbors = len(G.neighbors(eco_id))
        print >> foo, "%s\t%s\t%d\t%d\t%.2f" % (
                       eco_id, cur_name, num_neighbors, n_papers, num_neighbors/n_papers)
    foo.close() 

def term_count_in_papers(pmid_eco):
    papers_per_term = defaultdict(int)
    for pmid in pmid_eco:
        for term in pmid_eco[pmid]:
            term_id, term_name = term
            papers_per_term[term_id] += 1
    return papers_per_term
         
        

if __name__ == '__main__':
    eco_id_count, eco_name_count, eco_id_meaning, pmid_eco = \
        parse_eco_pmid(open(sys.argv[1]))

    
    G = graph_pmid_eco(pmid_eco,eco_id_count)
    

    get_n_neighbors(G, pmid_eco, eco_id_meaning)


            
            
        
    
    
