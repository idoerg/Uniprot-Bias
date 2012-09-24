#!/usr/bin/env python
import sys
from collections import defaultdict
from pylab import *
def parse_annot_file(annot_file):
    paper_sets_mfo = defaultdict(dict)
    paper_sets_bpo = defaultdict(dict)
    paper_sets_cco = defaultdict(dict)
    cluster_id = -1
    for inline in annot_file:
        if 'True' in inline or 'False' in inline: continue
        if inline.split()[0] == '#':
            # new cluster, First starts at 0
            cluster_id += 1
            paper_sets_mfo[cluster_id] = defaultdict(set)
            paper_sets_bpo[cluster_id] = defaultdict(set)
            paper_sets_cco[cluster_id] = defaultdict(set)
            continue
        (pub_id, prot_id, ontology, go_term, go_text) = inline.split('\t')
        if ontology == 'MFO':
            paper_sets_mfo[cluster_id][pub_id].add(go_term)
        elif ontology == 'BPO':
            paper_sets_bpo[cluster_id][pub_id].add(go_term)
        elif ontology == 'CCO':
            paper_sets_cco[cluster_id][pub_id].add(go_term)
        else:
            raise ValueError, "Ontology %s" % ontology
    return paper_sets_mfo, paper_sets_bpo, paper_sets_cco

def jaccard(insets):
    jdict = {}
    nclusters = len(insets)
    for cluster_id in range(nclusters):
        if len(insets[cluster_id]) < 2:
            continue
        intersect = insets[cluster_id].values()[0]
        union = set([])
        for pub_id in insets[cluster_id]:
            union |= insets[cluster_id][pub_id]
            intersect &= insets[cluster_id][pub_id]
        jdict[cluster_id] = float(len(intersect)) / float(len(union))
    return jdict
            
def jaccard_stats(jdict):
    if len(jdict) == 0:
        jmean = 0.0
        jstd = 0.0
        jsterr = 0.0
        jn = 0
    else:
        jmean = mean(jdict.values())
        jstd = std(jdict.values())
        jn = len(jdict)
        jsterr = jstd/sqrt(jn)
    return (jn, jmean, jstd, jsterr)
    
def multi_union(first, *others):
    return first.union(*others)

def multi_intersect(first, *others):
    return first.intersection(*others)

def consistency_one_organism(inpath):
    mfo_sets, bpo_sets, cco_sets = parse_annot_file(open(inpath))
    mfo_jdict = jaccard(mfo_sets)
    bpo_jdict = jaccard(bpo_sets)
    cco_jdict = jaccard(cco_sets)

    mfo_stats = jaccard_stats(mfo_jdict)
    bpo_stats = jaccard_stats(bpo_jdict)
    cco_stats = jaccard_stats(cco_jdict)

    return mfo_stats, bpo_stats, cco_stats

def all_annotation_consistencies(pathlist):
    outfile = open("foo.csv","w")

    outfile.write("file\tontology\tn\tmean\tstdv\tstderr\n")
    for inpath in pathlist:
        mfo_stats, bpo_stats, cco_stats = consistency_one_organism(inpath)
        outfile.write("%s\tmfo\t%d\t%.3f\t%.3f\t%.3f\n" % tuple([inpath] + list(mfo_stats)))
        outfile.write("%s\tbpo\t%d\t%.3f\t%.3f\t%.3f\n" % tuple([inpath] + list(bpo_stats)))
        outfile.write("%s\tcco\t%d\t%.3f\t%.3f\t%.3f\n" % tuple([inpath] + list(cco_stats)))
    outfile.close()
    
if __name__ == '__main__':
#    mfo_stats, bpo_stats, cco_stats = consistency_one_organism(sys.argv[1])
#    print "MFO n, mean, stdv, stderr", mfo_stats
#    print "BPO n, mean, stdv, stderr", bpo_stats
#    print "CCO n, mean, stdv, stderr", cco_stats
    all_annotation_consistencies(sys.argv[1:])
