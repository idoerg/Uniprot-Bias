#!/usr/bin/env python

from scipy import *
from pylab import *
from scipy import stats
import sp_tools as spt

onto_vec = ['biological_process', 'cellular_component',
            'molecular_function']

import sys
def go_depth_stats(infile):
    depth_vecs = {}
    for onto in onto_vec:
        depth_vecs[onto] = []
    for inline in infile:
        go_id, onto_type, depth = inline.strip().split()
        if onto_type not in depth_vecs:
            continue
        depth_vecs[onto_type].append(float(depth))
    print infile.name
    for i in depth_vecs:
        print i,mean(depth_vecs[i]), \
                median(depth_vecs[i]),std(depth_vecs[i])
                 
    return(depth_vecs)

def all_go_depth_stats(infiles):
    plot_data = {'biological_process':[], 'molecular_function':[], 
                 'cellular_component':[]}
    data1 = []
    data2 = []
    data3 = []
    all_depth_vecs = {}
    for infile in infiles:
        all_depth_vecs[infile] = go_depth_stats(open(infile))

    for infile in infiles:
        data1.append(all_depth_vecs[infile]['molecular_function'])
        data2.append(all_depth_vecs[infile]['biological_process'])
        data3.append(all_depth_vecs[infile]['cellular_component'])

    print "MFO ttest"
    print "singleton vs. very low", stats.ttest_ind(data1[0], data1[1])
    print "very low vs. low", stats.ttest_ind(data1[1], data1[2])
    print "low vs. top", stats.ttest_ind(data1[2], data1[3])
    print "BPO ttest"
    print "singleton vs. very low", stats.ttest_ind(data2[0], data2[1])
    print "very low vs. low", stats.ttest_ind(data2[1], data2[2])
    print "low vs. top", stats.ttest_ind(data2[2], data2[3])
    print "CCO ttest"
    print "singleton vs. very low", stats.ttest_ind(data3[0], data3[1])
    print "very low vs. low", stats.ttest_ind(data3[1], data3[2])
    print "low vs. top", stats.ttest_ind(data3[2], data3[3])
    figure()
    suptitle('molecular_function')
    xticks(arange(len(infiles)), infiles)
    boxplot(data1)

    figure()
    suptitle('biological_process')
    xticks(arange(len(infiles)), infiles)
    boxplot(data2)

    figure()
    suptitle('cellular_component')
    xticks(arange(len(infiles)), infiles)
    boxplot(data3)
    print mean(data3[0]), mean(data3[1]), mean(data3[2]), mean (data3[3])
    show()
        


if __name__ == '__main__':
    infiles = sys.argv[1:]
#    all_go_depth_stats(infiles)
    spt.go_depth_stats_all(infiles)

        

