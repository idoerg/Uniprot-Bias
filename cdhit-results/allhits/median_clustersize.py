#!/usr/bin/env python
import re
import sys
import numpy
def read_clstr(infile):
    cluster_rec = re.compile('^\>Cluster')
    nclust = 0
    cluster = {}
    for inrec in infile:
        if re.match(cluster_rec, inrec):
            nclust += 1
            cluster[nclust] = 0
            continue
        else:
            cluster[nclust] += 1
    return cluster

if __name__ == '__main__':
    cluster = read_clstr(open(sys.argv[1]))
    cluster_sizes = cluster.values()
    cluster_sizes_gt_1 = [i for i in cluster_sizes if i>1]
    print "STARTING %s" % (sys.argv[1]); sys.stdout.flush()
    print "*>=1*:n= %d; mean=%.2f median=%.2f" % (
           len(cluster_sizes), numpy.mean(cluster_sizes),
           numpy.median(cluster_sizes)) 
    sys.stdout.flush()
    print "*>1*: n= %d; mean=%.2f median=%.2f" % (
           len(cluster_sizes_gt_1),
           numpy.mean(cluster_sizes_gt_1),
           numpy.median(cluster_sizes_gt_1))
    sys.stdout.flush()
    print "*"*40; sys.stdout.flush()
