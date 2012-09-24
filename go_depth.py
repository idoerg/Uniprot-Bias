#!/usr/bin/env python
import GO.go_utils as gu
import getopt
import sys

opts, args = getopt.getopt(sys.argv[1:],'i:f:o:',
                           ['infile=','field=','outfile='])

gocon, goc = gu.open_go(user="idoerg", passwd="mingus", db="MyGO")
infile = None
go_acc_field = 3 # starts at 1, not 0! 
outfile = sys.stdout
for o, a in opts:
    if o in ('-i','--infile'):
        infile = a
    elif o in ('-o','--outfile'):
        outfile = open(a,"w")
    elif o in ('-f','--field'):
        go_acc_field = int(a)
print "infile", infile
if infile:
    for inline in file(infile):
        sp_id = inline.strip().split()[3]
        go_acc = inline.strip().split()[go_acc_field-1]
        go_level = gu.go_level(go_acc,goc)
        try:
            go_term_type = gu.go_acc_to_term_type(go_acc, goc)
        except IndexError:
            outfile.write("%s\t%s\t%.1f\t%s\n" % (go_acc, "NOTERMTYPE",
                                                  go_level, sp_id))
            continue
        outfile.write("%s\t%s\t%.1f\t%s\n" % (go_acc, go_term_type, go_level, sp_id))
    
else:
    go_acc = args[0]
    go_level = gu.go_level(go_acc,goc)
    go_term_type = gu.go_acc_to_term_type(go_acc, goc)
    outfile.write("%s\t%s\t%.1f\n" % (go_acc, go_term_type, go_level))
gocon.close()
if outfile != sys.stdout:
    outfile.close()
        






