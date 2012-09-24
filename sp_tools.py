
from Bio import SwissProt as SP
from Bio import Entrez, Medline
import matplotlib.pyplot as pyplot
from GO import go_utils as gu
import MySQLdb
from collections import *
from SQL.mysql_utils import mysql_query, mysql_query_noreturn


# edits: AMS Jan 4, 2012

#edit IF 2/17/2012
#edit IF 2-APR-2012 changed the efetch call 
# instead of id=idlist in you fetch line with id=",".join(idlist).

# Set of GO experimental evidence codes
EEC = set(['EXP','IDA','IPI','IMP','IGI','IEP'])

#Email address used for NCBI efetch tools
MY_EMAIL = 'idoerg@gmail.com'

#cPickle
#When using the pickeled data...
#    papers corresponds to any pickle that has 'papers' in the name
#    papers_prots corresponds to any pickle that has 'papers_prots' in the name
#    must load pickles first to use this code. e.g.
#        test1 = open('goa-pickles/goa_exp_papers.pik', 'rb')
#        data1 = cPickle.load(test1)

###########################################################
# mysqlConnect
###########################################################
def mysqlConnect():
    """This is where you set all the mysql login, db info"""
#    return MySQLdb.connect(host="mysql-dev.cgl.ucsf.edu", user="schnoes",
#                   passwd="schnoes", db="GeneOntology", port=13308)
    return MySQLdb.connect(host="localhost", user="idoerg",
                   passwd="mingus", db="MyGO")


def exp_in_papers(papers,papers_prots):
    # Annotations with experimental evidence codes only
    # exp_papers_prots: key: pubmed_id; value: protein count dictionary
    # A protein count dictionary: key: UniProtID; value: count of that protein.
    # The protein count is the number of different experimentally-evidenced GO codes this
    # protein receives for a given paper "p".
    #
    # exp_papers: key: pubmed_id; value: a list of "go_rec" records.
    # go_rec record is a dictionary. Keys (values): 'sp_id' (swissprot id); 
    # 'go_id': (GO ID); 'go_ec': (GO Evidence Code).
    # 
    # Can be used with SP & GOA data
    
    exp_papers = {}
    exp_papers_prots = {}
    for p in papers:
        for go_rec in papers[p]:
            if go_rec['go_ec'] in EEC:
                exp_papers.setdefault(p,[]).append(go_rec)
                if p not in exp_papers_prots:
                    exp_papers_prots[p] = {go_rec['sp_id']: 1}
                else:
                    exp_papers_prots[p][go_rec['sp_id']] = \
                    exp_papers_prots[p].get(go_rec['sp_id'],0)+1
                
    return exp_papers, exp_papers_prots

def not_exp_in_papers(papers,papers_prots):
    # Annotations with no experimental evidence codes
    # nexp_papers_prots: key: pubmed_id; value: protein count dictionary
    # A protein count dictionary: key: SwissProtID; value: count of that protein.
    # The protein count is the number of different experimentally-evidenced GO codes this
    # protein receives for a given paper "p".
    #
    # nexp_papers: key: pubmed_id; value: a list of "go_rec" records.
    # go_rec record is a dictionary. Keys (values): 'sp_id' (swissprot id); 
    # 'go_id': (GO ID); 'go_ec': (GO Evidence Code).
    #
    # Can be used with SP & GOA data
    
    nexp_papers = {}
    nexp_papers_prots = {}
    for p in papers:
        for go_rec in papers[p]:
            if go_rec['go_ec'] not in EEC:
                nexp_papers.setdefault(p,[]).append(go_rec)
                if p not in nexp_papers_prots:
                    nexp_papers_prots[p] = {go_rec['sp_id']: 1}
                else:
                    nexp_papers_prots[p][go_rec['sp_id']] = \
                    nexp_papers_prots[p].get(go_rec['sp_id'],0)+1
                
    return nexp_papers, nexp_papers_prots

def plot_papers_prots(papers_prots, exp_papers_prots, do_plot=True):
    # hist: key: number of proteins; value: number of papers describing that number of proteins
    #
    # Can be used with SP & GOA data
    
    hist  = {}
    exp_hist  = {}
    for p in papers_prots:
        n_prots = len(papers_prots[p].keys())
        hist[n_prots] = hist.get(n_prots,0) + 1
    for p in exp_papers_prots:
        n_prots = len(exp_papers_prots[p].keys())
        exp_hist[n_prots] = hist.get(n_prots,0) + 1
    # Make a list in  the format: [(n_papers, n_proteins), ...]
    hist_table = [(hist[i],i) for i in hist]
    hist_table.sort()
    hist_table.reverse()
    exp_hist_table = [(hist[i],i) for i in hist]
    exp_hist_table.sort()
    exp_hist_table.reverse()
    if do_plot:
        pyplot.figure()
        pyplot.loglog()
        pyplot.xlabel('papers')
        pyplot.ylabel('proteins')
        pyplot.plot([x[0] for x in hist_table],[y[1] for y in hist_table],'ob')
        pyplot.plot([x[0] for x in exp_hist_table],[y[1] for y in exp_hist_table],'xr')
        pyplot.show()
        
    return None
 
def ec_stats(papers):
    """Creat a dict that holds the types and number of uses of each evidence code associate 
    with that particular paper."""
    # FN Can be used with SP & GOA data
    
    ec_count = {}
    #ec_count[PMID] = { {Ev Code : # of times it was used} } 
    for p in papers:
        #Papers[PMID] = [a list of dicts: each dict containing 3 entries: swissProt ID entry (key='sp_id'), 
        # go_id entry (key = 'go_id'), GO evidence code (key = 'go_ec')]    
        ec_count[p] = {}
        for annot in papers[p]:
            ec = annot['go_ec']
            ec_count[p][ec] = ec_count[p].get(ec,0) + 1
    return ec_count

def prot_stats(papers):
    #
    # Can be used with SP & GOA data
    
    prot_count = {}
    for p in papers:
        prot_count[p] = {}
        for annot in papers[p]:
            prot = annot['sp_id']
            prot_count[p][prot] = prot_count[p].get(prot,0) + 1
    return prot_count

def top_go_terms(papers,outpath=None,top=20):
    """Determines the top GO terms annotated in the analysis set and 1) puts it in 
    the output list top_go and 2) writes it out to a tab delim file 'outpath'
    
    Note: this function is currently identical to top_ontology()"""    #
    # Can be used with SP & GOA data
    
    go_count = {}
    for p in papers:
        for rec in papers[p]:
            go_id = rec['go_id']
            go_count[go_id] = go_count.get(go_id,0) + 1
    # top_go = [(# of uses, GO ID)]
    top_go = [(i[1],i[0]) for i in go_count.items()]
    top_go.sort()
    if outpath:
        go_con = mysqlConnect()
        go_cur = go_con.cursor()
        f = open(outpath,"w")
        for i in top_go[-top:]:
            name = gu.go_acc_to_name(i[1],go_cur)
            f.write("%d\t%s\t%s\n" % (i[0], i[1], name))
        go_hist = {}
        for i in top_go:
            go_hist[i[0]] = go_hist.get(i[0],0) + 1
        go_hist_list = [(h[1],h[0]) for h in go_hist.items()] 
        go_hist_list.sort()
        fhist = open("hist_%s" % outpath, "w")
        for h in go_hist_list:
#            print h
            fhist.write("%d\t%d\n" % h)
        
        f.close()
        fhist.close()
        go_con.close()
    return top_go
        
def go_terms_per_paper(papers,outpath=None,top=20):
    #
    # Can be used with SP & GOA data
    
    go_count = {}
    
    go_con = mysqlConnect()
    go_cur = go_con.cursor()
    
    for p in papers:
        for rec in papers[p]:
            go_id = rec['go_id']
            try:
                name = gu.go_acc_to_name(go_id,go_cur)
            except IndexError: #sometimes the GO ID given is actually a synonym
                try:
                    name = gu.go_acc_to_synonym_name(go_id, go_cur)
                except IndexError: #sometimes it just doesn't work
                    print "problem with GO ID", go_id
                    name = ''
            gokey = (go_id, name)
            # go_count[PMID] = {{(GO ID, GO Term Text) : # times paper gives this annotaion}}
            if p in go_count:
                go_count[p][gokey] = go_count[p].get(gokey,0) + 1
            else:
                go_count[p] = {gokey: 1}
    go_con.close()
    return go_count

###########################################################
# go_terms_with_ec_per_paper
###########################################################
def go_terms_with_ec_per_paper(papers,outpath=None,top=20):
    """Create a dict that counts up how many times a specific (GO ID, GO Term Text, EvCode) 
    tuple occurs for each paper"""
    #
    # Can be used with SP & GOA data
    
    go_ec_count = {}
    
    go_con = mysqlConnect()
    go_cur = go_con.cursor()
    for p in papers:
        for rec in papers[p]:
            go_id = rec['go_id']
            go_ec = rec['go_ec']
            try:
                name = gu.go_acc_to_name(go_id,go_cur)
            except IndexError: #sometimes the GO ID given is actually a synonym
                try:
                    name = gu.go_acc_to_synonym_name(go_id, go_cur)
                except IndexError: #sometimes it just doesn't work
                    print "problem with GO ID", go_id
                    name = ''
            gokey = (go_id, name, go_ec)
            # go_ec_count[PMID] = {{(GO ID, GO Term Text, Ev Code) : # times paper gives this annotaion}}
            if p in go_ec_count:
                go_ec_count[p][gokey] = go_ec_count[p].get(gokey,0) + 1
            else:
                go_ec_count[p] = {gokey: 1}
    go_con.close()
    return go_ec_count


###########################################################
# go_terms_with_ec_per_paper
###########################################################
def print_paper_per_prots_go(papers_annots2_dict, all_tt_count, go_ec_count,
                             allEvCodes_dict, sortedProtsPerPaper_tuple, outpath, top=20):
    """Prints out all information that we have for each paper(PMID).
    papers_annots2_dict: dict of the top X papers, with title, year and journal name
    all_tt_count: this is a dict that gives how many term types each paper annotates
    go_ec_count: this is a dict that gives how many times a paper gives a specific (go ID, go Name, Ev code) annotation
    allEvCodes_dict: this is a dict that gives how many times a paper supports a given experimental ev Code (EEC global).
    sortedProtsPerPaper_tuple: this is a sorted named tuple (largest first) of all the papers, sorted by the number of proteins the paper annotates
    outpath: where the data will be printed out.
    top: the number of papers we want to print out. Sorted by the number of proteins the paper annotates.
    """
    outFile = open(outpath, 'w')
    #print out header line
    outFile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % 
                  ('Num Prots', 'Num Annots', 'PMID', 'Title', 'Year', 'Journal', 'MFO Annot', 'BPO Annot', 'CCO Annot', 
                   'Num EXP', 'Num IDA', 'Num IEP', 'Num IGI', 'Num IMP', 'Num IPI',   'GO ID', 'GO Name', 'Ev Code', 'Num Used'))
    out_list = []
    type_list = ['molecular_function', 'biological_process', 'cellular_component']
    for sortedPaper in sortedProtsPerPaper_tuple[0:top]:
        # add the number of proteins & the PMID
        out_list.append(str(sortedPaper.numProts)) #[0]
        out_list.append(str(sortedPaper.PMID)) # will be [2]
        # get paper info
        try:
            out_list.extend(papers_annots2_dict[sortedPaper.PMID][1:]) #add title, year & journal [3, 4, 5]
        except KeyError:
            print "problem with", sortedPaper.PMID
            continue
        out_list.insert(1, str(papers_annots2_dict[sortedPaper.PMID][0])) # add number of annotations, [1]
        # get GO type info, [6, 7, 8]
        for t in type_list:
            try:
                out_list.append(str(all_tt_count[sortedPaper.PMID][t]))
            except:
                out_list.append('0')
        #get ev code data [9, 10, 11, 12, 13, 14]
        for ec in sorted(EEC):
            try:
                out_list.append(str(allEvCodes_dict[(sortedPaper.PMID, ec)]))
            except:
                out_list.append('0')
        # get GO ID info
        goInfo_dict = go_ec_count[sortedPaper.PMID]
        for key, value in goInfo_dict.iteritems(): # [15, 16, 17, 18....]
            out_list.extend(list(key))
            out_list.append(str(value))
        outFile.write('\t'.join((out_list)))
        outFile.write('\n')
        out_list = []
        
    outFile.close()  
    return



###########################################################
# sort_papers_prots
###########################################################
def sort_papers_prots(papers_prots):
    """Sort the dictionary papers_prots according to the number of proteins annotated by 
    a particular paper (PMID). Return the sorted tuple (sorted highest to lowest) as a 
    named tuple."""
    ProtsPerPaper_collect = namedtuple('ProtsPerPaper_collect', 'numProts PMID')
    #creates a tuple of (numProts, PMID) sorted on numProts highest to lowest
    sortedProtsPerPaper_tuple = sorted([ProtsPerPaper_collect(len(value), key) for (key, value) in papers_prots.items()],
                                       reverse=True)
    return sortedProtsPerPaper_tuple

def top_ontology(papers,outpath=None,top=20):
    """Determines the top GO terms annotated in the analysis set and 1) puts it in 
    the output dict top_go and 2) writes it out to a tab delim file 'outpath'
    
    Note: this function is currently identical to top_go_terms()"""
    #
    # Can be used with SP & GOA data
    
    go_count = {}
    for p in papers:
        for rec in papers[p]:
            go_id = rec['go_id']
            go_count[go_id] = go_count.get(go_id,0) + 1
    top_go = [(i[1],i[0]) for i in go_count.items()]
    top_go.sort()
    if outpath:
        go_con = mysqlConnect()
        go_cur = go_con.cursor()
        f = open(outpath,"w")
        for i in top_go[-top:]:
            name = gu.go_acc_to_name(i[1],go_cur)
            f.write("%d\t%s\t%s\n" % (i[0], i[1], name))
        go_hist = {}
        for i in top_go:
            go_hist[i[0]] = go_hist.get(i[0],0) + 1
        go_hist_list = [(h[1],h[0]) for h in go_hist.items()] 
        go_hist_list.sort()
        fhist = open("hist_%s" % outpath, "w")
        for h in go_hist_list:
#            print h
            fhist.write("%d\t%d\n" % h)
        
        f.close()
        fhist.close()
        go_con.close()
    return top_go
 

def top_papers(papers,outpath=None,delim="\t", top=20):
    """This function fetches all the relevant PubMed info for each PMID in 'papers' and 
    1) puts it into a list and 2) outputs it to a file named in outpath."""
    #
    # Can be used with SP & GOA data
    
    papers_annots = [(len(papers[p]), p) for p in papers]
    papers_annots2 = []
        
    papers_annots.sort()
    idlist = [p[1] for p in papers_annots[-top:]]
    Entrez.email = "idoerg@gmail.com"
    h = Entrez.efetch(db="pubmed", id=",".join(idlist), 
                          rettype="medline", retmode="text")
    medrecs = list(Medline.parse(h))
    titles = [medrec.get("TI","?") for medrec in medrecs]
    years = [medrec.get("DP","?") for medrec in medrecs]
    journals = [medrec.get("JT", "?") for medrec in medrecs]
    for p, title, year, journal in zip(papers_annots[-top:], titles,years, journals):
        papers_annots2.append((p[0],p[1], title, year.split()[0].strip(), journal))
    if outpath:
        fout = open(outpath,"w")
        print >> fout, "num proteins\tpubmed ID\tTitle\tYear\tJournal"
        for p in papers_annots2:
            print >> fout, "%d\t%s\t%s\t%s\t%s" % p
        fout.close()
    #papers_annots2 = [(# all annotations, PMID, Title, Year, Journal)] 
    return papers_annots2

###########################################################
# top_papers_dict
###########################################################
def top_papers_dict(papers, outpath=None,delim="\t", top=None):
    """This function fetches all the relevent PubMed info for each PMID in 'papers' 
    (at the limit supplied in 'top') and 1) puts it into a dict."""
    #
    # Can be used with SP & GOA data
    
#    papers_annots = [(len(papers_prots[p]), p) for p in papers_prots]
    papers_annots = [(len(papers[p]), p) for p in papers]
    papers_annots2_dict = {}
        
    papers_annots.sort()
    if top is None:
        negTop = 0
    else:
        negTop = -top
    idlist = [p[1] for p in papers_annots[negTop:]]
    Entrez.email = MY_EMAIL
    h = Entrez.efetch(db="pubmed", id=",".join(idlist), 
                          rettype="medline", retmode="text")
    medrecs = list(Medline.parse(h))
    titles = [medrec.get("TI","?") for medrec in medrecs]
    years = [medrec.get("DP","?") for medrec in medrecs]
    journals = [medrec.get("JT", "?") for medrec in medrecs]
    for p, title, year, journal in zip(papers_annots[negTop:], titles,years, journals):
        #papers_annots2_dict[PMID] = [# of total annotations, Title, Year, Journal] 
        papers_annots2_dict[p[1]] = [len(papers[p[1]]), title, year.split()[0].strip(), journal]
    """if outpath:
        fout = open(outpath,"w")
        print >> fout, "num proteins\tpubmed ID\tTitle\tYear\tJournal"
        for p in papers_annots2:
            print >> fout, "%d\t%s\t%s\t%s\t%s" % p
        fout.close()
    """
    return papers_annots2_dict

def go_terms_hi_v_lo(papers,top=50):
    my_singleton_papers = {}
    my_very_lo_papers = {}
    my_lo_papers = {}
    my_top_papers = {}
    top_papers_list = top_papers_dict(papers,top=top)
    top_papers_pubids = top_papers_list.keys()
    for pubid in papers:
        if pubid in top_papers_pubids:
            my_top_papers[pubid] = papers[pubid]
        elif len(papers[pubid]) == 1:
            my_singleton_papers[pubid] = papers[pubid]
        elif len(papers[pubid]) <= 10:
            my_very_lo_papers[pubid] = papers[pubid]
        else:
            my_lo_papers[pubid] = papers[pubid]
    # very low: up to 10 annotations
    top_papers_go = open("top_papers_go.txt","w")
    lo_papers_go = open("lo_papers_go.txt","w")
    very_lo_papers_go = open("very_lo_papers_go.txt","w")
    singleton_papers_go = open("singleton_papers_go.txt","w")
    for hi_id in my_top_papers:
        for l in my_top_papers[hi_id]:
            print >> top_papers_go, "%s\t%s\t%s" % (hi_id,l['go_ec'], l['go_id'])
    for lo_id in my_lo_papers:
        for l in my_lo_papers[lo_id]:
            print >> lo_papers_go, "%s\t%s\t%s" % (lo_id,l['go_ec'], l['go_id'])
    for v_lo_id in my_very_lo_papers:
        for l in my_very_lo_papers[v_lo_id]:
            print >> very_lo_papers_go, "%s\t%s\t%s" % (v_lo_id,l['go_ec'], l['go_id'])
    for s_id in my_singleton_papers:
        for l in my_singleton_papers[s_id]:
            print >> singleton_papers_go, "%s\t%s\t%s" % (s_id,l['go_ec'], l['go_id'])
    top_papers_go.close()
    lo_papers_go.close()
    very_lo_papers_go.close()
    singleton_papers_go.close()

def countprot(pubid,paper_rec):
	prot_count = {}
	for annot_rec in paper_rec:
		prot_count[(pubid,annot_rec['sp_id'])] = None
	return len(prot_count)
		
def go_terms_hi_v_lo_withprots(papers,top=50):
    my_singleton_papers = {}
    my_very_lo_papers = {}
    my_lo_papers = {}
    my_hi_papers = {}
    my_top_papers = {}
    top_papers_list = top_papers_dict(papers,top=top)
    top_papers_pubids = top_papers_list.keys()
    for pubid in papers:
        protein_count = countprot(pubid,papers[pubid])
        if pubid in top_papers_pubids:
            my_top_papers[pubid] = papers[pubid]
        if protein_count == 1:
            my_singleton_papers[pubid] = papers[pubid]
        elif protein_count <= 10 and protein_count > 1:
            my_very_lo_papers[pubid] = papers[pubid]
        elif protein_count <= 100 and protein_count > 10:
            my_lo_papers[pubid] = papers[pubid]
        elif protein_count > 100:
            my_hi_papers[pubid] = papers[pubid]
	# singleton: 1 annotation
    # very low: 2 to 10 annotations
	# low: 11-100 annotations
	# high > 100 annotations
	# top: top 50
    top_papers_go = open("top_papers_go_pp.txt","w")
    hi_papers_go = open("hi_papers_go_pp.txt","w")
    lo_papers_go = open("lo_papers_go_pp.txt","w")
    very_lo_papers_go = open("very_lo_papers_go_pp.txt","w")
    singleton_papers_go = open("singleton_papers_go_pp.txt","w")
    for top_id in my_top_papers:
        for l in my_top_papers[top_id]:
            print >> top_papers_go, "%s\t%s\t%s\t%s" % (top_id,l['go_ec'],
                                                    l['go_id'], l['sp_id'])
    for hi_id in my_hi_papers:
        for l in my_hi_papers[hi_id]:
            print >> hi_papers_go, "%s\t%s\t%s\t%s" % (hi_id,l['go_ec'],
                                                    l['go_id'], l['sp_id'])

    for lo_id in my_lo_papers:
        for l in my_lo_papers[lo_id]:
            print >> lo_papers_go, "%s\t%s\t%s\t%s" % (lo_id,l['go_ec'],
                                                    l['go_id'], l['sp_id'])
    for v_lo_id in my_very_lo_papers:
        for l in my_very_lo_papers[v_lo_id]:
            print >> very_lo_papers_go, "%s\t%s\t%s\t%s" % (v_lo_id,l['go_ec'],
                                                    l['go_id'], l['sp_id'])
    for s_id in my_singleton_papers:
        for l in my_singleton_papers[s_id]:
            print >> singleton_papers_go, "%s\t%s\t%s\t%s" % (s_id,l['go_ec'],
                                                    l['go_id'], l['sp_id'])
    top_papers_go.close()
    hi_papers_go.close()
    lo_papers_go.close()
    very_lo_papers_go.close()
    singleton_papers_go.close()
from scipy import stats
import numpy as np

def go_info_content_stats_all(inpath_list):
    fignum=  0
    colors = ('b','g','r','c','m','o')
    titles = []
    tot_bp = []
    tot_mf = []
    tot_cc = []
    go_ic = {} # dictionary of information content
    # Deal with information content
    for inline in file(inpath_list[0]):
        go_id, ic, term_type = inline.split()
        go_ic[go_id] =  (float(ic), term_type)
    bp_ic_list, mf_ic_list, cc_ic_list = go_ic_stats(go_ic)
    for inpath in inpath_list[1:]:
        titles.append(inpath)
        mf_ic_list = []
        bp_ic_list = []
        cc_ic_list = []
        for inline in file(inpath):
            go_id, term_type, depth, prot_id = inline.split()
            if go_id not in go_ic:
                print "can't find ", go_id
                continue

            ic = go_ic[go_id][0]
            
            if term_type == 'molecular_function':
                mf_ic_list.append(ic)
            elif term_type == 'biological_process':
                bp_ic_list.append(ic)
            elif term_type == 'cellular_component':
                cc_ic_list.append(ic)
        tot_mf.append(mf_ic_list)
        tot_bp.append(bp_ic_list)
        tot_cc.append(cc_ic_list)
        
    uvec, pval_vec = ontology_vector_stats(tot_mf)
    print "MF Uvals", uvec
    print "MF pvals", pval_vec
    uvec, pval_vec = ontology_vector_stats(tot_bp)
    print "BP Uvals", uvec
    print "BP pvals", pval_vec
    uvec, pval_vec = ontology_vector_stats(tot_cc)
    print "CC Uvals", uvec
    print "CC pvals", pval_vec
    locs = ('MF','BP','CC', '   ',
            'MF','BP','CC', '   ',
            'MF','BP','CC', '   ',
            'MF','BP','CC', '   ',
            'MF','BP','CC')
    locs2 = ('MF=1','<10','<100', '>100','   ',
             'BP=1','<10','<100', '>100','   ',
             'CC=1','<10','<100', '>100')
    pyplot.title('Information Content Statstics')
    distros_to_plot = (tot_mf, tot_bp, tot_cc)
    for i in range(3):
        bp = pyplot.boxplot(distros_to_plot[i],
                 positions=range(5*i+1,5*i+5))
        pyplot.setp(bp['boxes'],color=colors[i])
        pyplot.setp(bp['medians'],color='black')
        pyplot.setp(bp['fliers'],color=colors[i])
        pyplot.setp(bp['whiskers'],color=colors[i])
        pyplot.setp(bp['caps'],color=colors[i])

#    for i in range(len(tot_bp)):
#        bp=pyplot.boxplot( (tot_mf[i],tot_bp[i],tot_cc[i]), bootstrap=5000,
#                        positions=range(4*i+1,4*i+4),notch=False )
#        pyplot.setp(bp['boxes'],color=colors[i])
#        pyplot.setp(bp['medians'],color='black')
#        pyplot.setp(bp['fliers'],color=colors[i])
#        pyplot.setp(bp['whiskers'],color=colors[i])
#        pyplot.setp(bp['caps'],color=colors[i])
#                        #positions=pos+(np.ones(4)*delta_pos[i]) )
#    pyplot.xticks(range(1,20),locs)
    pyplot.xticks(range(1,16),locs2)
    pyplot.axis('tight')
    pyplot.show()

def ontology_vector_stats(tot_onto):
    # Check for significant differences between the distribution 
    # of IC or GO-depth between the singleton and very low, very low and low, low and high, etc.
    #
    # Use Mann-Whitney U.
    pval_vec = []
    uvec = []
    for i in range(1,len(tot_onto)):
        u, pval  = stats.mannwhitneyu(tot_onto[i-1], tot_onto[i]) 
#        u, pval  = stats.ttest_ind(tot_onto[i-1], tot_onto[i]) 
        pval_vec.append(pval * 2) # to get a two-sided pval
        uvec.append(u)
    return uvec, pval_vec
        

def go_depth_stats_all(inpath_list):
    fignum=  0
    colors = ('b','g','r','c','m','o')
    titles = []
    tot_bp = []
    tot_mf = []
    tot_cc = []
    go_ic = {} # dictionary of information content, if we want it.
    # Deal with information content
    for inpath in inpath_list:
        fignum += 1
#        print "mean MF IC", np.mean(mf_ic_list)
#        print inpath
#        print "*************"
        bp_list, mf_list, cc_list = go_depth_stats(inpath)
            
            
        titles.append(inpath)
        tot_bp.append(bp_list)
        tot_mf.append(mf_list)
        tot_cc.append(cc_list)
#       pyplot.figure(fignum)
#       pyplot.boxplot((bp_list, mf_list, cc_list))
    print "tot fracs"
    bp_sum = sum([len(i) for i in tot_bp])
    mf_sum = sum([len(i) for i in tot_mf])
    cc_sum = sum([len(i) for i in tot_cc])
    for i, inpath in enumerate(inpath_list):
        print inpath
        print "BP fraction", float(len(tot_bp[i]))/bp_sum
        print "MF fraction", float(len(tot_mf[i]))/mf_sum
        print "CC fraction", float(len(tot_cc[i]))/cc_sum
        print "*******************"
        
    uvec, pval_vec = ontology_vector_stats(tot_mf)
    print "MF Uvals", uvec
    print "MF pvals", pval_vec
    uvec, pval_vec = ontology_vector_stats(tot_bp)
    print "BP Uvals", uvec
    print "BP pvals", pval_vec
    uvec, pval_vec = ontology_vector_stats(tot_cc)
    print "CC Uvals", uvec
    print "CC pvals", pval_vec
    locs = ('MF','BP','CC', '   ',
            'MF','BP','CC', '   ',
            'MF','BP','CC', '   ',
            'MF','BP','CC', '   ',
            'MF','BP','CC')
    locs2 = ('MF=1','<10','<100', '>100','   ',
             'BP=1','<10','<100', '>100','   ',
             'CC=1','<10','<100', '>100')
    pyplot.title('GO Depth Statistics')
    distros_to_plot = (tot_mf, tot_bp, tot_cc)
    for i in range(3):
        bp = pyplot.boxplot(distros_to_plot[i],
                 positions=range(5*i+1,5*i+5))
        pyplot.setp(bp['boxes'],color=colors[i])
        pyplot.setp(bp['medians'],color='black')
        pyplot.setp(bp['fliers'],color=colors[i])
        pyplot.setp(bp['whiskers'],color=colors[i])
        pyplot.setp(bp['caps'],color=colors[i])
#    for i in range(len(tot_bp)):
##        pyplot.figure(i)
##        pyplot.xticks(range(4*i+1,4*i+4),
#        bp=pyplot.boxplot( (tot_mf[i],tot_bp[i],tot_cc[i]), bootstrap=1000,
#                        positions=range(4*i+1,4*i+4),notch=False )
##        pyplot.xticks(locs,('BPO%d' % i,'MFO%d' % i,'CCO%d' % i))
##        pyplot.xticks(locs)
#        pyplot.setp(bp['boxes'],color=colors[i])
#        pyplot.setp(bp['medians'],color='black')
#        pyplot.setp(bp['fliers'],color=colors[i])
#        pyplot.setp(bp['whiskers'],color=colors[i])
#        pyplot.setp(bp['caps'],color=colors[i])
#                        #positions=pos+(np.ones(4)*delta_pos[i]) )
    pyplot.xticks(range(1,16),locs2)
    pyplot.axis('tight')
    pyplot.show()
def go_ic_stats(go_ic):
    mf_ic_list = []
    bp_ic_list = []
    cc_ic_list = []
    for go_id in go_ic:
        ic, term_type = go_ic[go_id]
        if term_type == 'molecular_function':
            mf_ic_list.append(ic)
        elif term_type == 'biological_process':
            bp_ic_list.append(ic)
        elif term_type == 'cellular_component':
            cc_ic_list.append(ic)
    print "MF information content mean, stdv", np.mean(mf_ic_list), np.std(mf_ic_list)
    print "BP information content mean, stdv", np.mean(bp_ic_list), np.std(bp_ic_list)
    print "CC information content mean, stdv", np.mean(cc_ic_list), np.std(cc_ic_list)
    return mf_ic_list, bp_ic_list, cc_ic_list

def go_depth_stats(inpath):
    bp_list = []
    mf_list = []
    cc_list = []
    infile = file(inpath)
    for inline in file(inpath):
        inrec = inline.strip().split()
        if inrec[1] == 'biological_process':
            bp_list.append(float(inrec[2]))
        elif inrec[1] == 'molecular_function': 
            mf_list.append(float(inrec[2]))
        elif inrec[1] == 'cellular_component':
            cc_list.append(float(inrec[2]))
 
    n_sum = len(bp_list) + len(mf_list) + len(cc_list)

    print "mf n, frac", len(mf_list), len(mf_list)/float(n_sum)
    print "mf mean", np.mean(mf_list)
    print "mf stdv", np.std(mf_list)
        
    print "bp n, frac", len(bp_list), len(bp_list)/float(n_sum)
    print "bp mean", np.mean(bp_list)
    print "bp stdv", np.std(bp_list)

    print "cc n, frac", len(cc_list), len(cc_list)/float(n_sum)
    print "cc mean", np.mean(cc_list)
    print "cc stdv", np.std(cc_list)
        
    return bp_list, mf_list, cc_list 

###########################################################
# term_types_all_papers
###########################################################
def ev_codes_all_papers(papers):
    """Calculate the number of times a paper gives a certain experimental evidence code.
    Possible evidence codes: EEC global defined above.
    """
    allEvCodes_dict = {}
    for pmid, go_annot_list in papers.iteritems():
        for go_annot in go_annot_list:
            go_ec = go_annot['go_ec']
            if go_ec in EEC:
                #allEvCodes_dict[(PMID, evidence code)] = # times evidence code supported by that paper
                allEvCodes_dict[(pmid, go_ec)] = allEvCodes_dict.get((pmid, go_ec), 0) + 1
    return allEvCodes_dict


###########################################################
# term_types_all_papers
###########################################################
def term_types_all_papers(papers):
    """Count up how many times each paper annotates to a certain GO group:
    Possible term types:
    biological_process
    molecular_function
    cellular_component
    """
    all_tt_count = {}
    for pmid, annot_list in papers.iteritems():
        all_tt_count[pmid] = term_types_in_paper(annot_list)
    # all_tt_count[pmid] = { {term type : term type counts} }
    return all_tt_count



def term_types_in_paper(paper):
    """For each paper, count how often different term types appear.
    Possible term types:
    biological_process
    molecular_function
    cellular_component
    """
    # term types (ontologies) in a given paper
    tt_count = {}
    go_con = mysqlConnect()
    go_cursor = go_con.cursor()
    for prec in paper:
        go_id = prec['go_id']
        try:
            term_type = gu.go_acc_to_term_type(go_id, go_cursor)
        except IndexError:
            print go_id,"may be deprecated"
            continue
        # tt_count[term_type] = count of that term
        tt_count[term_type] = tt_count.get(term_type,0) + 1
    go_con.close()
    return tt_count # count of term types
        
    

get_term_type_sql = """
SELECT
    term_type
FROM 
    term
WHERE
    term.acc = "%s"
"""

def go_sum_in_papers(papers):
    # Frequency of GO terms in specific papers
    # Key: GO ID; Value: count of that GO ID.
    #
    # Can be used with SP & GOA data  
   
    go_count = {} # all GO terms
    exp_go_count = {} # GO terms with experimental evidence codes only
    for p in papers:
        for go_rec in papers[p]:
            # All evidence codes
            go_id = go_rec['go_id']
            go_count[go_id] = go_count.get(go_id,0) + 1
            # Experimental evidence codes only
            if go_rec['go_ec'] in EEC:
                exp_go_count[go_id] = exp_go_count.get(go_id,0) + 1
         
    
    return go_count, exp_go_count
            
def scalarp(item,scalar_example):
    return (type(item) == type(scalar_example) or
                item is None)

def flatten(sequence,scalarp,result=None):
    if result is None: result = []
    for item in sequence:
        if scalarp(item,''):
            result.append(item)
        else:
            flatten(item, scalarp, result)
    return result

acc_by_term_name_sql = \
"""
SELECT
    acc
FROM
    term 
WHERE
    name = '%s';
"""
find_ancestors_sql = """
    SELECT p.acc
    FROM
        graph_path
        INNER JOIN
        term AS t ON (t.id = graph_path.term2_id)
    INNER JOIN
        term AS p ON (p.id = graph_path.term1_id)
    WHERE t.acc = '%s'
    GROUP BY p.acc ORDER BY p.acc;
"""
def go_count_by_ontology(go_count):
    id_list = go_count.keys()
    go_count_mf = {}
    go_count_bp = {}
    go_count_cc = {}
    
    go_con = mysqlConnect()
    go_cursor = go_con.cursor()

    mf_root = mysql_query(acc_by_term_name_sql,("molecular_function",),go_cursor)
    cc_root = mysql_query(acc_by_term_name_sql,("cellular_component",),go_cursor)
    bp_root = mysql_query(acc_by_term_name_sql,("biological_process",),go_cursor)

    # split go_count dict into three dictionaries, mf, bp, cc
    obsolete_go = open("obsolete_go.txt","w")
    for go_id in go_count:
        try:
            term_type = flatten(
            mysql_query(get_term_type_sql,(go_id,),go_cursor),
                                        scalarp)[0]
        except IndexError:
            obsolete_go.write("%s\n" % go_id)
            continue
        if term_type == 'molecular_function':
            go_count_mf[go_id] = go_count_mf.get(go_id,0) + 1
        elif term_type == 'biological_process':
            go_count_bp[go_id] = go_count_bp.get(go_id,0) + 1
        elif term_type == 'cellular_component':
            go_count_cc[go_id] = go_count_cc.get(go_id,0) + 1
        else:
            obsolete_go.write("no term_type found for %s\n" % go_id)

        # Now, we propagate the count

        ancestor_list = flatten(mysql_query(find_ancestors_sql,
                               (go_id,),go_cursor),scalarp)
        for anc_go_id in ancestor_list:
            if anc_go_id != go_id:
                if term_type == 'molecular_function':
                    go_count_mf[anc_go_id] = go_count_mf.get(anc_go_id,0) + 1
                elif term_type == 'biological_process':
                    go_count_bp[anc_go_id] = go_count_bp.get(anc_go_id,0) + 1
                elif term_type == 'cellular_component':
                    go_count_cc[anc_go_id] = go_count_cc.get(anc_go_id,0) + 1
    go_con.close()
    obsolete_go.close()
    return go_count_mf, go_count_bp, go_count_cc

from math import log 
def go_info_content_by_papers(go_count,n_papers):
    # Information content as -log2(go_frequency)
    # Frequency is term_count/paper_count
    go_ic = {}
#    denom = float(go_count['all'])
    denom = float(n_papers)
    for go_id in go_count:
        go_ic[go_id] = -(log(go_count[go_id]/denom) / log(2.))
    return go_ic
def go_info_content_by_terms(go_count):
    # Information content as -log2(go_frequency)
    go_ic = {}
    denom = float(go_count['all'])
    for go_id in go_count:
        go_ic[go_id] = -(log(go_count[go_id]/denom) / log(2.))
    return go_ic
def go_info_content_all(go_count_mf, go_count_bp, go_count_cc):
    go_ic = {}
    go_ic_mf = go_info_content_by_terms(go_count_mf)
    go_ic_bp = go_info_content_by_terms(go_count_bp)
    go_ic_cc = go_info_content_by_terms(go_count_cc)
    # put in a general dictionary
    for i in go_ic_mf:
        go_ic[i] = go_ic_mf[i]
    for i in go_ic_bp:
        go_ic[i] = go_ic_bp[i]
    for i in go_ic_cc:
        go_ic[i] = go_ic_cc[i]
    outfile = open("papers_go_info_content.txt","w")
    # print it
    for i in go_ic_mf:
        outfile.write("%s\t%f\t%s\n" % (i,go_ic_mf[i], 'molecular_function'))
    for i in go_ic_bp:
        outfile.write("%s\t%f\t%s\n" % (i,go_ic_bp[i], 'biological_process'))
    for i in go_ic_cc:
        outfile.write("%s\t%f\t%s\n" % (i,go_ic_cc[i], 'cellular_component'))
    outfile.close()
     
def go_info_content_by_papers_all(go_count_mf, go_count_bp, go_count_cc,n_papers):
    go_ic = {}
    go_ic_mf = go_info_content(go_count_mf)
    go_ic_bp = go_info_content(go_count_bp)
    go_ic_cc = go_info_content(go_count_cc)
    # put in a general dicitonary
    for i in go_ic_mf:
        go_ic[i] = go_ic_mf[i]
    for i in go_ic_bp:
        go_ic[i] = go_ic_bp[i]
    for i in go_ic_cc:
        go_ic[i] = go_ic_cc[i]
    outfile = open("papers_go_info_content.txt","w")
    # print it
    for i in go_ic_mf:
        outfile.write("%s\t%f\t%s\n" % (i,go_ic_mf[i], 'molecular_function'))
    for i in go_ic_bp:
        outfile.write("%s\t%f\t%s\n" % (i,go_ic_bp[i], 'biological_process'))
    for i in go_ic_cc:
        outfile.write("%s\t%f\t%s\n" % (i,go_ic_cc[i], 'cellular_component'))
    outfile.close()
    
def go_in_papers(sp_path):
    # Returns: papers: key: pubmed_id; value: list of go_rec records
    # go_rec record is a dictionary. Keys (values): 'sp_id' (swissprot id); 
    # 'go_id': (GO ID); 'go_ec': (GO Evidence Code).
    
    # To be used with SP data, not GOA
    
    papers = {}
    go_ids = {}
    sp_recs = {}
    papers_prots = {}
    sph = open(sp_path)
    for sp_rec in SP.parse(sph):
        cur_go_recs = get_go_evidence_codes(sp_rec)
#        print cur_go_recs
        if not cur_go_recs: 
            continue
        cur_papers = get_papers(sp_rec)
        for paper in cur_papers:
            if paper not in papers_prots:
                papers_prots[paper] = {sp_rec.entry_name: 1}
            else:
                papers_prots[paper][sp_rec.entry_name] = \
                    papers_prots[paper].get(sp_rec.entry_name,0)+1
            for cur_go_rec in cur_go_recs:
                d1 = dict(sp_id=sp_rec.entry_name,
                          go_id=cur_go_rec[0],
                          go_ec=cur_go_rec[1])
                papers.setdefault(paper,[]).append(d1)
    return papers, papers_prots        
    
###########################################################
# printDict
###########################################################
def printDict(generic_dict, fileName):
    """Just does a quick print of generic_dict to file fileName)"""
    outFile = open(fileName, "w")

    for key, value in generic_dict.iteritems():
        outFile.write('Key: ' + str(key) + '\n')
        outFile.write("Value: " + str(value) + '\n')
    outFile.write('\n')
    outFile.close()


def go_in_papers_goa(goa_path):
    """Extract the GO data from the Uniprot GOA download"""
    papers = {}
    papers_prots = {}
    for inline in file(goa_path):
        if inline[0] == '!': continue
        db, db_object_id, db_object_symbol, qualifier, go_id, \
        db_reference, evidence, withit, aspect, \
        db_object_name, synonym, db_object_type, \
        taxon_id, date, assigned_by, \
        annotation_extension, gene_product_form_id = inline.rstrip('\n').split('\t')
        key_id = "%s:%s" % (db, db_object_id)

        if db_reference[:4] != "PMID": # only take the PMIDs, don't care about anything else
            continue
        paper = db_reference.split(':')[1] # Paper = the PMID
        if paper not in papers_prots:
            #papers_prots holds how many proteins each paper annotates
            #papers_prots[PMID] = {{key=Uniprot ID(key_id), value=# of times this paper annotates this protein}}
            #it is possible for one paper to produce more than one GO annotation for the same protien
            #Example: Unpript ID = Q9H1C4 and PMID = 19006693 
            papers_prots[paper] = {key_id: 1}
        else:
            papers_prots[paper][key_id] = \
                 papers_prots[paper].get(key_id,0)+1
    
        #Papers[PMID] = [a list of dicts: each dict containing 3 entries: swissProt ID entry (key='sp_id'), 
        # go_id entry (key = 'go_id'), GO evidence code (key = 'go_ec')]
        d1 = dict(sp_id=key_id,
                  go_id=go_id,
                  go_ec=evidence)
        papers.setdefault(paper,[]).append(d1)
        
        printDict(papers, 'papers.txt')
        printDict(papers_prots, 'papers_prots.txt')
    return papers, papers_prots
        

def get_go_evidence_codes(sp_rec):
    # isolate go evidence codes from sp_rec
    go_ids = []
    for xref in sp_rec.cross_references:
#        print xref
        if xref[0] == 'GO' and len(xref) == 4:
            go_ids.append((xref[1], xref[3][:3]))
    return go_ids

def get_papers(sp_rec):
    # get all pubs association with this sp_rec
    refs = {}
    for ref in sp_rec.references:
        for refid in ref.references:
            if refid[0] == 'PubMed':
                refs[refid[1]] = ref.title
    return refs
    
from sets import Set
def redundant_annotations(go_papers_dict):
    go_con, go_cur = gu.open_go()
    ancestors_found = {}
    to_remove = {}
    gpd_leaves_only = {}
    for pmid in go_papers_dict.keys(): #[:1000]:
        if len(go_papers_dict[pmid]) < 2:
            continue
        for i in range(len(go_papers_dict[pmid])):
            for j in range(len(go_papers_dict[pmid])):
                if j >= i: continue
                go_id_1 = go_papers_dict[pmid][i]['go_id']
                go_id_2 = go_papers_dict[pmid][j]['go_id']
                if go_id_1 == go_id_2: continue
                sp_id_1 = go_papers_dict[pmid][i]['sp_id']
                sp_id_2 = go_papers_dict[pmid][j]['sp_id']
                if sp_id_1 != sp_id_2: continue
                if gu.is_ancestor(go_id_1, go_id_2, go_cur):
                    ancestors_found[(pmid,go_id_1,go_id_2)] = \
                        ancestors_found.get((pmid,go_id_1,go_id_2),0) + 1
                    to_remove.setdefault(pmid,Set([])).add(j)
                    #ancestors_found.setdefault(pmid,[]).append(go_id_1,go_id_2)
                elif gu.is_ancestor(go_id_2, go_id_1, go_cur):
                    ancestors_found[(pmid,go_id_2,go_id_1)] = \
                        ancestors_found.get((pmid,go_id_2,go_id_1),0) + 1
                    to_remove.setdefault(pmid,Set([])).add(i)
                    #ancestors_found.setdefault(pmid,[]).append(go_id_2,go_id_1)
    go_con.close()
    for pmid in go_papers_dict:
        if pmid not in to_remove:
            gpd_leaves_only[pmid] = go_papers_dict[pmid]
            continue
        else:
            gpd_leaves_only[pmid] = []
            for i in range(len(go_papers_dict[pmid])):
                if i not in to_remove[pmid]:
                    gpd_leaves_only[pmid].append(go_papers_dict[pmid][i])
            
            

    return ancestors_found, to_remove,gpd_leaves_only

def correlate_goic_npapers():
    mf_depth  ={}
    bp_depth = {}
    cc_depth = {}

    mf_ic = {}
    bp_ic = {}
    cc_ic = {}
    all_ic = {}
    go_depth = {}
    papers2go = defaultdict(list)
    papers2nprots = {}

    # number of proteins per paper
    print "reading papers-nprots..."
    for inline in file("papers-nprots.txt"):
        pub_id, nprots = inline.split()
        papers2nprots[pub_id] = float(nprots)
    # information content of each go term
    print "reading papers_go_info_content..."
    for inline in file("papers_go_info_content.txt"):
        go_id, ic_s, onto_type = inline.split()
        ic = float(ic_s)
        all_ic[go_id] = ic
        if onto_type == "molecular_function":
            mf_ic[go_id] = ic
        elif onto_type == "biological_process":
            bp_ic[go_id] = ic
        elif onto_type == "cellular_component":
            cc_ic[go_id] = ic
    # depth for each go term
    print "reading all_papers_go_depth_pp..."
    for inline in file("all_papers_go_depth_pp.txt"):
        go_id, onto_type, d, prot = inline.split()
        depth = float(d)
        go_depth[go_id] = depth
        if onto_type == "molecular_function":
            mf_depth[go_id] = depth
        elif onto_type == "biological_process":
            bp_depth[go_id] = depth
        elif onto_type == "cellular_component":
            cc_depth[go_id] = depth
    print mf_depth.items()[:5],len(mf_depth)
    # go terms in each paper
    for inline in file("all_papers_go_pp_uniq.txt"):
        pub_id, go_id, prot_id = inline.split()
        papers2go[pub_id].append(go_id)

    # data points: go_depth vs. mean number of proteins in paper
    xyvec = []
    xvec = []
    yvec = []
    xydict = defaultdict(list)
    for pub_id in papers2nprots:
        x = papers2nprots[pub_id]
        go_depth_list = []
        for go_id in papers2go[pub_id]:
            if go_id in bp_depth:
#            if go_id in mf_depth:
                go_depth_list.append(go_depth[go_id])
#        y = np.mean(go_depth_list)
        for y in go_depth_list:
            xydict[x].append(y)

    for x in xydict:
        yvec = []
        for y in xydict[x]:
            yvec.append(y)
        xyvec.append( (x,np.mean(yvec)) )
#        xyvec.append( (x, xydict[x]) )
    xyvec.sort()

    # data points: go_ic vs. mean number of proteins in paper
    xyvec_ic = []
    xvec = []
    yvec = []
    xydict = defaultdict(list)
    for pub_id in papers2nprots:
        x = papers2nprots[pub_id]
        go_ic_list = []
        for go_id in papers2go[pub_id]:
            if go_id in bp_ic:
#            if go_id in mf_ic:
                go_ic_list.append(all_ic[go_id])
#        y = np.mean(go_depth_list)
        for y in go_ic_list:
            xydict[x].append(y)
    for x in xydict:
        yvec = []
        for y in xydict[x]:
            yvec.append(y)
        xyvec_ic.append( (x, np.mean(yvec)) )
    xyvec_ic.sort()
    return xyvec, xyvec_ic
