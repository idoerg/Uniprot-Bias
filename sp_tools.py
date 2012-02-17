
from Bio import SwissProt as SP
from Bio import Entrez, Medline
import matplotlib.pyplot as pyplot
from GO import go_utils as gu
import MySQLdb
import collections


# edits: AMS Jan 4, 2012

#edit IF 2/17/2012
# Set of GO experimental evidence codes
EEC = set(['EXP','IDA','IPI','IMP','IGI','IEP'])

#Email address used for NCBI efetch tools
MY_EMAIL = 'schnoes@gmail.com'

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
    return MySQLdb.connect(host="mysql-dev.cgl.ucsf.edu", user="schnoes",
                   passwd="schnoes", db="GeneOntology", port=13308)


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
        out_list.extend(papers_annots2_dict[sortedPaper.PMID][1:]) #add title, year & journal [3, 4, 5]
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
    ProtsPerPaper_collect = collections.namedtuple('ProtsPerPaper_collect', 'numProts PMID')
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
    """This function fetches all the relevent PubMed info for each PMID in 'papers' and 
    1) puts it into a list and 2) outputs it to a file named in outpath."""
    #
    # Can be used with SP & GOA data
    
    papers_annots = [(len(papers[p]), p) for p in papers]
    papers_annots2 = []
        
    papers_annots.sort()
    idlist = [p[1] for p in papers_annots[-top:]]
    Entrez.email = MY_EMAIL
    h = Entrez.efetch(db="pubmed", id=idlist, 
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
def top_papers_dict(papers, papers_prots, outpath=None,delim="\t", top=None):
    """This function fetches all the relevent PubMed info for each PMID in 'papers' 
    (at the limit supplied in 'top') and 1) puts it into a dict."""
    #
    # Can be used with SP & GOA data
    
    papers_annots = [(len(papers_prots[p]), p) for p in papers_prots]
    papers_annots2_dict = {}
        
    papers_annots.sort()
    if top is None:
        negTop = 0
    else:
        negTop = -top
    idlist = [p[1] for p in papers_annots[negTop:]]
    Entrez.email = MY_EMAIL
    h = Entrez.efetch(db="pubmed", id=idlist, 
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
        term_type = gu.go_acc_to_term_type(go_id, go_cursor)
        # tt_count[term_type] = count of that term
        tt_count[term_type] = tt_count.get(term_type,0) + 1
    go_con.close()
    return tt_count # count of term types
        
    
def go_freq_in_papers(papers):
    # Frequency of GO terms in specific papers
    # Key: GO ID; Value: count of that GO ID.
    #
    # Can be used with SP & GOA data  
   
    go_ids = {} # all GO terms
    exp_go_ids = {} # GO terms with experimental evidence codes only
    for p in papers:
        for go_rec in papers[p]:
            # All evidence codes
            go_id = go_rec['go_id']
            go_ids[go_id] = go_ids.get(go_id,0) + 1
            # Experimental evidence codes only
            if go_rec['go_ec'] in EEC:
                exp_go_ids[go_id] = exp_go_ids.get(go_id,0) + 1
         
    return go_ids, exp_go_ids
            
    
    
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
    



