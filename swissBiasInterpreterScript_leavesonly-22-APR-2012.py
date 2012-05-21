import cPickle
import sp_tools

#Experimental
papersExp_handle = open('goa_exp_papers_lo.pik', 'rb')
papersExp_dict = cPickle.load(papersExp_handle)
papers_protsExp_handle = open('goa_exp_papers_prots_lo.pik', 'rb')
papers_protsExp_dict = cPickle.load(papers_protsExp_handle)

# going for the top fifty papers
top = 50
papers_annots2_dict = sp_tools.top_papers_dict(papersExp_dict, papers_protsExp_dict, top=top)
print "really long time... (1)"
all_tt_count = sp_tools.term_types_all_papers(papersExp_dict) #takes a really long time
print "...done"
go_ec_count = sp_tools.go_terms_with_ec_per_paper(papersExp_dict, top=top) # this takes a bit of time too
allEvCodes_dict = sp_tools.ev_codes_all_papers(papersExp_dict)
sortedProtsPerPaper_tuple = sp_tools.sort_papers_prots(papers_protsExp_dict)
sp_tools.print_paper_per_prots_go(papers_annots2_dict, all_tt_count, go_ec_count, allEvCodes_dict, 
                         sortedProtsPerPaper_tuple, "allExpPaperInfoTop50_lo.txt", top=top)

# Not experimental
"""papersNoExp_handle = open('goa-pickles/goa_not_exp_papers.pik', 'rb')
papersNoExp_dict = cPickle.load(papersNoExp_handle)
papers_protsNoExp_handle = open('goa-pickles/goa_not_exp_papers_prots.pik', 'rb')
papers_protsNoExp_dict = cPickle.load(papers_protsNoExp_handle)

# going for the top fifty papers
top = 50
papers_annots2_dict_no = sp_tools.top_papers_dict(papersNoExp_dict, papers_protsNoExp_dict, top=top)
all_tt_count_no = sp_tools.term_types_all_papers(papersNoExp_dict) #takes a really long time
go_ec_count_no = sp_tools.go_terms_with_ec_per_paper(papersNoExp_dict, top=top) # this takes a bit of time too
allEvCodes_dict_no = sp_tools.ev_codes_all_papers(papersNoExp_dict)
sortedProtsPerPaper_tuple_no = sp_tools.sort_papers_prots(papers_protsNoExp_dict)
sp_tools.print_paper_per_prots_go(papers_annots2_dict_no, all_tt_count_no, go_ec_count_no, allEvCodes_dict_no, 
                         sortedProtsPerPaper_tuple_no, "allNoExpPaperInfoTop50.txt", top=top)
"""
# 
# # All
# papersAllExp_handle = open('goa_papers.pik', 'rb')
# papersAllExp_dict = cPickle.load(papersAllExp_handle)
# papers_protsAllExp_handle = open('goa_papers_prots.pik', 'rb')
# papers_protsAllExp_dict = cPickle.load(papers_protsAllExp_handle)
# 
# # going for the top fifty papers
# top = 50
# papers_anallts2_dict_all = sp_tools.top_papers_dict(papersAllExp_dict, papers_protsAllExp_dict, top=top)
# print "really long time... (2)"
# all_tt_count_all = sp_tools.term_types_all_papers(papersAllExp_dict) #takes a really long time
# print "...done"
# go_ec_count_all = sp_tools.go_terms_with_ec_per_paper(papersAllExp_dict, top=top) # this takes a bit of time too
# allEvCodes_dict_all = sp_tools.ev_codes_all_papers(papersAllExp_dict)
# sortedProtsPerPaper_tuple_all = sp_tools.sort_papers_prots(papers_protsAllExp_dict)
# sp_tools.print_paper_per_prots_go(papers_anallts2_dict_all, all_tt_count_all, go_ec_count_all, allEvCodes_dict_all, 
#                          sortedProtsPerPaper_tuple_all, "AllPaperInfoTop50_lo.txt", top=top)
# 
# 
# 
# 
# 
# 
# 
