#!/usr/bin/env python

import sp_tools as spt
import sys

# Get a list of files. Each file has data in the following format:
# First file: (papers_go_info_content.txt)
# GO_ID   Information_content TERM_TYPE

# Next files: (top_papers_go_depth_pp.txt)
#GO_ID   TERM_TYPE     DEPTH    PROT_ID


if __name__ == '__main__':
    infiles = sys.argv[1:]
    spt.go_info_content_stats_all(infiles)

