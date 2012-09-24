#!/usr/bin/env python

import sp_tools as spt
import sys

if __name__ == '__main__':
    infiles = sys.argv[1:]
    spt.go_info_content_stats_all(infiles)

