#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 16:41:44 2016
version 0.1
Minor version changes:

TODO:

USAGE:
 $
"""

import argparse
import os
import re
DEBUG = False
#%%
#define inputs
if DEBUG:
    fastq = "3123-1_1_trimmed.fastq.gz"
    query = "3123"
    escape = False
else:
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("fastq", help="input fastq")
    parser.add_argument("query", help="sd")
    parser.add_argument("-e", "--escape", default="F", help="escape or no")
    args = parser.parse_args()
    fastq = args.fastq
    query = args.query
    escape = args.escape

if escape in ["T", "t", "true", "True"]:
    escape = True
    query = re.escape(query)
else:
    escape = False

def search_fastq(fastq, query):
    """
    Basic python solution for fastq average length
    ideas from will: https://www.biostars.org/p/1758/
    """
    tot = 0
    if os.path.splitext(fastq)[-1] == ".gz":
        import gzip
        with gzip.open(fastq) as handle:
            for i, line in enumerate(handle):
                if i % 4 == 0:
                    try:
                        print(re.search(query, line).string)
                        tot += 1
                    except AttributeError:
                        continue
    elif os.path.splitext(fastq)[-1] == ".fastq":
        with open(fastq) as handle:
            for i, line in enumerate(handle):
                if i % 4 == 0:
                    try:
                        print(re.search(query, line).string)
                        tot += 1
                    except AttributeError:
                        continue
    else:
        raise ValueError("Must be either a fastq or fastq.gz")
    print("Total of %i hits for query %s" %
          (tot, query))
if __name__ == "__main__":
    search_fastq(fastq, query)
