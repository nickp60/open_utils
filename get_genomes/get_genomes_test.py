#!/usr/bin/env python
"""
version 0.1
minor revisions:
 -initial version
""" 
#TODO implement disabling of sortTestMethodsUsing
import unittest
import os
import urllib2
from Bio import Entrez
import sys


sys.dont_write_bytecode = True


# change dir to dir of script
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

from get_genomes import parse_accession_list, fetch_and_write_seqs


class GetGenomesTestCase(unittest.TestCase):
    """Tests for `get_genomes.py`."""


    def test_accessionlist_is_greaterthanzero(self):
        """Make sure this parses out 2 accessions"""        
        self.assertTrue(2==len(parse_accession_list("accessions.txt")))


    def test_entrez_canreach_ncbi(self):
        """test that ncbi can be reached ok
        """
        def make_urlcall(acc="LC075482.1"):
            try:
                return(Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text"))
            except urllib2.HTTPError:
                return(1)

        test=make_urlcall()
        self.assertTrue(test!=1)


    def test_fetchandwriteseqs_makes_output(self):
        """Does the output have the right number of lines?
        """
        with open(outputpath, "r") as doc:
            test_doc=doc.readlines()
        self.assertTrue(len(test_doc)==56) # determined emperically


    def test_zcleanup_if_done(self):
        print(outputpath)
        os.remove(outputpath)  # bye bye temp file from tests. Begins with z cause of the sorting :(. see sortTestMethodsUsing if you are feeling ambitions
        self.assertTrue(not os.path.isfile(outputpath))
        

if __name__ == '__main__':
    acc_list=parse_accession_list("accessions.txt")
    outputpath = fetch_and_write_seqs(acc_list, os.path.join(dname,""))
    unittest.main()

