#!/usr/bin/env python
"""
version 0.5

#TODO:
minor revisions:
 - FTP and argparse
"""
import argparse
from Bio import Entrez
import sys
import os
import time
import subprocess
## dd/mm/yyyy format
datetimetag = time.strftime("%Y%m%d%I%M%S")
#


def get_args(DEBUG=False):
    dummy_args = argparse.Namespace(inputlist="accessions.txt",  # "batch or single
                                    outputdir=os.getcwd(),
                                    outfmt='fasta',
                                    email="Joe@salad.edu",
                                    concat=True)
    if DEBUG:  # ie, if running interactively
        print("Using Dummy Arguments")
        return(dummy_args)
    parser = argparse.ArgumentParser(
        description=str("A script to fetch nucleotide sequences " +
                        "from NCBI when given a file containing " +
                        "NCBI accession numbers."))
    parser.add_argument("-i", "--inputlist", action="store", type=str,
                        default="accessions.txt",
                        help="file containing list of accessions or " +
                        "ftp addresses; default: %(default)s")
    parser.add_argument("-q","--quick_fetch", action="store",
                        dest="quick_fetch", default='',
                        help="if used, this will just get single accession "+
                        "provided as argument. Ignores any provided accession"+
                        " list")
    parser.add_argument("-o", "--outputdir", action="store", dest="outputdir",
                        help="where to put the output; default: %(default)s",
                        default=os.getcwd(), type=str)
    parser.add_argument("-f", "--outfmt", action="store", dest="outfmt",
                        help="fasta or gb; default: %(default)s",
                        default='fasta', type=str)
    parser.add_argument("-e", "--email", action="store", dest="email",
                        help="email address",
                        default='joe_smith@mail.gov', type=str)
    parser.add_argument("--concat", action="store_true", dest="concat", default=False,
                        help="output as single file; default: %(default)s")
    parser.add_argument("--DEBUG", action="store_true", dest="DEBUG", default=False,
                        help="DEBUG mode for testing; default: %(default)s")
    args = parser.parse_args()
    if args.DEBUG:  # ie, if running debug from cline
        print("Using Dummy Arguments")
        return(dummy_args)
    else:
        return(args)


def parse_accession_list(pathtofile):
    accessions = []
    for line in open(pathtofile):
        li = line.strip()
        if not li.startswith("#"):
            accessions.append(line.rstrip())
    return(accessions)


def fetch_and_write_seqs(accessions, destination, outfmt='fasta', concat=False):
    """
    now should be able to detect and handle ftp calls to NCBI only
    takes an accession list
    if first item listed is an ftp address, it will run FTP calls for all the
    rest of the items. if not, the files will be fetched from entrex and if
    concat, smooshed into one file. Returns 0 if succesful
    """
    if outfmt == 'fasta':
        ftpsuf = "_genomic.fna.gz"  # ftp suffix
        rettype = 'fasta'
    elif outfmt == "gb":
        ftpsuf = "_genomic.gbff.gz"  # ftp suffix
        rettype = "gbwithparts"
    else:
        print("only supports fasta and gb")
        sys.exit(1)
    if accessions[0].startswith("ftp://"):
        if concat:
            raise ValueError(str("FTP input can only output individual files for" +
                                 "now; turn off --concat flag"))
        ftpcmds = []
        for i in accessions:
            filebasename = str(i.split('/')[-1] + ftpsuf)
            fastapath = str(i + "/" + filebasename)
            cmd = "wget %s -O %s" % (fastapath,
                                     os.path.join(destination, filebasename))
            ftpcmds.append(cmd)
        try:
            for j in ftpcmds:
                print("Fetching item %i of %i" % (ftpcmds.index(j) + 1,
                                                  len(ftpcmds)))
                subprocess.call(j, shell=sys.platform != "Win32")
            return(0)
        except:
            print("Somethings gone sour")
            sys.exit(1)
    if concat:
        for i in accessions:
            print("fetching %s as genomic %s" % (i, outfmt))
            out_handle = open(str(destination.strip() + datetimetag +
                                  "get_genomes_result." + outfmt), "a")
            sequence_handle = Entrez.efetch(db="nucleotide", id=i, rettype=rettype, retmode="text")
            for line in sequence_handle:
                out_handle.write(line)
            out_handle.close()
            sequence_handle.close()
    else:
        for i in accessions:
            print("fetching %s as %s" % (i, outfmt))
            out_handle = open(str(destination.strip() + i + "." + outfmt), "w")
            sequence_handle = Entrez.efetch(db="nucleotide", id=i, rettype=rettype, retmode="text")
            for line in sequence_handle:
                out_handle.write(line)
            out_handle.close()
            sequence_handle.close()
    return(0)


#%%
##########################################################################

if __name__ == '__main__':
    args = get_args()
    Entrez.email = args.email
    output_dir_path = os.path.join(args.outputdir, "")
    if not os.path.isdir(output_dir_path):
        print("creating %s" % output_dir_path)
        os.mkdir(output_dir_path)
    if args.quick_fetch is not "":
        accession_list = [args.quick_fetch]
    else:
        accession_list = parse_accession_list(args.inputlist)
    outputpath = fetch_and_write_seqs(accession_list,  output_dir_path,
                                      outfmt=args.outfmt, concat=args.concat)
    print("Outputting results to %s" % (output_dir_path))
