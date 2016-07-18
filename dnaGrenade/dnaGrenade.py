#!/usr/bin/env python
"""
version 0.3
Minor version changes:
 - add output option for loci that fall beneath length cutoff

TODO:
- make work for protein too
- only works with phred numbers.  is this a good thing?
Fragments a DNA (multi)fasta file into pseudoreads for assembly help.
- seeded randomness?  Should this even be a thing?

Input:
-DNA multifasta (.fasta)
-read length (int)
-depth of coverage (int)
Optional Input:
-seed (int)
-output directory (path)
-quality score
-minimum sequence length
Output:
-DNA multifasta

USAGE:
 $ python dnaGrenade.py multi.fasta --read_length 250 --depth 10 --seed 27 --output
          ~/grenade/ --quality 38
"""
import os
import random
#import sys
import datetime
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

DEBUG = False
#%%
#define inputs
if DEBUG:
    input_seqs_path = "genes.fasta"
    target_read_len = 250
    target_quality = 35
    target_depth = 10
    seed = 27
    output = os.path.expanduser("~/Desktop/")
    minimum_len = 500
    shorties = True
    depth_shorties = 2
else:
    parser = argparse.ArgumentParser(description="This script takes a multifasta and fragments the\
                                                  entries according to the parameters provided,\
                                                  outputing the results as a Illumina 1.8+ fastq\
                                                  file with quality provided.")
    parser.add_argument("fasta_file", help="fasta containing CDS's or genes or something")
    parser.add_argument("-q", "--quality", help="Phred Quality from 1-40 as", default=35, type=int)
    parser.add_argument("-o", "--output", help="directory in which to place the output files",
                        default=os.getcwd())
    parser.add_argument("-r", "--read_length", default=250, type=int,
                        help="desired length of pseudoreads; default is 250")
    parser.add_argument("-d", "--depth", help="depth of pseudoread coverage",
                        default=10, type=int)
    parser.add_argument("-s", "--seed", help="seed", default=27, type=int)
    parser.add_argument("-m", "--minimum", help="minimum sequence size", default=500, type=int)
    parser.add_argument('--shorties', dest='shorties', help="use this flag to make pseudoreads for\
                        loci that fall beneath length cutoff; default is False",
                        action='store_true', default=False)
    parser.add_argument("--depth_shorties", help="depth of coverage for --shorties",
                        default=0, type=int)
    args = parser.parse_args()

    input_seqs_path = args.fasta_file
    target_read_len = args.read_length
    target_quality = args.quality
    target_depth = args.depth
    seed = args.seed
    minimum_len = args.minimum
    output = args.output
    shorties = args.shorties
    depth_shorties = args.depth_shorties


date = str(datetime.datetime.now().strftime('%Y%m%d'))
# Throw error for --shorties with no depth and depth with no --shorties
if (shorties and depth_shorties == 0) or (not shorties and depth_shorties != 0):
    raise ValueError("if generating pseudoreads of short loci with --shorties, must set a coverage \
    depth! Short sequences will have fewer true reads, and therefore must have different \
    pseudocoverage from the longer reads")


def generate_break_coords(seq_len, depth_of_coverage, target_read_lenngth):
    """ Random corrdinate generator for a given sequence length, depth and read length
    could add line at end to add end of seq to each set to always capture the end fragments
    """
    coords = {}  # recipient dictionary to return
    for i in range(0, target_depth):  # For each level of depth
        point = random.randint(0, target_read_len)  # get a random point from 0 to read_len
        these_coords = []  # recipient list for coords.  could add 0 here to get gene start
        while point <= seq_len:  # until the end of the sequence,
            these_coords.append(point)  # add point to the coord list
            point = point + target_read_len  # move up by one read length, and repeat
        coords[i] = these_coords  # key is depth level, item is coord list
    return(coords)
#%%
# testcoords = generate_break_coords(seq_len, depth, target_read_len)


def breakumup(sequence, coords):
    """takes a string, splitting coords in a dictionary;returns a list
    """
    seq_list = []
    for i in range(0, len(coords)):
        for j in range(0, int(round(len(sequence)/target_read_len))):
            try:
                seq_list.append(sequence[coords[i][j]:coords[i][j+1]])
            except IndexError:
                pass  # so that this will run all the pairs and nothing more
    return(seq_list)
    # breakumup(sequence=seq, coords=testcoords)
#%%


def convert_to_SeqRecord(seqs, name):
    seqrecordlist = []
    for i in range(0, len(seqs)):
        rec = SeqRecord(seq=seqs[i], id=str(name + str(i)), name=name,
                        description="pseudo")
        seqrecordlist.append(rec)
    return(seqrecordlist)
#%%


def main():
    output_handle = os.path.join(output,
                                 str("pseudoreads_cov_" + str(target_depth) +
                                     "_qual_" + str(target_quality) +
                                     "_len_" + str(target_read_len) + ".fastq"))
    outfile = open(output_handle, "w")
    counter = 0
    shortiescounter = 0
    skipcounter = 0
    for index, rec in enumerate(SeqIO.parse(open(input_seqs_path, "r"), 'fasta')):
        if len(rec.seq) >= minimum_len:  # for debugging
            coordinates = generate_break_coords(len(rec.seq), depth_of_coverage=target_depth,
                                                target_read_lenngth=target_read_len)
            sequences_to_write = breakumup(sequence=rec.seq, coords=coordinates)
            seqRecords = convert_to_SeqRecord(seqs=sequences_to_write, name=str(rec.name+"_"))
            for record in seqRecords:
                record.letter_annotations["phred_quality"] = [target_quality] * len(record)
                SeqIO.write(record, outfile, "fastq")
                counter = counter + 1
        elif len(rec.seq) < minimum_len and shorties:
            coordinates = []
            for i in range(0, depth_shorties):
                coordinates.append([0, len(rec.seq)])
            sequences_to_write = breakumup(sequence=rec.seq, coords=coordinates)
            seqRecords = convert_to_SeqRecord(seqs=sequences_to_write, name=str(rec.name+"_"))
            for record in seqRecords:
                record.letter_annotations["phred_quality"] = [target_quality] * len(record)
                SeqIO.write(record, outfile, "fastq")
                shortiescounter = shortiescounter + 1
        else:
            skipcounter = skipcounter + 1
    outfile.close()
    print("wrote %i fragments to %s" % (counter, output_handle))
    print("wrote %i pseudoreads for short fragments to %s" % (shortiescounter, output_handle))
    print("skipped %i short sequences" % skipcounter)


if __name__ == '__main__':
    main()
