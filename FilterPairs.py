import argparse
import os
import sys
from time import time
from collections import defaultdict

parser = argparse.ArgumentParser(description="Given a HiC-Pro allValidPairs file, and a black list of pair ID, generate two files: fastq.allValidPairs and fastq.allValidPairs.removed")
parser.add_argument("-f","--file", type=str, help="fastq.allValidPairs, need to be filtered")
parser.add_argument("-b", "--black",type=str, help="black list of ID for valid pairs")

args = parser.parse_args()

def LoadBlackList(file):
    start_time=time()
    if file:
        f=open(file)
    blackID=defaultdict(int)
    while True:
        line=f.readline()
        #if line.startswith(('#','track','browser')):continue
        if not line:
            break
        cols=line.strip().split("\t")
        ID=cols[0]
        blackID[ID]=1
    print >> sys.stderr, 'time cost for LoadBlackList:', time()-start_time
    return blackID

def FilterPairs(file,blackID):
    start_time=time()
    if file:
        f=open(file)
        outfile1="fastq.allValidPairs"
        outfile2="fastq.AVP.removed"
        of1=open(outfile1,'w')
        of2=open(outfile2,'w')
    while True:
        line=f.readline()
        #if line.startswith(('#','track','browser')):continue
        if not line:
            break
        cols=line.strip().split("\t")
        id=cols[0]
        if id in blackID:
            of2.write(line)
        else:
            of1.write(line)
    print >> sys.stderr, 'time cost for FilsterRef:', time()-start_time


blackID=LoadBlackList(args.black)

FilterPairs(args.file,blackID)


