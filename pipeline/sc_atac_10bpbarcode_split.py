#!/usr/bin/env python2

import argparse
import subprocess
import sys
sys.path.append('~/bin/Levenshtein/')
import Levenshtein
import gzip
import io
import cStringIO
io_method = cStringIO.StringIO

parser = argparse.ArgumentParser(description='Fix/normalize barcodes in scATAC FASTQ headers.')
parser.add_argument('-1','--Read1', dest='input1', required=True)
parser.add_argument('-2','--Read2', dest='input2', required=True)
parser.add_argument('-O1','--output1', dest='output1', required=True)
parser.add_argument('-O2','--output2', dest='output2', required=True)
parser.add_argument('-L','--log', dest='logfile', required=True)
parser.add_argument('-X','--nextseq', dest='nextseq', action="store_true")
parser.add_argument('-Z','--gzip', dest='gzip', action="store_true")
args = parser.parse_args()

def submitter(cmd):
    p = subprocess.Popen(cmd, shell=True); p.wait()

def editcheck(barc,reflist):
    try:
        reflist[barc]; eddist = '0'
    except KeyError:
        winner = '_CTF' + '_'*(len(barc)-4)
        winner_ed = 10; runnerup_ed = 10
        for barcode in reflist.keys():
            curred = Levenshtein.distance(barc,barcode)
            if curred <= winner_ed:
                runnerup_ed = winner_ed; winner = barcode; winner_ed = curred
            if (curred > winner_ed) and (curred < runnerup_ed):
                runnerup_ed = curred
        if winner_ed > 3: winner = '_CTF' + '_'*(len(barc)-4)
        if runnerup_ed - winner_ed < 2: winner = '_AMBIG' + '_'*(len(barc)-6)
        barc = winner; eddist = str(winner_ed)
    return (barc, eddist)

complements = {'A':'T', 'T':'A', 'G':'C', 'C':'G','N':'N'}
def reverse_complement(x):
    return ''.join([complements.get(z,'N') for z in x[::-1]])

# (same barcode dictionaries as before)
nex_i7 = {"ATTACTCG":"","TCCGGAGA":"","CGCTCATT":"","GAGATTCC":"","ATTCAGAA":"","GAATTCGT":"","CTGAAGCT":"","TAATGCGC":"","CGGCTATG":"","TCCGCGAA":"","TCTCGCGC":"","AGCGATAG":""}
# (pcr_i7, pcr_i5, nex_i5 omitted here for brevity – use your full dicts)
pcr_i7 = {...}
pcr_i5 = {...}
nex_i5 = {"TATAGCCT":"","ATAGAGGC":"","CCTATCCT":"","GGCTCTGA":"","AGGCGAAG":"","TAATCTTA":"","CAGGACGT":"","GTACTGAC":""}

FULL_BARC = {}
for i in nex_i7.keys():
    for j in pcr_i7.keys():
        for k in pcr_i5.keys():
            for l in nex_i5.keys():
                FULL_BARC[i + j + k + l] = ""

totreads = exactmatch = editmatch = failed = readcount = 0
outfasta1 = open(args.output1,'w')
outfasta2 = open(args.output2,'w')

if1 = gzip.open(args.input1,'rb')
infasta1 = io.BufferedReader(if1)
if2 = gzip.open(args.input2,'rb')
infasta2 = io.BufferedReader(if2)

reads1 = []; reads2 = []

while True:
    h1 = infasta1.readline()
    if not h1: break
    seq1 = infasta1.readline().strip()
    plus1 = infasta1.readline()
    qual1 = infasta1.readline().strip()

    h2 = infasta2.readline()
    seq2 = infasta2.readline().strip()
    plus2 = infasta2.readline()
    qual2 = infasta2.readline().strip()

    h1s = h1.strip()
    # Original header pattern: "@<id> <machine>:..._<i7>-<i5> length=50"
    # Extract the last whitespace-separated token from the header *ID* part
    # or fall back to the second field and split on '_' then '-'.
    try:
        # token like "..._I7-I5"
        token = h1s.split()[1]
        tags = token.split('_', 1)[1]
    except Exception:
        # if unexpected, put everything after first space
        tags = h1s.split()[-1]
    parts = tags.split('-')
    if len(parts) == 2:
        i7, i5 = parts[0], parts[1]
    else:
        i7, i5 = "", ""

    seqbarc = i7 + i5
    if args.nextseq:
        seqbarca = seqbarc[0:18]
        seqbarcb = reverse_complement(seqbarc[18:36])
        seqbarc = seqbarca + seqbarcb

    try:
        FULL_BARC[seqbarc]
        reads1.append('@' + seqbarc + ':' + str(totreads) + '#0000/1\n' + seq1 + '\n+\n' + qual1 + '\n')
        reads2.append('@' + seqbarc + ':' + str(totreads) + '#0000/2\n' + seq2 + '\n+\n' + qual2 + '\n')
        readcount += 1; totreads += 1; exactmatch += 1
    except KeyError:
        b1 = seqbarc[0:8]; b2 = seqbarc[8:18]; b3 = seqbarc[18:28]; b4 = seqbarc[28:36]
        b1ed = editcheck(b1,nex_i7); b2ed = editcheck(b2,pcr_i7)
        b3ed = editcheck(b3,pcr_i5); b4ed = editcheck(b4,nex_i5)
        seqbarc2 = b1ed[0] + b2ed[0] + b3ed[0] + b4ed[0]
        seqed = b1ed[1] + b2ed[1] + b3ed[1] + b4ed[1]
        reads1.append('@' + seqbarc2 + ':' + str(totreads) + '#' + seqed + '/1\n' + seq1 + '\n+\n' + qual1 + '\n')
        reads2.append('@' + seqbarc2 + ':' + str(totreads) + '#' + seqed + '/2\n' + seq2 + '\n+\n' + qual2 + '\n')
        totreads += 1; readcount += 1
        if ('CTF' in seqbarc2) or ('AMBIG' in seqbarc2): failed += 1
        else: editmatch += 1

    if readcount == 250000:
        outfasta1.writelines(reads1); outfasta2.writelines(reads2)
        readcount = 0; reads1 = []; reads2 = []

if readcount > 0:
    outfasta1.writelines(reads1); outfasta2.writelines(reads2)

if1.close(); if2.close()
outfasta1.close(); outfasta2.close()

logout = open(args.logfile,'w')
print >> logout, 'total=' + str(totreads) + '\texact=' + str(round(float(exactmatch)*100/max(1,totreads),2)) + '%\tby_ed=' + str(round(float(editmatch)*100/max(1,totreads),2)) + '%\tfail=' + str(round(float(failed)*100/max(1,totreads),2)) + '%'
logout.close()

if args.gzip:
    submitter('gzip -f ' + args.output1 + '; gzip -f ' + args.output2)

