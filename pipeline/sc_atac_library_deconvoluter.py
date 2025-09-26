#!/usr/bin/env python2

import sys, os, gzip
import pysam

# Usage: python sc_atac_library_deconvoluter.py [Input Bam] [Index table] [Output prefix] [Output extension]
if len(sys.argv) != 5:
    sys.exit('Usage: python sc_atac_library_deconvoluter.py [Input Bam file] [Input Index table] [Output prefix] [Output extension]')

inbam_path = sys.argv[1]
index_path = sys.argv[2]
out_prefix = sys.argv[3]
out_ext = sys.argv[4]

inbam = pysam.Samfile(inbam_path,'rb')

# open index (gz or plain)
def openg(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'r')

inindex = openg(index_path)

libdic = {}
outdic = {}

for line in inindex:
    liner = line.strip().split()
    if not liner: continue
    # first column: barcode; second: experiment/group name (may include underscores; keep as-is)
    libdic[liner[0]] = liner[1]
    if liner[1] not in outdic:
        out_path = "%s.%s%s" % (out_prefix, liner[1], out_ext)
        outdic[liner[1]] = pysam.Samfile(out_path, 'wb', template=inbam)

inindex.close()

readdic = {}
for aln in inbam.fetch():
    try:
        currlib = libdic[aln.qname.split(':')[0]]
        outdic[currlib].write(aln)
        readdic[currlib] = readdic.get(currlib, 0) + 1
    except KeyError:
        # barcode not in index table
        continue

inbam.close()

for name in readdic.keys():
    print name + "\t" + str(readdic[name])
    outdic[name].close()

# close any outputs that might be empty but were opened
for name, fh in outdic.items():
    try:
        fh.close()
    except:
        pass

