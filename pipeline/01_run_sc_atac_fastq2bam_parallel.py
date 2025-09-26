#!/usr/bin/env python2

import os, sys, argparse, glob, subprocess
from multiprocessing.dummy import Pool as ThreadPool

def parse_args():
    p = argparse.ArgumentParser(description="Parallel runner for sc_atac_fastq2bam.py")
    p.add_argument("fastq_dir", help="Dir with *_1.fastq.gz and *_2.fastq.gz")
    p.add_argument("outdir", help="Output BAM dir")
    p.add_argument("-j","--jobs", type=int, default=8)
    p.add_argument("--python", default="python2")
    p.add_argument("--tool", default="sc_atac_fastq2bam.py")
    p.add_argument("--extra", default="", help="Extra args to pass")
    return p.parse_args()

def pair_fastqs(fastq_dir):
    fastq_dir = os.path.abspath(fastq_dir)
    r1s = sorted(glob.glob(os.path.join(fastq_dir, "*_1.fastq.gz")))
    pairs = []
    for r1 in r1s:
        r2 = r1.replace("_1.fastq.gz", "_2.fastq.gz")
        if os.path.isfile(r2):
            pairs.append((r1, r2))
    return pairs

def run_one(task):
    r1, r2, outdir, py, tool, extra = task
    srr = os.path.basename(r1).split("_")[0]
    if not os.path.isdir(outdir):
        try: os.makedirs(outdir)
        except OSError: pass
    cmd = [py, tool,
           "-R1", r1, "-R2", r2,
           "-O", outdir, "-P", srr]
    if extra: cmd.extend(extra.split())
    print("[RUN]", " ".join(cmd))
    return subprocess.call(cmd)

def main():
    a = parse_args()
    pairs = pair_fastqs(a.fastq_dir)
    if not pairs:
        sys.stderr.write("No *_1.fastq.gz pairs found in %s\n" % a.fastq_dir)
        sys.exit(2)
    pool = ThreadPool(a.jobs)
    try:
        rcodes = pool.map(lambda pr: run_one(pr + (a.outdir, a.python, a.tool, a.extra)),
                          pairs)
    finally:
        pool.close(); pool.join()
    ok = sum(1 for r in rcodes if r == 0)
    print("[SUMMARY] %d/%d succeeded" % (ok, len(rcodes)))

if __name__ == "__main__":
    main()

