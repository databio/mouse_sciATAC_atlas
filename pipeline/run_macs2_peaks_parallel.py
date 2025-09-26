#!/usr/bin/env python2
from __future__ import print_function
import os, sys, argparse, glob, subprocess
from multiprocessing.dummy import Pool as ThreadPool

def parse_args():
    p = argparse.ArgumentParser(
        description="Run MACS2 callpeak on many BAMs in parallel (Python 2.7)."
    )
    p.add_argument("bam_dir", help="Directory with input BAMs")
    p.add_argument("outdir", help="Directory to write MACS2 outputs")
    p.add_argument("-g", "--genome", default="mm",
                   help='MACS2 -g/--gsize (e.g. "mm", "hs", or numeric)')
    p.add_argument("-j", "--jobs", type=int, default=8,
                   help="Parallel workers (default: 8)")
    p.add_argument("--pattern", default="*.nodups.bam",
                   help='Comma-separated globs for BAMs (default: "*.nodups.bam")')
    p.add_argument("--format", dest="bam_format", default="BAMPE",
                   choices=["BAM", "BAMPE"],
                   help="MACS2 -f format (BAMPE for paired-end)")
    p.add_argument("--qval", default="0.05", help="MACS2 -q (FDR) (default: 0.05)")
    p.add_argument("--keep-dup", default="all", help='MACS2 --keep-dup (default: "all")')
    p.add_argument("--call-summits", action="store_true", help="Pass --call-summits")
    p.add_argument("--extra", default="", help="Extra args appended verbatim")
    return p.parse_args()

def find_bams(bam_dir, pattern_csv):
    bam_dir = os.path.abspath(bam_dir)
    pats = [p.strip() for p in pattern_csv.split(",") if p.strip()]
    out = []
    for pat in pats:
        out.extend(sorted(glob.glob(os.path.join(bam_dir, pat))))
    # de-dupe preserve order
    seen, uniq = set(), []
    for b in out:
        if b not in seen:
            uniq.append(b); seen.add(b)
    return uniq

def run_one(task):
    (bam, outdir, gsize, bam_format, qval, keepdup, call_summits, extra) = task
    base = os.path.splitext(os.path.basename(bam))[0]
    od = os.path.join(outdir, base)
    if not os.path.isdir(od):
        try: os.makedirs(od)
        except OSError: pass
    cmd = [
        "macs2", "callpeak",
        "-t", bam,
        "-f", bam_format,
        "-g", str(gsize),
        "-q", str(qval),
        "--keep-dup", str(keepdup),
        "-n", base,
        "--outdir", od,
        "-B", "--SPMR"
    ]
    if call_summits: cmd.append("--call-summits")
    if extra:
        cmd.extend(extra.split())
    print("[RUN]", " ".join(cmd))
    try:
        rc = subprocess.call(cmd)
        return (bam, rc)
    except Exception as e:
        sys.stderr.write("ERROR running macs2: %s\n" % str(e))
        return (bam, 1)

def main():
    a = parse_args()
    bams = find_bams(a.bam_dir, a.pattern)
    if not bams:
        sys.stderr.write("No BAMs found under %s with %s\n" % (a.bam_dir, a.pattern))
        sys.exit(2)
    if not os.path.isdir(a.outdir):
        try: os.makedirs(a.outdir)
        except OSError: pass
    tasks = [(b, a.outdir, a.genome, a.bam_format, a.qval, a.keep_dup, a.call_summits, a.extra)
             for b in bams]
    pool = ThreadPool(a.jobs)
    try:
        results = pool.map(run_one, tasks)
    finally:
        pool.close(); pool.join()
    ok = sum(1 for _, rc in results if rc == 0)
    print("[SUMMARY] %d/%d succeeded" % (ok, len(results)))

if __name__ == "__main__":
    main()

