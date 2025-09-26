#!/usr/bin/env python2

"""
Run sc_atac_window_counter.py over many BAMs in parallel.

You can use either a *window BED* (5kb windows) or per-sample MACS2
narrowPeak files (default). For peaks, this script will try to find
<peaks_root>/<prefix>/<prefix>_peaks.narrowPeak where
prefix = basename(BAM).split('.bam')[0].

Example (peaks):
  python2 run_window_counter_parallel.py \
    /scratch/.../bam \
    /scratch/.../srr2gsm.tsv \
    /scratch/.../indextable \
    /scratch/.../peaks \
    /scratch/.../matrix_out \
    -j 12
"""

import os, sys, glob, subprocess
from multiprocessing.dummy import Pool as ThreadPool

def read_srr2gsm(path):
    m = {}
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith('#'): continue
            parts = [p for p in ln.replace(',', ' ').split() if p]
            if len(parts) >= 2:
                m[parts[0]] = parts[1]
    return m

def pick_indextable(indextable_dir, gsm):
    pats = [os.path.join(indextable_dir, "%s*.indextable.txt" % gsm),
            os.path.join(indextable_dir, "%s*.indextable.txt.gz" % gsm)]
    hits = []
    for p in pats: hits.extend(sorted(glob.glob(p)))
    return hits[0] if hits else None

def peak_for_bam(peak_root, bam):
    pref = os.path.splitext(os.path.basename(bam))[0]
    d = os.path.join(peak_root, pref)
    cand = os.path.join(d, "%s_peaks.narrowPeak" % pref)
    return cand

def run_one(t):
    (bam, srr2gsm, indextable_dir, region_arg, outdir, use_windows, pyexe, tool_path) = t

    srr = os.path.basename(bam).split('.')[0]
    gsm = srr2gsm.get(srr)
    if not gsm:
        return (bam, 2, "No GSM for %s" % srr)

    idx = pick_indextable(indextable_dir, gsm)
    if not idx:
        return (bam, 3, "No indextable for %s" % gsm)

    if use_windows:
        bed = region_arg
        if not os.path.isfile(bed):
            return (bam, 4, "Missing windows BED")
        out = os.path.join(outdir, "%s.windowmatrix.txt" % os.path.splitext(os.path.basename(bam))[0])
    else:
        bed = peak_for_bam(region_arg, bam)
        if not os.path.isfile(bed):
            return (bam, 5, "Missing peaks: %s" % bed)
        out = os.path.join(outdir, "%s.peakmatrix.txt" % os.path.splitext(os.path.basename(bam))[0])

    if not os.path.isdir(outdir):
        try: os.makedirs(outdir)
        except OSError: pass

    cmd = [pyexe, tool_path, bam, idx, bed, out, "True"]
    print("[RUN]", " ".join(cmd))
    try:
        rc = subprocess.call(cmd)
        return (bam, rc, "exit %d" % rc)
    except Exception as e:
        return (bam, 1, "ERR: %s" % str(e))

def main():
    import argparse
    p = argparse.ArgumentParser(description="Parallel sc_atac_window_counter.py (Python 2.7).")
    p.add_argument("bam_dir")
    p.add_argument("srr2gsm")
    p.add_argument("indextable_dir")
    p.add_argument("region_arg", help="Either peaks_root (dir) OR windows.bed (file)")
    p.add_argument("outdir")
    p.add_argument("-j","--jobs", type=int, default=8)
    p.add_argument("--windows", action="store_true", help="Treat region_arg as windows BED")
    p.add_argument("--pattern", default="*.nodups.bam")
    p.add_argument("--python", default="python2")
    args = p.parse_args()

    tool_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sc_atac_window_counter.py")
    if not os.path.isfile(tool_path):
        sys.exit("sc_atac_window_counter.py not found")

    bams = []
    for pat in args.pattern.split(','):
        bams.extend(sorted(glob.glob(os.path.join(args.bam_dir, pat.strip()))))
    if not bams:
        sys.exit("No BAMs found")

    srr2gsm = read_srr2gsm(args.srr2gsm)

    tasks = [(b, srr2gsm, args.indextable_dir, args.region_arg, args.outdir,
              args.windows, args.python, tool_path) for b in bams]

    pool = ThreadPool(args.jobs)
    try:
        res = pool.map(run_one, tasks)
    finally:
        pool.close(); pool.join()

    ok = sum(1 for _, rc, _ in res if rc == 0)
    print("\n[SUMMARY] %d/%d completed" % (ok, len(res)))
    for b, rc, msg in res:
        tag = "OK" if rc == 0 else "ERR"
        print("[%s] %s -> %s" % (tag, os.path.basename(b), msg))

if __name__ == "__main__":
    main()

