#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Run sc_atac_bam2matrix.py over many BAMs in parallel (Python 2.7),
using the original per-GSM indextable files.

USAGE
-----
python2 run_sc_atac_bam2matrix_parallel.py \
  /path/to/bam_dir \
  /path/to/srr2gsm.tsv \
  /path/to/indextable_dir \
  /path/to/windows.bed \
  /path/to/matrix_out \
  -j 10 -C auto --pattern "*.nodups.bam"
"""

from __future__ import print_function
import argparse
import os
import glob
import subprocess
from multiprocessing.dummy import Pool as ThreadPool
import tempfile
import shutil

def gunzip_to_temp(gz_path, workdir):
    """Decompress .gz to a temp .txt inside workdir; return the .txt path."""
    # make a per-run temp dir under the output dir (or fallback to system tmp)
    base_dir = workdir if os.path.isdir(workdir) else None
    tmpdir = tempfile.mkdtemp(prefix="idx_", dir=base_dir)
    txt_path = os.path.join(tmpdir, os.path.basename(gz_path).replace(".gz",""))
    # decompress using Python to avoid shell deps
    import gzip
    with gzip.open(gz_path, "rb") as fin, open(txt_path, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    return txt_path


def read_srr2gsm(path):
    """Read two-column SRR->GSM map (whitespace or comma separated)."""
    m = {}
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith('#'):
                continue
            parts = [p for p in ln.replace(',', ' ').split() if p]
            if len(parts) >= 2:
                m[parts[0]] = parts[1]
    return m

def pick_indextable(indextable_dir, gsm):
    """Find GSM*.indextable.txt[.gz] in indextable_dir."""
    pats = [
        os.path.join(indextable_dir, "%s*.indextable.txt" % gsm),
        os.path.join(indextable_dir, "%s*.indextable.txt.gz" % gsm),
    ]
    hits = []
    for pat in pats:
        hits.extend(sorted(glob.glob(pat)))
    # de-duplicate preserving order
    seen = set()
    uniq = []
    for p in hits:
        if p not in seen:
            uniq.append(p); seen.add(p)
    if not uniq:
        return None
    if len(uniq) > 1:
        # Just use the first; warn on stderr
        try:
            import sys
            sys.stderr.write("[WARN] Multiple indextables for %s; using %s\n" % (gsm, uniq[0]))
        except Exception:
            pass
    return uniq[0]

def find_bams(bam_dir, pattern_csv):
    bam_dir = os.path.abspath(bam_dir)
    patterns = [p.strip() for p in pattern_csv.split(",") if p.strip()]
    found = []
    for pat in patterns:
        found.extend(sorted(glob.glob(os.path.join(bam_dir, pat))))
    # de-duplicate while preserving order
    seen = set()
    uniq = []
    for b in found:
        if b not in seen:
            uniq.append(b)
            seen.add(b)
    return uniq

def srr_from_bam(bam_path):
    """
    Extract SRR accession from BAM filename.
    Works for names like:
      SRR6819229.sort.Liver.nodups.bam
      SRR6819229.split.q10.sort.bam
    """
    base = os.path.basename(bam_path)
    i = base.find('.')
    return base[:i] if i > 0 else base

def run_one(args_tuple):
    bam, srr2gsm, indextable_dir, windows_bed, outdir, cutoff, pyexe, tool_path = args_tuple

    if not os.path.isfile(windows_bed):
        return (bam, 1, "windows BED not found: %s" % windows_bed)

    if not os.path.isdir(outdir):
        try:
            os.makedirs(outdir)
        except OSError:
            pass  # ok if it exists

    srr = srr_from_bam(bam)
    gsm = srr2gsm.get(srr)
    if not gsm:
        return (bam, 2, "SRR→GSM mapping missing for %s" % srr)

    index_table = pick_indextable(indextable_dir, gsm)
    if not index_table:
        return (bam, 3, "No indextable for %s under %s" % (gsm, indextable_dir))
    # If the index table is gzipped, decompress to a temp .txt for tools that use open(...)
    
    if index_table.endswith(".gz"):
        try:
            # use the output dir as a parent for temp files (keeps permissions sane)
            index_table_txt = gunzip_to_temp(index_table, outdir)
        except Exception as e:
            return (bam, 4, "Failed to gunzip %s: %s" % (index_table, str(e)))
    else:
        index_table_txt = index_table

    prefix = os.path.splitext(os.path.basename(bam))[0]  # keep full BAM stem
    cmd = [
        pyexe, tool_path,
        "-B", bam,
        "-I", index_table_txt,
        "-O", outdir,         # OUTDIR is a directory
        "-P", prefix,
        "-C", str(cutoff),
        "-W", windows_bed,    # WINDOW BED is a file
    ]
    print("[RUN] " + " ".join(cmd))
    try:
        ret = subprocess.call(cmd)
        if ret != 0:
            return (bam, ret, "Exited with code %d" % ret)
        return (bam, 0, "OK (%s)" % os.path.basename(index_table))
    except Exception as e:
        return (bam, 1, "Error: %s" % str(e))

def parse_args():
    p = argparse.ArgumentParser(
        description="Run sc_atac_bam2matrix.py over many BAMs in parallel (Python 2.7) using per-GSM indextables."
    )
    p.add_argument("bam_dir", help="Directory containing input BAM files")
    p.add_argument("srr2gsm", help="Two-column SRR→GSM mapping file")
    p.add_argument("indextable_dir", help="Directory with GSM*.indextable.txt[.gz]")
    p.add_argument("windows_bed", help="Path to windows BED (genome 'universe')")
    p.add_argument("outdir", help="Directory to write outputs")
    p.add_argument("-j", "--jobs", type=int, default=8, help="Parallel workers (default: 8)")
    p.add_argument("-C", "--cutoff", default="auto",
                   help='Read-depth cutoff passed to sc_atac_bam2matrix.py (default: "auto"; or give a number)')
    p.add_argument("--pattern", default="*.nodups.bam",
                   help='Comma-separated glob(s) to select BAMs (default: "*.nodups.bam")')
    p.add_argument("--python", default="python2",
                   help='Python executable to run sc_atac_bam2matrix.py (default: "python2")')
    return p.parse_args()

def main():
    args = parse_args()

    # Resolve path to sc_atac_bam2matrix.py next to this runner
    script_dir = os.path.dirname(os.path.abspath(__file__))
    tool_path = os.path.join(script_dir, "sc_atac_bam2matrix.py")
    if not os.path.isfile(tool_path):
        print("[ERROR] sc_atac_bam2matrix.py not found at: %s" % tool_path)
        return

    bams = find_bams(args.bam_dir, args.pattern)
    if not bams:
        print("[WARN] No BAMs found in %s with patterns: %s" % (args.bam_dir, args.pattern))
        return

    srr2gsm = read_srr2gsm(args.srr2gsm)

    print("[INFO] Found %d BAM(s)." % len(bams))
    print("[INFO] Output -> %s" % os.path.abspath(args.outdir))
    print("[INFO] Windows BED -> %s" % os.path.abspath(args.windows_bed))
    print("[INFO] Index dir -> %s" % os.path.abspath(args.indextable_dir))

    job_args = [
        (bam, srr2gsm, args.indextable_dir, args.windows_bed, args.outdir,
         args.cutoff, args.python, tool_path)
        for bam in bams
    ]

    pool = ThreadPool(args.jobs)
    try:
        results = pool.map(run_one, job_args)
    finally:
        pool.close()
        pool.join()

    # Summary
    ok = [r for r in results if r[1] == 0]
    bad = [r for r in results if r[1] != 0]
    print("\n[SUMMARY] %d succeeded, %d failed." % (len(ok), len(bad)))
    if bad:
        for bam, code, msg in bad:
            print("  - %s: %s" % (bam, msg))

if __name__ == "__main__":
    main()

