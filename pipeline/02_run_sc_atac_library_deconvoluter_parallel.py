#!/usr/bin/env python2

"""
Deconvolute *.split.q10.sort.bam using the *per-GSM* indextable
(chosen via an SRR→GSM mapping file). Output BAMs are named by the
*second column* of each indextable (the script itself handles that).

Example:
  python2 run_sc_atac_library_deconvoluter_parallel.py \
    /scratch/.../bam \
    /scratch/.../srr2gsm.tsv \
    /scratch/.../indextable \
    -j 20 --python python2
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
    # de-dup
    seen, out = set(), []
    for h in hits:
        if h not in seen:
            out.append(h); seen.add(h)
    return out[0] if out else None

def derive_prefix(bam_path):
    base = os.path.basename(bam_path)
    suf = ".split.q10.sort.bam"
    if not base.endswith(suf): return (None, None)
    return (os.path.dirname(os.path.abspath(bam_path)), base[:-len(suf)])

def already_done(bam_dir, srr, out_ext):
    pat = os.path.join(bam_dir, "%s.sort.*%s" % (srr, out_ext))
    return len(glob.glob(pat)) > 0

def run_one(t):
    (bam_path, srr2gsm, indextable_dir, force, out_ext, pyexe, deconv_py) = t
    bam_dir, srr = derive_prefix(bam_path)
    if srr is None:
        return (bam_path, False, "SKIP (not *.split.q10.sort.bam)")
    gsm = srr2gsm.get(srr)
    if not gsm:
        return (srr, False, "NO_GSM_MAP")

    idx = pick_indextable(indextable_dir, gsm)
    if not idx:
        return (srr, False, "NO_INDEX")

    if (not force) and already_done(bam_dir, srr, out_ext):
        return (srr, True, "SKIP (exists)")

    out_prefix = os.path.join(bam_dir, "%s.sort" % srr)
    cmd = [pyexe, deconv_py, bam_path, idx, out_prefix, out_ext]
    print("[RUN]", " ".join(cmd))
    try:
        subprocess.check_call(cmd)
        return (srr, True, "OK (%s)" % os.path.basename(idx))
    except subprocess.CalledProcessError as e:
        return (srr, False, "FAILED (%d)" % e.returncode)

def main():
    import argparse
    p = argparse.ArgumentParser(description="Parallel deconvolution using per-GSM indextables.")
    p.add_argument("bam_dir")
    p.add_argument("srr2gsm_map")
    p.add_argument("indextable_dir")
    p.add_argument("-j","--jobs", type=int, default=8)
    p.add_argument("--force", action="store_true")
    p.add_argument("--out-ext", default=".nodups.bam")
    p.add_argument("--python", default="python2")
    args = p.parse_args()

    deconv_py = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sc_atac_library_deconvoluter.py")
    if not os.path.isfile(deconv_py):
        sys.exit("sc_atac_library_deconvoluter.py not found next to this script")

    bams = sorted(glob.glob(os.path.join(args.bam_dir, "*.split.q10.sort.bam")))
    if not bams:
        sys.exit("No *.split.q10.sort.bam in %s" % args.bam_dir)

    srr2gsm = read_srr2gsm(args.srr2gsm_map)

    tasks = [(b, srr2gsm, args.indextable_dir, args.force, args.out_ext, args.python, deconv_py)
             for b in bams]

    pool = ThreadPool(args.jobs)
    try:
        res = pool.map(run_one, tasks)
    finally:
        pool.close(); pool.join()

    ok = sum(1 for _, s, _ in res if s)
    print("\n[SUMMARY] %d OK, %d ERR" % (ok, len(res)-ok))
    for srr, s, msg in res:
        print("[%s] %s -> %s" % ("OK" if s else "ERR", srr, msg))

if __name__ == "__main__":
    main()

