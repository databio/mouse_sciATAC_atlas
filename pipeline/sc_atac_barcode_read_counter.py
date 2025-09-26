import sys
import pysam
import os, subprocess

def ensure_bam_index(bam_path):
    """
    Ensure a BAM index exists. If not, run `samtools index <bam>`.
    Checks both .bai and .csi variants.
    """
    bai1 = bam_path + '.bai'
    bai2 = os.path.splitext(bam_path)[0] + '.bai'
    csi  = bam_path + '.csi'

    if os.path.exists(bai1) or os.path.exists(bai2) or os.path.exists(csi):
        return

    sys.stderr.write('[INFO] No BAM index found for %s; running: samtools index\n' % bam_path)
    try:
        subprocess.check_call(['samtools', 'index', bam_path])
    except OSError:
        sys.stderr.write('[ERROR] samtools not found in PATH. Please load a module or add it to PATH.\n')
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        sys.stderr.write('[ERROR] samtools index failed (exit %d) for %s\n' % (e.returncode, bam_path))
        sys.exit(1)

if len(sys.argv) != 4:
    sys.exit('Usage: python sc_atac_samespecies_individual_readcounter.py [Input Bam file] [Input Index table (NA if no table)] [Out file]')

inputbam = sys.argv[1]
indextable = sys.argv[2]
outfile = sys.argv[3]

totalct = {}
descriptor = {}
descdic = {}

# ---- minimal robust parsing of the index table ----
if indextable != 'NA':
    with open(indextable, 'r') as descer:
        for line in descer:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            liner = line.split()
            # if the line has only a barcode, assign default group "bkgd"
            if len(liner) >= 2:
                descdic[liner[0]] = liner[1]
            elif len(liner) == 1:
                descdic[liner[0]] = 'bkgd'
# ---------------------------------------------------

ensure_bam_index(inputbam)

print 'Counting total reads...'
bamfile = pysam.Samfile(inputbam, 'rb')
for read in bamfile.fetch():
    # safer split (in case there is no ':')
    qn = read.qname
    tagger = qn.split(':', 1)[0] if ':' in qn else qn
    try:
        totalct[tagger] += 1
    except KeyError:
        totalct[tagger] = 1
        try:
            descriptor[tagger] = descdic[tagger]
        except KeyError:
            descriptor[tagger] = 'bkgd'

bamfile.close()

outter = open(outfile, 'w')
print >> outter, 'Barcode\tExperiment\tReadCount'
for tag in sorted(totalct.keys()):
    print >> outter, tag + '\t' + descriptor[tag] + '\t' + str(totalct[tag])
outter.close()

