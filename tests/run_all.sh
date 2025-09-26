# Convert fastq files to deduplicated bam files

./run_sc_atac_fastq2bam_parallel.py \
  -F /scratch/bx2ur/mouse_atac_atlas/fastq \
  -O /scratch/bx2ur/mouse_atac_atlas/bam \
  -G /scratch/bx2ur/mouse_atac_atlas/mm10/0f10d83b1050c08dd53189986f60970b92a315aa7a16a6f1 \
  -j 20

# Deconvolute bam file into barcodes from different experiments

python2 run_sc_atac_library_deconvoluter_parallel.py \
  /scratch/bx2ur/mouse_atac_atlas/bam \
  ./srr2gsm.txt \
  ./indextable \
  -j 20 


# Determine read depth cutoff for cells and then generate window matrix

python2 run_sc_atac_bam2matrix_parallel.py \
  /scratch/bx2ur/mouse_atac_atlas/bam \
  ./srr2gsm.txt \
  ./indextable \
  ./mm10.5kb.windows.bed \
  /scratch/bx2ur/mouse_atac_atlas/matrix_out \
  -j 20 \
  -C auto \
  --pattern "*.nodups.bam"



# run peak calling

python2 run_macs2_peaks_parallel.py \
  /scratch/bx2ur/mouse_atac_atlas/bam \
  /scratch/bx2ur/mouse_atac_atlas/peak_loose_qval \
  -g mm \
  -j 20 \
  --pattern "*.nodups.bam" \
  --format BAMPE \
  --qval 0.05 \
  --call-summits

# generate peak matrix for defined peaks

python2 run_window_counter_parallel.py \
  /scratch/bx2ur/mouse_atac_atlas/bam \
  /scratch/bx2ur/mouse_atac_atlas/matrix_out/peak_matrix_loose_qval \
  --index ./indextable/master.indextable.txt \
  --peaks-root /scratch/bx2ur/mouse_atac_atlas/peak_loose_qval \
  -j 20
