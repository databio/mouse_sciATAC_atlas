# sc_atac_cell_cutoff.R
# Usage: Rscript sc_atac_cell_cutoff.R <outdir/prefix> <cutoff|auto>

args <- commandArgs(TRUE)
if (length(args) < 2) stop("Usage: Rscript sc_atac_cell_cutoff.R <outdir/prefix> <cutoff|auto>")
currexpt <- args[1]
cellfloor_arg <- args[2]

report_path <- paste0(currexpt, ".readcount.report.txt")
if (!file.exists(report_path)) stop("Missing report: ", report_path)

# Expect columns: Barcode\tExperiment\tReadCount
report <- read.table(report_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
if (!all(c("Barcode","Experiment","ReadCount") %in% names(report))) {
  stop("Report must have columns: Barcode, Experiment, ReadCount")
}

# Drop background rows labeled 'bkgd' in Experiment
bkgd.ind <- grep("bkgd", report$Experiment)
nobkgd <- if (length(bkgd.ind) > 0) report[-bkgd.ind, ] else report

# --- guard: handle empty or non-numeric counts --------------------------------
if (nrow(nobkgd) == 0) {
  # write empty index table and a tiny PDF, then exit quietly
  out_index <- paste0(currexpt, ".readdepth.cells.indextable.txt")
  write.table(matrix(nrow = 0, ncol = 2), out_index,
              row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  pdf(paste0(currexpt, ".results.hists.pdf"), height = 4, width = 8)
  plot.new(); title("scATAC read-depth cutoff"); mtext("No non-background barcodes", 3, 0.5)
  dev.off()
  quit(save = "no", status = 0)
}

counts <- suppressWarnings(as.numeric(nobkgd$ReadCount))
valid <- !is.na(counts)
nobkgd <- nobkgd[valid, , drop = FALSE]
counts <- counts[valid]

if (length(counts) == 0) {
  out_index <- paste0(currexpt, ".readdepth.cells.indextable.txt")
  write.table(matrix(nrow = 0, ncol = 2), out_index,
              row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  pdf(paste0(currexpt, ".results.hists.pdf"), height = 4, width = 8)
  plot.new(); title("scATAC read-depth cutoff"); mtext("All counts NA", 3, 0.5)
  dev.off()
  quit(save = "no", status = 0)
}
# -----------------------------------------------------------------------------

# Determine cutoff
if (tolower(cellfloor_arg) == "auto") {
  x <- sort(counts)
  if (length(x) < 10) {
    cellfloor <- ifelse(length(x) == 0, 0, max(1, floor(median(x))))
  } else {
    lx <- log10(pmax(1, x))
    med <- median(lx)
    hi <- x[lx >= med]
    cellfloor <- if (length(hi) == 0) max(1, floor(median(x))) else floor(quantile(hi, 0.1, names = FALSE))
  }
} else {
  cellfloor <- as.numeric(cellfloor_arg)
  if (is.na(cellfloor)) stop("Cutoff must be numeric or 'auto'.")
}

# Subset to cells and write indextable for downstream
keep <- which(counts >= cellfloor)
cells <- nobkgd[keep, , drop = FALSE]
out_index <- paste0(currexpt, ".readdepth.cells.indextable.txt")
write.table(
  cbind(as.character(cells$Barcode), as.character(cells$Experiment)),
  out_index, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE
)

# Quick PDF summary (safe when counts is empty or very small)
pdf(paste0(currexpt, ".results.hists.pdf"), height = 8, width = 10)
par(mfrow = c(1,2))
if (length(counts) > 0) {
  hist(log10(pmax(1, counts)), breaks = 60, col = "mediumseagreen",
       main = "All (no bkgd)", xlab = "Reads per barcode (log10)", las = 1, xlim = c(0,7))
  abline(v = log10(max(1, cellfloor)), lty = "dashed", lwd = 2)
  text(x = 1, y = par("usr")[4]*0.9, adj = 0,
       labels = paste0("Cutoff = ", cellfloor, "\nCells = ", nrow(cells)))
} else {
  plot.new(); title("All (no bkgd)"); mtext("No counts", 3, 0.5)
}

labs <- sort(unique(nobkgd$Experiment))
plot.new(); title("Per-experiment summaries"); legend("center", legend = labs, bty = "n")
dev.off()

