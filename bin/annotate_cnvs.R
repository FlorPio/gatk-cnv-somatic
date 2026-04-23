#!/usr/bin/env Rscript

# ============================================================
# SOMATIC CNV ANALYSIS WITH GENE ANNOTATION
# ============================================================
# Description: Filters CNVs/LOH and annotates cancer-relevant genes
#              using MANE Select transcripts and an external gene list.

library(optparse)
library(dplyr)
library(tidyr)
library(GenomicRanges)

# ============================================================
# PARAMETERS
# ============================================================
option_list <- list(
  make_option(c("--input_dir"), type="character", default="model_segments_results",
              help="Directory containing .modelFinal.seg files"),
  make_option(c("--mane"), type="character", default=NULL,
              help="MANE Select exon annotation file (tab-delimited)"),
  make_option(c("--genes_list"), type="character", default=NULL,
              help="File with gene symbols to annotate (one per line, # for comments)"),
  make_option(c("--output_dir"), type="character", default="CNV_annotated_results",
              help="Output directory"),
  make_option(c("--min_size_kb"), type="numeric", default=50,
              help="Minimum segment size in kb [default: %default]"),
  make_option(c("--maf_threshold"), type="numeric", default=0.35,
              help="LOH threshold (MinorAF < threshold) [default: %default]"),
  make_option(c("--log2cr_amp"), type="numeric", default=0.3,
              help="Amplification threshold (Log2CR > threshold) [default: %default]"),
  make_option(c("--log2cr_del"), type="numeric", default=-0.3,
              help="Deletion threshold (Log2CR < threshold) [default: %default]"),
  make_option(c("--diploid_cn"), type="numeric", default=2,
              help="Normal diploid copy number [default: %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

# ============================================================
# LOAD GENE LIST FROM EXTERNAL FILE
# ============================================================
if (!is.null(opt$genes_list) && file.exists(opt$genes_list)) {
  cat("Loading gene list from:", opt$genes_list, "\n")
  genes_raw <- readLines(opt$genes_list)
  # Remove comments and empty lines
  genes_list <- trimws(genes_raw[!grepl("^#", genes_raw) & nchar(trimws(genes_raw)) > 0])
  cat("  Genes loaded:", length(genes_list), "\n")
} else {
  stop("ERROR: --genes_list is required. Provide a text file with one gene symbol per line.")
}

# Create output directory
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# LOAD DATA
# ============================================================
cat("================================================\n")
cat("SOMATIC CNV ANALYSIS WITH GENE ANNOTATION\n")
cat("================================================\n\n")

# Find .modelFinal.seg files
input_files <- list.files(opt$input_dir,
                          pattern = ".*\\.modelFinal\\.seg$",
                          full.names = TRUE)

if (length(input_files) == 0) {
  stop("No .modelFinal.seg files found in: ", opt$input_dir)
}

cat("Files found:\n")
for (f in input_files) {
  cat("  -", basename(f), "\n")
}
cat("\n")

# Load MANE Select annotation
if (is.null(opt$mane) || !file.exists(opt$mane)) {
  stop("ERROR: --mane annotation file is required and must exist.")
}

cat("Loading MANE annotations...\n")
MANE_file <- read.delim(opt$mane, sep = "\t", header = TRUE)
MANE_G <- makeGRangesFromDataFrame(
  MANE_file,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end",
  keep.extra.columns = TRUE
)

# ============================================================
# EVENT CLASSIFICATION FUNCTION (based on absolute CN)
# ============================================================
classify_event_by_cn <- function(cn_absolute, minor_af,
                                  diploid_cn = 2, maf_th = 0.35) {
  event <- case_when(
    cn_absolute < diploid_cn & minor_af < maf_th ~ "DEL+LOH",
    cn_absolute < diploid_cn & minor_af >= maf_th ~ "DELETION",
    cn_absolute > diploid_cn & minor_af < maf_th ~ "AMP+LOH",
    cn_absolute > diploid_cn & minor_af >= maf_th ~ "AMPLIFICATION",
    cn_absolute == diploid_cn & minor_af < maf_th ~ "CN-LOH",
    cn_absolute == diploid_cn & minor_af >= maf_th ~ "NEUTRAL",
    TRUE ~ "UNKNOWN"
  )
  return(event)
}

# ============================================================
# TUMOR PURITY ESTIMATION FUNCTION
# ============================================================
# Based on Minor Allele Fraction (MAF) of LOH segments
# Purity = 1 - (2 x weighted_avg_MAF)
estimate_purity <- function(seg_data, min_segments = 1) {

  seg_classified <- seg_data %>%
    mutate(
      LOG2CR = LOG2_COPY_RATIO_POSTERIOR_50,
      MINOR_AF = MINOR_ALLELE_FRACTION_POSTERIOR_50,
      SIZE_KB = (END - START) / 1000
    )

  seg_classified <- seg_classified %>%
    mutate(
      HAS_LOH = case_when(
        is.na(MINOR_AF) ~ FALSE,
        MINOR_AF < 0.35 ~ TRUE,
        TRUE ~ FALSE
      )
    )

  loh_segments <- seg_classified %>%
    filter(
      HAS_LOH == TRUE,
      !is.na(MINOR_AF),
      MINOR_AF > 0,
      SIZE_KB >= 50
    )

  if (nrow(loh_segments) < min_segments) {
    warning("No LOH segments detected for purity estimation (n=", nrow(loh_segments), ")")
    return(list(
      purity = NA,
      n_segments = 0,
      avg_maf = NA,
      weighted_avg_maf = NA,
      method = "no_LOH_detected"
    ))
  }

  loh_segments <- loh_segments %>%
    mutate(WEIGHT = SIZE_KB)

  weighted_avg_maf <- weighted.mean(
    loh_segments$MINOR_AF,
    loh_segments$WEIGHT
  )

  avg_maf <- mean(loh_segments$MINOR_AF)

  # Formula: Purity = 1 - (2 x avg_MAF)
  # Normal heterozygous tissue: MAF = 0.5
  # Pure tumor with LOH: MAF -> 0
  # Mixture: MAF = (1-p)*0.5 + p*0 = 0.5*(1-p), so p = 1 - 2*MAF
  purity_estimate <- 1 - (2 * weighted_avg_maf)
  purity_estimate <- max(0.1, min(1.0, purity_estimate))

  return(list(
    purity = round(purity_estimate, 4),
    n_segments = nrow(loh_segments),
    avg_maf = round(avg_maf, 4),
    weighted_avg_maf = round(weighted_avg_maf, 4),
    method = "LOH_based"
  ))
}

# ============================================================
# ABSOLUTE COPY NUMBER CALCULATION
# ============================================================
calculate_absolute_cn <- function(log2cr, purity, diploid_cn = 2) {

  if (length(purity) == 1) {
    purity <- rep(purity, length(log2cr))
  }

  # Observed CN (tumor + normal mixture)
  ratio_observed <- 2^log2cr
  cn_observed <- diploid_cn * ratio_observed

  # If no valid purity, return observed CN without correction
  cn_tumor <- ifelse(
    is.na(purity) | purity < 0.2 | purity > 1.0,
    cn_observed,
    # Purity-corrected CN:
    (cn_observed - (1 - purity) * diploid_cn) / purity
  )

  cn_absolute <- round(cn_tumor)
  cn_absolute <- pmax(0, pmin(20, cn_absolute))

  return(cn_absolute)
}

# ============================================================
# PROCESS EACH SAMPLE
# ============================================================
all_results <- list()
sample_names <- c()
purity_results <- list()

for (i in seq_along(input_files)) {
  file_path <- input_files[i]
  sample_name <- sub(".modelFinal.seg", "", basename(file_path))
  sample_names <- c(sample_names, sample_name)

  cat("Processing:", sample_name, "\n")

  seg_data <- read.delim(
    file_path,
    sep = "\t",
    header = TRUE,
    comment.char = "@",
    stringsAsFactors = FALSE
  )

  # Estimate tumor purity
  cat("  - Estimating tumor purity...\n")
  purity_info <- estimate_purity(seg_data)
  sample_purity <- purity_info$purity

  # Store purity info
  purity_results[[i]] <- data.frame(
    SAMPLE = sample_name,
    PURITY = ifelse(is.na(sample_purity), NA, sample_purity),
    PURITY_PERCENT = ifelse(is.na(sample_purity), NA, round(sample_purity * 100, 2)),
    N_LOH_SEGMENTS = purity_info$n_segments,
    AVG_MAF = purity_info$avg_maf,
    WEIGHTED_AVG_MAF = purity_info$weighted_avg_maf,
    METHOD = purity_info$method
  )

  if (is.na(sample_purity)) {
    cat("    Warning: Could not estimate purity (no LOH regions detected)\n")
    cat("    CN will be calculated WITHOUT purity correction (observed CN)\n")
  } else {
    cat("    Estimated purity:", round(sample_purity * 100, 2), "%\n")
    cat("    LOH segments analyzed:", purity_info$n_segments, "\n")
  }

  # Filter and classify using purity-corrected absolute CN
  seg_filtered <- seg_data %>%
    mutate(
      SIZE_KB = (END - START) / 1000,
      LOG2CR = LOG2_COPY_RATIO_POSTERIOR_50,
      MINOR_AF = MINOR_ALLELE_FRACTION_POSTERIOR_50,
      PURITY = sample_purity,
      CN_ABSOLUTE = calculate_absolute_cn(LOG2CR, sample_purity, opt$diploid_cn),
      EVENT_TYPE = classify_event_by_cn(CN_ABSOLUTE, MINOR_AF, opt$diploid_cn, opt$maf_threshold),
      SAMPLE = sample_name
    ) %>%
    filter(
      SIZE_KB >= opt$min_size_kb,
      EVENT_TYPE != "NEUTRAL"
    )

  # Convert to GRanges for overlap
  seg_G <- makeGRangesFromDataFrame(
    seg_filtered,
    seqnames.field = "CONTIG",
    start.field = "START",
    end.field = "END",
    keep.extra.columns = TRUE
  )

  # Annotate with MANE
  overlaps <- findOverlaps(seg_G, MANE_G)

  if (length(overlaps) > 0) {
    query_hits <- queryHits(overlaps)
    subject_hits <- subjectHits(overlaps)

    annotated <- seg_filtered[query_hits, ] %>%
      mutate(
        GENE = MANE_file$gene_symbol[subject_hits],
        EXON = MANE_file$exon_number[subject_hits]
      ) %>%
      filter(GENE %in% genes_list) %>%
      distinct() %>%
      arrange(CONTIG, START)

    all_results[[i]] <- annotated
  } else {
    cat("  Warning: No overlaps found with annotated genes\n")
    all_results[[i]] <- data.frame()
  }
}

# ============================================================
# COMBINE RESULTS
# ============================================================
if (length(all_results) > 0) {
  combined_results <- bind_rows(all_results)

  # Save combined results
  output_combined <- file.path(opt$output_dir, "CNVs_annotated_all_samples.txt")
  write.table(
    combined_results,
    file = output_combined,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  cat("\nResults saved to:", output_combined, "\n")

  # Save purity estimates
  purity_table <- bind_rows(purity_results)
  output_purity <- file.path(opt$output_dir, "tumor_purity_estimates.txt")
  write.table(
    purity_table,
    file = output_purity,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  cat("Purity estimates saved to:", output_purity, "\n")

  # Summary by sample
  summary_by_sample <- combined_results %>%
    group_by(SAMPLE, EVENT_TYPE) %>%
    summarise(
      N_EVENTS = n(),
      N_GENES = n_distinct(GENE),
      TOTAL_SIZE_MB = sum(SIZE_KB) / 1000,
      .groups = "drop"
    ) %>%
    arrange(SAMPLE, EVENT_TYPE)

  output_summary <- file.path(opt$output_dir, "summary_by_sample.txt")
  write.table(
    summary_by_sample,
    file = output_summary,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  cat("Sample summary saved to:", output_summary, "\n")

  # Gene x Sample matrix
  gene_matrix <- combined_results %>%
    select(SAMPLE, GENE, EVENT_TYPE) %>%
    distinct() %>%
    pivot_wider(
      names_from = SAMPLE,
      values_from = EVENT_TYPE,
      values_fill = list(EVENT_TYPE = "NORMAL"),
      values_fn = function(x) paste(unique(x), collapse = ";")
    ) %>%
    arrange(GENE)

  output_matrix <- file.path(opt$output_dir, "gene_matrix.txt")
  write.table(
    gene_matrix,
    file = output_matrix,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  cat("Gene x sample matrix saved to:", output_matrix, "\n")

  # Top affected genes
  top_genes <- combined_results %>%
    group_by(GENE) %>%
    summarise(
      N_SAMPLES_AFFECTED = n_distinct(SAMPLE),
      EVENTS = paste(unique(EVENT_TYPE), collapse = ", "),
      AVG_LOG2CR = mean(LOG2CR),
      AVG_MINOR_AF = mean(MINOR_AF),
      .groups = "drop"
    ) %>%
    arrange(desc(N_SAMPLES_AFFECTED), GENE)

  output_top_genes <- file.path(opt$output_dir, "top_affected_genes.txt")
  write.table(
    top_genes,
    file = output_top_genes,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  cat("Top affected genes saved to:", output_top_genes, "\n")

  # Individual sample files
  for (sample in unique(combined_results$SAMPLE)) {
    sample_data <- combined_results %>%
      filter(SAMPLE == sample) %>%
      select(-SAMPLE)

    output_individual <- file.path(opt$output_dir, paste0(sample, "_annotated.txt"))
    write.table(
      sample_data,
      file = output_individual,
      quote = FALSE,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE
    )
  }
  cat("Individual sample files saved\n")

  # Console summary
  cat("\n================================================\n")
  cat("ANALYSIS SUMMARY\n")
  cat("================================================\n\n")

  cat("Samples processed:", length(sample_names), "\n")
  cat("Unique genes affected:", n_distinct(combined_results$GENE), "\n")
  cat("Total annotated events:", nrow(combined_results), "\n\n")

  cat("Tumor purity estimates:\n")
  purity_display <- purity_table %>%
    select(SAMPLE, PURITY_PERCENT, N_LOH_SEGMENTS, WEIGHTED_AVG_MAF, METHOD)
  print(purity_display, row.names = FALSE)
  cat("\n")

  cat("Event distribution:\n")
  print(table(combined_results$EVENT_TYPE))

  cat("\nAbsolute copy number distribution:\n")
  print(table(combined_results$CN_ABSOLUTE))

  cat("\nTop 10 most affected genes:\n")
  print(head(top_genes, 10))

} else {
  cat("\nWarning: No events found to annotate\n")
}

cat("\n================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("================================================\n")
