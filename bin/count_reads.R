#!/usr/bin/env Rscript
library(optparse)
library(Rsubread)
library(ape)
library(stringr)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(reshape2)
library(tibble)


option_list <- list(
    make_option(c("-m", "--metadata"),
        type = "character", default = NULL,
        help = "sample metadata tsv file", metavar = "character"
    ),
    make_option(c("-g", "--gff"),
        type = "character", default = NULL,
        help = "GFF annotation file for the reference strain",
        metavar = "character"
    ),
    make_option(c("-p", "--is_paired"),
        type = "character", default = NULL,
        help = "are the reads paired-end? default = FALSE",
        metavar = "character"
    ),
    make_option(c("-s", "--strandedness"),
        type = "character", default = NULL,
        help = "read strandedness. default = reverse",
        metavar = "character"
    ),
    make_option(c("-t", "--threads"),
        type = "numeric", default = 1,
        help = "number of threads to use. default = 1.",
        metavar = "numeric"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

meta_f <- opt$metadata
gff_f <- opt$gff
ispaired <- if (opt$is_paired == "TRUE") TRUE else FALSE
strandedness <- opt$strandedness
threads <- opt$threads



## ------------------------------------------------------------------------------
## Read data
## ------------------------------------------------------------------------------
meta_tab <- read.table(meta_f, header = TRUE, sep = "\t")

total_counts_list <- lapply(meta_tab$sample, function(x) {
    mapped_count <- read.table(paste0(x, ".counts"), header = FALSE)
    colnames(mapped_count) <- "mapped"
    mapped_count$sample <- x
    mapped_count
})
merged_total_counts <- as.data.frame(do.call(rbind, total_counts_list))



## ------------------------------------------------------------------------------
## Read genome annotation
## ------------------------------------------------------------------------------
ref_annot <- ape::read.gff(gff_f, na.strings = c(".", "?"), GFF3 = TRUE)

ref_annot <- subset(ref_annot, type %in% c("CDS","tRNA"))

gene_biotypes <- ref_annot$type

gene_attr <- stringr::str_split(ref_annot$attributes, ";")

locus_tags <- unlist(lapply(gene_attr, function(x) {
  x[grepl("^locus_tag", x)]
}))

common_gene_names <- unlist(lapply(gene_attr, function(x) {
  x <- x[grepl("^Name=", x)]
  x[identical(x, character(0))] <- ""
  x
}))

gene_lengths <- (ref_annot$end - ref_annot$start) + 1

ref_gene_df <- data.frame(
  locus_tag = locus_tags,
  biotype = gene_biotypes,
  gene_name = common_gene_names,
  gene_length = gene_lengths
)

ref_gene_df$locus_tag <- sub("^.*=", "", ref_gene_df$locus_tag)
ref_gene_df$biotype <- sub("CDS", "protein_coding", ref_gene_df$biotype)
ref_gene_df$gene_name <- sub("^.*=", "", ref_gene_df$gene_name)


write.table(
  ref_gene_df, "ref_gene_df.tsv",
  col.names = TRUE, row.names = FALSE,
  sep = "\t", quote = FALSE
)


## ------------------------------------------------------------------------------
## Count reads mapping to genes
## ------------------------------------------------------------------------------
bamfilesCount <- paste0(meta_tab$sample, ".bam")

# 0 (unstranded), 1 (stranded) and 2 (reversely stranded)
strand <- switch(as.character(strandedness),
    "unstranded" = 0,
    "forward" = 1,
    "reverse" = 2,
    stop("Invalid input")
)

gene_counts <- Rsubread::featureCounts(
    bamfilesCount,
    annot.ext = gff_f,
    isGTFAnnotationFile = TRUE,
    GTF.featureType = "gene",
    GTF.attrType = "locus_tag",
    nthreads = threads,
    countMultiMappingReads = TRUE,
    fraction = TRUE, ## assign fractional counts to multimappers
    isPairedEnd = ispaired,
    strandSpecific = strand
)
colnames(gene_counts$counts) <- gsub(".bam", "", colnames(gene_counts$counts))
colnames(gene_counts$counts) <- gsub("\\.", "_", colnames(gene_counts$counts))


counts_mat <- gene_counts$counts
counts_mat <- tibble::rownames_to_column(as.data.frame(counts_mat), "feature_id")


write.table(
    counts_mat, "gene_counts.tsv",
    col.names = TRUE, row.names = FALSE,
    sep = "\t", quote = FALSE
)

## protein-coding genes only
gene_counts_pc <- counts_mat[ref_gene_df$biotype == "protein_coding", ]
write.table(
    gene_counts_pc, "gene_counts_pc.tsv",
    col.names = TRUE, row.names = FALSE,
    sep = "\t", quote = FALSE
)
