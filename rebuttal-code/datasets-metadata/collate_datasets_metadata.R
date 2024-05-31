library(tidyverse)

file_names <- commandArgs(trailingOnly = TRUE)

sample_ids <- strsplit(file_names, split="/")
sample_ids <- unlist(lapply(sample_ids, "[", 3))

all_metadata_fields <- c()

for (i in 1:length(file_names)){
  fname <- file_names[[i]]
  sample <- sample_ids[[i]]

  print(sprintf("---------%s-%s----------", sample, i))


  # read in the metadata file
  metadata <- tryCatch({
    read.table(fname, sep="\t")
  }, error=function(e){
    read.table(fname, sep="\t", fill=TRUE)
  })


  colnames(metadata) <- metadata[1,]
  metadata <- metadata[-1,]

  all_metadata_fields <- c(all_metadata_fields, colnames(metadata))

  # # grab the relevant metadata fields
  # organism_part <- metadata[grepl("organism part", colnames(metadata))]
  # single_cell_isolation <- metadata[grepl("single cell isolation", colnames(metadata))]
  # instrument <- metadata[grepl("INSTRUMENT_MODEL", colnames(metadata))]
  # disease <-  metadata[grepl("disease", colnames(metadata))]
  # assay_name <-  metadata[grepl("Assay Name", colnames(metadata))]
  # tech_type <- metadata[grepl("Technology Type", colnames(metadata))]

  # print(unique(organism_part))
  # print(unique(single_cell_isolation))
  # print(unique(instrument))
  # print(unique(disease))
  # print(unique(assay_name))
  # print(unique(tech_type))
  #
  # organism_part <- paste(organism_part, collapse=", ")
  # single_cell_isolation <- paste(single_cell_isolation, collpase=", ")
  # instrument <- paste(instrument, collapse=", ")
  # disease <- paste(disease, collapse=", ")
  # assay_name <- paste(assay_name, collapse=", ")
  # tech_type <- paste(tech_type, collapse=", ")
}

all_metadata_fields <- unique(unlist(all_metadata_fields))
# take out the really long and annoying fields

remove <- c("technical replicate group", "FASTQ_URI", "Assay Name",
  "ENA_RUN", "Source Name", "read2 file", "Scan Name","Extract Name",
  "read1 file", "ENA_EXPERIMENT", "BioSD_SAMPLE", "Sample_title", "Sample_source_name", "ENA_SAMPLE", "RUN",
  "GSA_EXPERIMENT", "GSA_SAMPLE", "single cell identifier", "individual", "cell type_orig",
  "HCA", "file", "seq_id", "replicate", "Array Data File",
  "single cell library construction", "BAM_URI", "SUBMITTED_FILE_NAME", "Performer")

inds_remove <- lapply(remove, function(rem, list){
  grep(rem, list)
}, all_metadata_fields)
inds_remove <- unlist(inds_remove)
all_metadata_fields <- all_metadata_fields[-inds_remove]

print(all_metadata_fields)

all_metadata_per_dataset <- list()
for (i in 1:length(file_names)){
  fname <- file_names[[i]]
  sample <- sample_ids[[i]]

  print(sprintf("---------%s-%s----------", sample, i))


  # read in the metadata file
  metadata <- tryCatch({
    read.table(fname, sep="\t")
  }, error=function(e){
    read.table(fname, sep="\t", fill=TRUE)
  })

  colnames(metadata) <- metadata[1,]
  metadata <- metadata[-1,]

  metadata <- metadata[,unique(colnames(metadata))]

  print(head(metadata))

  metadata_summary <- lapply(all_metadata_fields, function(meta_name, meta_df){
    meta_name <- gsub("\\[", "\\\\[", meta_name)
    meta_name <- gsub("\\]", "\\\\]", meta_name)
    meta_col <- meta_df[,grepl(meta_name, colnames(meta_df))]
    meta_col <- unique(meta_col)
    print(meta_name)
    print(head(meta_col))
    meta_col <- paste(meta_col, collapse=", ")
  }, metadata)

  names(metadata_summary) <- all_metadata_fields
  #print(metadata_summary)

  metadata_summary_df <- do.call(cbind, metadata_summary)
  colnames(metadata_summary_df) <- all_metadata_fields

  all_metadata_per_dataset[[i]] <- metadata_summary_df
  names(all_metadata_per_dataset)[[i]] <- sample
}

collated_meta_df <- do.call(rbind, all_metadata_per_dataset)
rownames(collated_meta_df) <- names(all_metadata_per_dataset)
print(collated_meta_df)
write.csv(collated_meta_df, "collated_metadata.csv")
