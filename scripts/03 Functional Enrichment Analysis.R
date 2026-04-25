library(gprofiler2)
library(here)
library(dplyr)

# --- SETUP ---
# Path for input (from previous scripts) and output
input_dir <- here("results")
output_dir <- here("results", "enrichment")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define the gene lists you want to analyze
# You can add more files here (e.g., "up_list.txt", "down_list.txt")
input_files <- c("degs_list.txt", "dmgs_list.txt", "cor_gene_list.txt")

# Define the GO ontologies to test
ontologies <- c("GO:BP", "GO:MF", "GO:CC")

# --- EXECUTION LOOP ---
for (file_name in input_files) {
  
  file_path <- file.path(input_dir, file_name)
  
  # Check if the file exists before proceeding
  if (!file.exists(file_path)) {
    message(paste("Skipping:", file_name, "- File not found."))
    next
  }
  
  # Load gene list
  message(paste("Processing enrichment for:", file_name))
  gene_list <- readLines(file_path)
  
  # Clean list name for the output file naming (e.g., "degs_list.txt" -> "degs")
  list_label <- sub("_list\\.txt$", "", file_name)
  
  for (onto in ontologies) {
    
    # Run g:Profiler enrichment
    gost_res <- gost(
      query = gene_list,
      organism = "hsapiens",
      sources = onto,
      correction_method = "fdr",
      significant = FALSE # Show all results to allow manual filtering
    )
    
    # Check if results were found
    if (is.null(gost_res) || is.null(gost_res$result)) {
      message(paste("  No results for", onto))
      next
    }
    
    # Extract and format results
    results <- gost_res$result
    results <- results[order(results$p_value), ]
    
    # Select top 20 terms with p-value < 0.05
    top20_significant <- results %>%
      filter(p_value < 0.05) %>%
      head(20)
    
    # Handle list-type columns (like 'parents') for flat file saving (TSV)
    top20_clean <- as.data.frame(lapply(top20_significant, function(x) {
      if (is.list(x)) sapply(x, paste, collapse = ",") else x
    }))
    
    # Save results if significant terms exist
    if (nrow(top20_clean) > 0) {
      onto_label <- sub("GO:", "", onto) # e.g., GO:BP -> BP
      output_filename <- paste0(list_label, "_", onto_label, "_results.tsv")
      
      write.table(
        top20_clean, 
        file.path(output_dir, output_filename), 
        sep = "\t", 
        row.names = FALSE, 
        col.names = TRUE, 
        quote = FALSE
      )
      message(paste("  Saved results for", onto))
    }
  }
}

message("Functional enrichment analysis complete.")