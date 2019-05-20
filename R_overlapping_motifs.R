get_overlap_nd <- function(data, col) {
  real = data[data$sim_id == "real",]
  sims = data[data$sim_id != "real",]
  nd = (real[[col]] - mean(sims[[col]])) / mean(sims[[col]])
  return(nd)
}

get_overlap_p <- function(data, col, tail = "l") {
  real = data[data$sim_id == "real",]
  sims = data[data$sim_id != "real",]
  
  if(tail == "g") {
    p = (nrow(sims[sims[[col]] >= real[[col]],]) + 1) / (nrow(sims) + 1)
  } else {
    p = (nrow(sims[sims[[col]] <= real[[col]],]) + 1) / (nrow(sims) + 1)
  }
  return(p)
}

get_stats <- function(filepath, title) {
  data = read.csv(filepath, head = T)
  output_data = c(
      title,
      data$non_overlap_stop_density[data$sim_id == "real"],
      get_overlap_nd(data, "non_overlap_stop_density"),
      get_overlap_p(data, "non_overlap_stop_density", tail = "g"),
      data$overlap_stop_density[data$sim_id == "real"],
      get_overlap_nd(data, "overlap_stop_density"),
      get_overlap_p(data, "overlap_stop_density", tail = "l")
    )
  return(output_data)
}

output_file = "clean_run/tests/motif_overlaps/motif_overlaps_outputs.csv"

rownames = c("", "Non-overlap density:", "Non-overlap ND:", "Non-overlap P (tail-greater):", "Overlap density:", "Overlap ND:", "Overlap P (tail-less):")
data1 = get_stats("clean_run/tests/motif_overlaps/protein_coding_int3_motif_overlap_density.csv", "Protein Coding INT3")
data2 = get_stats("clean_run/tests/motif_overlaps/protein_coding_RESCUE_motif_overlap_density.csv", "Protein Coding RESCUE")
data3 = get_stats("clean_run/tests/motif_overlaps/lincrna_int3_motif_overlap_density.csv", "LincRNA INT3")
data4 = get_stats("clean_run/tests/motif_overlaps/lincrna_RESCUE_motif_overlap_density.csv", "LincRNA RESCUE")
dataframe = data.frame(rownames, data1, data2, data3, data4)
write.table(dataframe, file = output_file,row.names=FALSE, col.names = FALSE, sep = ",")

