#!/usr/bin/Rscript --vanilla


### Initiating user input ###
args <- commandArgs(trailingOnly = TRUE)
sam_path <- args[1]
exons_path <- args[2]


### Padding the SAM file into a data frame ###
sam <- sapply(strsplit(readLines(sam_path, warn = FALSE), split = "\t"), FUN = unlist)
col_max <- max(sapply(sam, FUN = length))
sam_pad <- as.data.frame(do.call(what = rbind, args = lapply(sam, function(row){ 
  if(length(row) < col_max){ 
    c(row, rep(NA, length.out = col_max - length(row)))
  }
  else{ 
    row 
  }
})))


### Identifying SAM file reads that start within 5 bp of a W1 or W1' exon start coordinate, aka candidate reads for W0 correction ###
exons <- read.table(exons_path, header = TRUE, sep = "\t")
total_reads <- subset(sam_pad, !grepl(sam_pad$V1, pattern = "^@"))
colnames(total_reads) <- c("qname", "flag", "rname", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq", "qual", paste0("V", 12:(col_max)))
candidate_reads <- total_reads[apply(total_reads, MARGIN = 1, function(read){ 
  any(abs(as.numeric(read["pos"]) - exons$start[exons$exon == "W1" | exons$exon == "W1_prime"]) <= 5) 
}), ]


### Extracting the left soft clip + 10 bp of each candidate reads ###
candidate_left_clip_sizes <- as.numeric(ifelse(grepl(candidate_reads$cigar, pattern = "^(\\d+)S", perl = TRUE), yes = sub(candidate_reads$cigar, pattern = "^(\\d+)S.*", replacement = "\\1"), no = 0)) 
candidate_seqs <- substr(candidate_reads$seq, start = 1, stop = candidate_left_clip_sizes + 10)


### Searching for W0 and W2 motifs in the extracted sequence of each candidate read ###
candidate_ACAAT_or_ACAAAT_statuses <- grepl(substr(candidate_seqs, start = pmax(1, candidate_left_clip_sizes - 11), stop = nchar(candidate_seqs)), pattern = "ACAAT|ACAAAT")
candidate_AGGAGT_or_CAGGAGT_statuses <- grepl(substr(candidate_seqs, start = pmax(1, candidate_left_clip_sizes - 24), stop = candidate_left_clip_sizes), pattern = "AGGAGT") | grepl(substr(candidate_seqs, start = 1, stop = candidate_left_clip_sizes), pattern = "CAGGAGT")
candidate_CAGGG_statuses <- grepl(substr(candidate_seqs, start = pmax(1, candidate_left_clip_sizes - 11), stop = nchar(candidate_seqs)), pattern = "CAGGG")


### Filtering the candidate reads based on the motif search results ###
W0_reads <- candidate_reads[(candidate_ACAAT_or_ACAAAT_statuses | candidate_AGGAGT_or_CAGGAGT_statuses) & !candidate_CAGGG_statuses, ]
W2_reads <- candidate_reads[!candidate_ACAAT_or_ACAAAT_statuses & !candidate_AGGAGT_or_CAGGAGT_statuses & candidate_CAGGG_statuses, ]
mixed_reads <- candidate_reads[(candidate_ACAAT_or_ACAAAT_statuses | candidate_AGGAGT_or_CAGGAGT_statuses) & candidate_CAGGG_statuses, ]
mystery_reads <- candidate_reads[!candidate_ACAAT_or_ACAAAT_statuses & !candidate_AGGAGT_or_CAGGAGT_statuses & !candidate_CAGGG_statuses, ]


### Unmasking, if present, the left soft clip of each candidate read ###
unmasked_candidate_reads <- do.call(what = rbind, args = lapply(seq_len(nrow(candidate_reads)), function(candidate_read){
  if(candidate_left_clip_sizes[candidate_read] > 0){
    candidate_reads[candidate_read,]["cigar"] <- sub(candidate_reads[candidate_read,]["cigar"], pattern = paste0("^", candidate_left_clip_sizes[candidate_read], "S"), replacement = paste0(candidate_left_clip_sizes[candidate_read], "M"))
    double_match_status <- grepl(candidate_reads[candidate_read,]["cigar"], pattern = "^(\\d+)M(\\d+)M", perl = TRUE)
    if(double_match_status == TRUE){
      sum_of_matches <- sum(as.numeric(unlist(regmatches(candidate_reads[candidate_read,]["cigar"], m = gregexpr(regmatches(candidate_reads[candidate_read,]["cigar"], m = regexpr(candidate_reads[candidate_read,]["cigar"], pattern = "^(\\d+)M(\\d+)M", perl = TRUE)), pattern = "\\d+"[[1]])))))
      candidate_reads[candidate_read,]["cigar"] <- sub(candidate_reads[candidate_read,]["cigar"], pattern = "^(\\d+)M(\\d+)M", replacement = paste0(sum_of_matches, "M"))
    }
    candidate_reads[candidate_read,]["pos"] <- as.numeric(candidate_reads[candidate_read,]["pos"]) - candidate_left_clip_sizes[candidate_read]
  }
  candidate_reads[candidate_read,] 
}))


### Attempting to add back W0 to each candidate read identified as containing W0 motif(s) in the previous motif search ###
corrected_candidate_reads <- do.call(what = rbind, args = lapply(seq_len(nrow(unmasked_candidate_reads)), function(unmasked_candidate_read){
  if((candidate_ACAAT_or_ACAAAT_statuses[unmasked_candidate_read] == TRUE | candidate_AGGAGT_or_CAGGAGT_statuses[unmasked_candidate_read]) & candidate_CAGGG_statuses[unmasked_candidate_read] == FALSE){
    motif_end <- -1
    W0_end <- -1
    distance_to_W0 <- -Inf
    if(grepl(substr(candidate_seqs[unmasked_candidate_read], start = pmax(1, candidate_left_clip_sizes[unmasked_candidate_read] - 11), stop = nchar(candidate_seqs[unmasked_candidate_read])), pattern = "ACAAT")){
      W0_end <- regexpr(substr(candidate_seqs[unmasked_candidate_read], pmax(1, start = candidate_left_clip_sizes[unmasked_candidate_read] - 11), stop = nchar(candidate_seqs[unmasked_candidate_read])), pattern = "ACAAT")[1] + 4
      if(candidate_left_clip_sizes[unmasked_candidate_read] - 12 > 0){
        W0_end <- W0_end + (candidate_left_clip_sizes[unmasked_candidate_read] - 12)
      }
    }
    else if(grepl(substr(candidate_seqs[unmasked_candidate_read], start = pmax(1, candidate_left_clip_sizes[unmasked_candidate_read] - 11), stop = nchar(candidate_seqs[unmasked_candidate_read])), pattern = "ACAAAT")){
      W0_end <- regexpr(substr(candidate_seqs[unmasked_candidate_read], start = pmax(1, candidate_left_clip_sizes[unmasked_candidate_read] - 11), stop = nchar(candidate_seqs[unmasked_candidate_read])), pattern = "ACAAAT")[1] + 5
      if(candidate_left_clip_sizes[unmasked_candidate_read] - 12 > 0){
        W0_end <- W0_end + (candidate_left_clip_sizes[unmasked_candidate_read] - 12)
      }
    }
    else if(grepl(substr(candidate_seqs[unmasked_candidate_read], start = pmax(1, candidate_left_clip_sizes[unmasked_candidate_read] - 24), stop = candidate_left_clip_sizes[unmasked_candidate_read]), pattern = "AGGAGT")){
      motif_end <- regexpr(substr(candidate_seqs[unmasked_candidate_read], start = pmax(1, candidate_left_clip_sizes[unmasked_candidate_read] - 24), stop = candidate_left_clip_sizes[unmasked_candidate_read]), pattern = "AGGAGT")[1] + 5
      if(candidate_left_clip_sizes[unmasked_candidate_read] - 25 > 0){
        motif_end <- motif_end + (candidate_left_clip_sizes[unmasked_candidate_read] - 25)
      }
    }
    else if(grepl(substr(candidate_seqs[unmasked_candidate_read], start = 1, stop = candidate_left_clip_sizes[unmasked_candidate_read]), pattern = "CAGGAGT")){
      motif_end <- regexpr(substr(candidate_seqs[unmasked_candidate_read], start = 1, stop = candidate_left_clip_sizes[unmasked_candidate_read]), pattern = "CAGGAGT")[1] + 6
    }
    if(motif_end != -1 & motif_end + 4 < nchar(candidate_seqs[unmasked_candidate_read])){
      W0_end <- regexpr(substr(candidate_seqs[unmasked_candidate_read], start = motif_end + 5, stop = nchar(candidate_seqs[unmasked_candidate_read])), pattern = "AT")[1] + 1
      if(W0_end > 0){
        W0_end <- W0_end + motif_end + 4
      }
    }
    if(W0_end > 0){
      distance_to_W0 <- (as.numeric(unmasked_candidate_reads[unmasked_candidate_read,]["pos"]) + W0_end - 1) - max(exons$end[exons$exon == "W0" & exons$end < (as.numeric(unmasked_candidate_reads[unmasked_candidate_read,]["pos"]) + W0_end - 1)], default = -Inf)
    }
    if(distance_to_W0 != -Inf){
      match_size <- as.numeric(sub(unmasked_candidate_reads[unmasked_candidate_read,]["cigar"], pattern = "^(\\d+)M.*", replacement = "\\1"))
      cigar_minus_match <- substr(unmasked_candidate_reads[unmasked_candidate_read, "cigar"], start = regexpr("[^0-9]", unmasked_candidate_reads[unmasked_candidate_read, "cigar"])[1] + 1, stop = nchar(unmasked_candidate_reads[unmasked_candidate_read, "cigar"]))
      if(W0_end < match_size){
        unmasked_candidate_reads[unmasked_candidate_read, ]["cigar"] <- paste0(W0_end, "M", distance_to_W0, "N", match_size - W0_end, "M", cigar_minus_match)
        unmasked_candidate_reads[unmasked_candidate_read, ]["pos"] <- as.numeric(unmasked_candidate_reads[unmasked_candidate_read, ]["pos"]) - distance_to_W0
      }
      else{
        unmasked_candidate_reads[unmasked_candidate_read, ]["cigar"] <- paste0(match_size, "M", distance_to_W0 - (W0_end - match_size), "N", cigar_minus_match)
        unmasked_candidate_reads[unmasked_candidate_read, ]["pos"] <- as.numeric(unmasked_candidate_reads[unmasked_candidate_read, ]["pos"]) - distance_to_W0 + (W0_end - match_size)
      }
    }
  }
  unmasked_candidate_reads[unmasked_candidate_read,]
}))


### Updating the original list of total reads with the abovementioned unmasking and correcting ### 
total_reads[!is.na(match(total_reads$qname, table = corrected_candidate_reads$qname)), ] <- corrected_candidate_reads[na.omit(match(total_reads$qname, table = corrected_candidate_reads$qname)), ]
corrected_total_reads <- total_reads


### Preparing the updated total reads and (un-updated) filtered candidate reads for outputting as individual SAM files ###
output_df_list <- list(corrected_total_reads[order(as.numeric(corrected_total_reads$pos)), ], W0_reads, W2_reads, mixed_reads, mystery_reads)
output_df_list <- lapply(output_df_list, function(output_df){ 
  colnames(output_df) <- paste0("V", 1:ncol(output_df)); output_df 
})
header <- subset(sam_pad, grepl(sam_pad$V1, pattern = "^@"))
output_sam_pad_list <- lapply(output_df_list, function(output_df){ 
  rbind(header, output_df) 
})


### Outputting the updated total reads and (un-updated) filtered candidate reads as individual SAM files ###
output_sam_names <- c(paste0("updated_", basename(sam_path)), "candidates_W0_id.sam", "candidates_W2_id.sam", "candidates_mixed_id.sam", "candidates_unassigned_id.sam")
output_sam_list <- lapply(seq_along(output_sam_pad_list), function(output_sam_pad){ 
  output_sam <- apply(output_sam_pad_list[[output_sam_pad]], MARGIN = 1, function(row){ 
    unpadded_row_string <- paste(row[!is.na(row)], collapse = "\t") 
  }); write.table(output_sam, file = output_sam_names[output_sam_pad], sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE) 
})
