#!/usr/bin/Rscript --vanilla

suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(stringi)))
suppressMessages(suppressWarnings(library(tidyverse)))


### Initiate user input ###
args <- commandArgs(trailingOnly = TRUE)
sam_path <- args[1]
significance_window <- as.numeric(args[2])
mininum_significant_read_count <- as.numeric(args[3])
noise_window <- as.numeric(args[4])
clip_window <- as.numeric(args[5])


### Open and store the SAM file ###
sam_conn <- file(sam_path, open = "r")
sam_file <- list()
for(i in readLines(sam_conn, warn = FALSE)){
  split_i <- unlist(strsplit(i, split = "\t"))
  sam_file <- c(sam_file, list(split_i))
}
close(sam_conn)


### Calculate the maximum column number of the SAM file ###
col_max <- 0
for(i in 1:length(sam_file)){
  if(length(sam_file[[i]]) > col_max){
    col_max <- length(sam_file[[i]])
  }
}


### Pad the SAM file rows with NAs based on the maximum column number ###
sam_pad <- data.frame(matrix(nrow = 0, ncol = col_max))
for(i in 1:length(sam_file)){
  if(length(sam_file[[i]]) < col_max){
    sam_file[[i]] <- c(sam_file[[i]], rep(NA, length.out = col_max - length(sam_file[[i]])))
  }
  sam_pad <- rbind(sam_pad, sam_file[[i]])
}
colnames(sam_pad) <- c("qname", "flag", "rname", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq", "qual", paste0("V", 12:(col_max)))
sam_header <- subset(sam_pad, grepl(sam_pad$qname, pattern = "^@"))
sam_reads <- subset(sam_pad, !grepl(sam_pad$qname, pattern = "^@"))


### Extract and store the 5' end position for plus strand reads and the 3' end position for minus strand reads ###
position_data <- data.frame(matrix(nrow = 0, ncol = 5))
for(i in 1:nrow(sam_reads)){
    position_data <- rbind(position_data, c(sam_reads$qname[i], sam_reads$flag[i], sam_reads$pos[i], sam_reads$cigar[i], as.character(sam_reads$seq[i])))
}
colnames(position_data) <- c("qname", "dna_strand", "bam_read_position", "read_cigar", "read_sequence")
position_data$dna_strand[position_data$dna_strand == 0] <- "+"
position_data$dna_strand[position_data$dna_strand == 16] <- "-"


### Extract and store the 5' end position for minus strand reads and the 3' end position for plus strand reads ###
reverse_bam_read_position <- c()
for(i in 1:nrow(position_data)){
    length_on_genome <- 0
    cigar <- unlist(strsplit(position_data$read_cigar[i], split = ""))
    for(j in 1:length(cigar)){
        if(cigar[j] == "N" | cigar[j] == "M" | cigar[j] == "D"){
            reverse_event_length <- ""
            event_length <- 0
            for(k in (j-1):1){
                if(cigar[k] != "S" & cigar[k] != "H" & cigar[k] != "M" & cigar[k] != "I" & cigar[k] != "D" & cigar[k] != "N"){
                    reverse_event_length <- paste(reverse_event_length, cigar[k], sep = "")
                }
                else{
                    break
                }
            }
            event_length <- as.numeric(stri_reverse(reverse_event_length))
            length_on_genome <- length_on_genome + event_length
        }
    }
    new_read_position <- (length_on_genome - 1) + as.numeric(position_data$bam_read_position[i])
    reverse_bam_read_position <- c(reverse_bam_read_position, new_read_position)
}
position_data <- add_column(position_data, reverse_bam_read_position, .after = 3)


### Separate plus strand and minus strand reads then separate into 5' ends and ends 3' ends ###
plus_position_data <- subset(position_data, position_data$dna_strand == "+", select = -dna_strand)
colnames(plus_position_data) <- c("qname", "read_start", "read_end", "read_cigar", "read_sequence")
minus_position_data <- subset(position_data, position_data$dna_strand == "-", select = -dna_strand)
colnames(minus_position_data) <- c("qname", "read_end", "read_start", "read_cigar", "read_sequence")
minus_position_data <- minus_position_data[, colnames(plus_position_data)]
combined_position_data <- list(plus_start_position_data = subset(plus_position_data, select = -read_end), plus_end_position_data = subset(plus_position_data, select = -read_start), minus_start_position_data = subset(minus_position_data, select = -read_end), minus_end_position_data = subset(minus_position_data, select = -read_start))
combined_position_data[[1]] <- combined_position_data[[1]][order(as.numeric(combined_position_data[[1]]$read_start)), ]
combined_position_data[[2]] <- combined_position_data[[2]][order(as.numeric(combined_position_data[[2]]$read_end)), ]
combined_position_data[[3]] <- combined_position_data[[3]][order(as.numeric(combined_position_data[[3]]$read_start)), ]
combined_position_data[[4]] <- combined_position_data[[4]][order(as.numeric(combined_position_data[[4]]$read_end)), ]


### Define significant (reproducible) 5' ends and 3' ends based on significance_window and mininum_significant_read_count input ### 
combined_position_comparison_data <- vector(mode = "list", length = 4)
unassigned_qnames <- list(plus_start = position_data$qname, plus_end = position_data$qname, minus_start = position_data$qname, minus_end = position_data$qname)
for(i in 1:length(combined_position_data)){
    comparison_data <- data.frame(matrix(nrow = 0, ncol = 3))
    last_query_index <- 0
    for(j in 1:nrow(combined_position_data[[i]])){
        significant_read_count <- 0
        significant_read_qnames <- c()
        if(j == last_query_index + 1){
            reference_position <- as.numeric(combined_position_data[[i]][j, 2])
            for(k in 1:nrow(combined_position_data[[i]])){
                query_position <- as.numeric(combined_position_data[[i]][k, 2])
                if((query_position >= reference_position - significance_window) & (query_position <= reference_position + significance_window)){
                    significant_read_count <- significant_read_count + 1
                    significant_read_qnames <- paste(significant_read_qnames, combined_position_data[[i]]$qname[k], sep = ",")
                    last_query_index <- k
                }
            }
            if(significant_read_count >= mininum_significant_read_count){
                comparison_data <- rbind(comparison_data, c(reference_position, significant_read_count, significant_read_qnames))
                significant_read_qnames <- unlist(strsplit(significant_read_qnames, split = ","))
                unassigned_qnames[[i]] <- unassigned_qnames[[i]][!(unassigned_qnames[[i]] %in% significant_read_qnames)]
            }
        }
    }
    comparison_data <- cbind(comparison_data, noisy_read_count = 0, noisy_read_qnames = "", clipped_read_count = 0, clipped_read_qnames = "")
    colnames(comparison_data) <- c("genomic_position", "significant_read_count", "significant_read_qnames", "noisy_read_count", "noisy_read_qnames", "clipped_read_count", "clipped_read_qnames")
    comparison_data <- comparison_data[order(as.numeric(comparison_data$genomic_position)), ]
    combined_position_comparison_data[[i]] <- comparison_data
}


### Redefinine 5' and 3' ends mistaken as noise in the previous stage as significant based on the noise_window input ###
for(i in 1:length(combined_position_data)){
    assigned_noisy_qnames <- c()
    for(j in 1:nrow(combined_position_data[[i]])){
        for(k in 1:length(unassigned_qnames[[i]])){
            if(unassigned_qnames[[i]][k] == combined_position_data[[i]]$qname[j]){ 
                for(l in 1:nrow(combined_position_comparison_data[[i]])){
                    if((as.numeric(combined_position_data[[i]][j, 2]) >= as.numeric(combined_position_comparison_data[[i]]$genomic_position[l]) - noise_window) & (as.numeric(combined_position_data[[i]][j, 2]) <= as.numeric(combined_position_comparison_data[[i]]$genomic_position[l]) + noise_window)){
                        combined_position_comparison_data[[i]]$noisy_read_count[l] <- as.numeric(combined_position_comparison_data[[i]]$noisy_read_count[l]) + 1
                        combined_position_comparison_data[[i]]$noisy_read_qnames[l] <- paste(combined_position_comparison_data[[i]]$noisy_read_qnames[l], unassigned_qnames[[i]][k], sep = ",")
                        assigned_noisy_qnames <- c(assigned_noisy_qnames, unassigned_qnames[[i]][k])
                    }
                }
            }
        }
    }
    if(length(assigned_noisy_qnames) > 0){
        unassigned_qnames[[i]] <- unassigned_qnames[[i]][!(unassigned_qnames[[i]] %in% assigned_noisy_qnames)]
    }
    combined_position_comparison_data[[i]]$noisy_read_qnames[combined_position_comparison_data[[i]]$noisy_read_qnames == ""] <- "NA"
    combined_position_comparison_data[[i]] <- combined_position_comparison_data[[i]][order(as.numeric(combined_position_comparison_data[[i]]$genomic_position)), ]
}


### Redefinine 5' and 3' ends mistaken as noise in the previous stages as significant based on the clip_window input ###
for(i in 1:length(combined_position_data)){
    assigned_clipped_qnames <- c()
    for(j in 1:nrow(combined_position_data[[i]])){
        for(k in 1:length(unassigned_qnames[[i]])){
            if(unassigned_qnames[[i]][k] == combined_position_data[[i]]$qname[j]){
                reverse_left_clip <- ""
                reverse_right_clip <- ""
                left_S_index = 0
                right_S_index = 0
                left_clip_length <- 0
                right_clip_length <- 0
                new_start_position <- 0
                new_end_position <- 0
                next_clip_is_right <- FALSE
                cigar <- unlist(strsplit(combined_position_data[[i]]$read_cigar[j], split = ""))
                for(l in 1:length(cigar)){
                    if(i == 1 | i == 4){
                        if(cigar[l] == "S"){
                            left_S_index <- l 
                            for(m in (l - 1):1){
                                reverse_left_clip <- paste(reverse_left_clip, cigar[m], sep = "")
                            }
                            left_clip_length <- as.numeric(stri_reverse(reverse_left_clip))
                            if(i == 1){
                                new_start_position <- as.numeric(combined_position_data[[i]][j, 2]) - left_clip_length
                            }
                            else if(i == 4){
                                new_end_position <- as.numeric(combined_position_data[[i]][j, 2]) - left_clip_length
                            }
                            for(m in 1:nrow(combined_position_comparison_data[[i]])){
                                if(i == 1){
                                    if((new_start_position >= as.numeric(combined_position_comparison_data[[i]]$genomic_position[m]) - clip_window) & (new_start_position <= as.numeric(combined_position_comparison_data[[i]]$genomic_position[m]) + clip_window)){
                                        combined_position_comparison_data[[i]]$clipped_read_count[m] <- as.numeric(combined_position_comparison_data[[i]]$clipped_read_count[m]) + 1
                                        combined_position_comparison_data[[i]]$clipped_read_qnames[m] <- paste(combined_position_comparison_data[[i]]$clipped_read_qnames[m], combined_position_data[[i]]$qname[j], sep = ",")
                                        assigned_clipped_qnames <- c(assigned_clipped_qnames, unassigned_qnames[[i]][k])
                                        position_data$bam_read_position[position_data$qname == combined_position_data[[i]]$qname[j]] <- as.character(new_start_position)
                                        position_data$read_cigar[position_data$qname == combined_position_data[[i]]$qname[j]] <- paste(as.character(left_clip_length), "M", substring(combined_position_data[[i]]$read_cigar[j], (left_S_index + 1), length(cigar)), sep = "")
                                    }
                                }
                                else if(i == 4){
                                    if((new_end_position >= as.numeric(combined_position_comparison_data[[i]]$genomic_position[m]) - clip_window) & (new_end_position <= as.numeric(combined_position_comparison_data[[i]]$genomic_position[m]) + clip_window)){
                                        combined_position_comparison_data[[i]]$clipped_read_count[m] <- as.numeric(combined_position_comparison_data[[i]]$clipped_read_count[m]) + 1
                                        combined_position_comparison_data[[i]]$clipped_read_qnames[m] <- paste(combined_position_comparison_data[[i]]$clipped_read_qnames[m], combined_position_data[[i]]$qname[j], sep = ",")
                                        assigned_clipped_qnames <- c(assigned_clipped_qnames, unassigned_qnames[[i]][k])
                                        position_data$bam_read_position[position_data$qname == combined_position_data[[i]]$qname[j]] <- as.character(new_end_position)
                                        position_data$read_cigar[position_data$qname == combined_position_data[[i]]$qname[j]] <- paste(as.character(left_clip_length), "M", substring(combined_position_data[[i]]$read_cigar[j], (left_S_index + 1), length(cigar)), sep = "")
                                    }
                                }
                            }
                            break
                        }
                        else if(cigar[l] == "M" | cigar[l] == "I" | cigar[l] == "D" | cigar[l] == "N" | cigar[l] == "H"){
                            break
                        }
                    }
                    else if(i == 2 | i == 3){
                        if(cigar[l] == "M" | cigar[l] == "I" | cigar[l] == "D" | cigar[l] == "N" | cigar[l] == "H"){
                            next_clip_is_right <- TRUE
                        }
                        else if(cigar[l] == "S" & next_clip_is_right == TRUE){
                            right_S_index <- l
                            for(m in (l -1):1){
                                if(cigar[m] != "M" & cigar[m] != "I" & cigar[m] != "D" & cigar[m] != "N" & cigar[m] != "S" & cigar[m] != "H"){
                                    reverse_right_clip <- paste(reverse_right_clip, cigar[m], sep = "") 
                                }
                                else{
                                    break
                                }
                            }
                            right_clip_length <- as.numeric(stri_reverse(reverse_right_clip)) 
                            if(i == 2){
                                new_end_position <- as.numeric(combined_position_data[[i]][j, 2]) + right_clip_length 
                            }
                            else if(i == 3){
                                new_start_position <- as.numeric(combined_position_data[[i]][j, 2]) + right_clip_length 
                            }
                            for(m in 1:nrow(combined_position_comparison_data[[i]])){
                                if(i == 2){
                                    if((new_end_position >= as.numeric(combined_position_comparison_data[[i]]$genomic_position[m]) - clip_window) & (new_end_position <= as.numeric(combined_position_comparison_data[[i]]$genomic_position[m]) + clip_window)){
                                        combined_position_comparison_data[[i]]$clipped_read_count[m] <- as.numeric(combined_position_comparison_data[[i]]$clipped_read_count[m]) + 1
                                        combined_position_comparison_data[[i]]$clipped_read_qnames[m] <- paste(combined_position_comparison_data[[i]]$clipped_read_qnames[m], combined_position_data[[i]]$qname[j], sep = ",")
                                        assigned_clipped_qnames <- c(assigned_clipped_qnames, unassigned_qnames[[i]][k])
                                        position_data$reverse_bam_read_position[position_data$qname == combined_position_data[[i]]$qname[j]] <- as.character(new_end_position)
                                        position_data$read_cigar[position_data$qname == combined_position_data[[i]]$qname[j]] <- paste(substring(combined_position_data[[i]]$read_cigar[j], 1, (right_S_index - 1)), "M", sep = "") 
                                    }
                                }
                                if(i == 3){
                                    if((new_start_position >= as.numeric(combined_position_comparison_data[[i]]$genomic_position[m]) - clip_window) & (new_start_position <= as.numeric(combined_position_comparison_data[[i]]$genomic_position[m]) + clip_window)){
                                        combined_position_comparison_data[[i]]$clipped_read_count[m] <- as.numeric(combined_position_comparison_data[[i]]$clipped_read_count[m]) + 1
                                        combined_position_comparison_data[[i]]$clipped_read_qnames[m] <- paste(combined_position_comparison_data[[i]]$clipped_read_qnames[m], combined_position_data[[i]]$qname[j], sep = ",")
                                        assigned_clipped_qnames <- c(assigned_clipped_qnames, unassigned_qnames[[i]][k])
                                        position_data$reverse_bam_read_position[position_data$qname == combined_position_data[[i]]$qname[j]] <- as.character(new_start_position)
                                        position_data$read_cigar[position_data$qname == combined_position_data[[i]]$qname[j]] <- paste(substring(combined_position_data[[i]]$read_cigar[j], 1, (right_S_index - 1)), "M", sep = "") 
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if(length(assigned_clipped_qnames) > 0){
        unassigned_qnames[[i]] <- unassigned_qnames[[i]][!(unassigned_qnames[[i]] %in% assigned_clipped_qnames)]
    }
    combined_position_comparison_data[[i]]$clipped_read_qnames[combined_position_comparison_data[[i]]$clipped_read_qnames == ""] <- "NA"
    combined_position_comparison_data[[i]] <- combined_position_comparison_data[[i]][order(as.numeric(combined_position_comparison_data[[i]]$genomic_position)), ]
}


### Separate reads into a dual assigned significant (reproducible 5' and 3' end), 5' end only, 3' end only, and an unassigned read sets ###
combined_assigned_qnames <- list(plus_start_assigned_qnames = c(), plus_end_assigned_qnames = c(), minus_start_assigned_qnames = c(), minus_end_assigned_qnames = c())
plus_strand_read_sets <- vector(mode = "list", length = 4)
minus_strand_read_sets <- vector(mode = "list", length = 4)
for(i in 1:length(combined_position_comparison_data)){
    for(j in 1:nrow(combined_position_comparison_data[[i]])){
        significant_read_count <- unlist(strsplit(combined_position_comparison_data[[i]]$significant_read_qnames[j], ","))
        noisy_read_count <- unlist(strsplit(combined_position_comparison_data[[i]]$noisy_read_qnames[j], ","))
        clipped_read_count <- unlist(strsplit(combined_position_comparison_data[[i]]$clipped_read_qnames[j], ","))
        combined_assigned_qnames[[i]] <- c(combined_assigned_qnames[[i]], significant_read_count, noisy_read_count, clipped_read_count)
    }
    combined_assigned_qnames[[i]] <- combined_assigned_qnames[[i]][!combined_assigned_qnames[[i]] == "NA" & !combined_assigned_qnames[[i]] == ""]
    if(i == 2){
        plus_qname_identification_data <- rbindlist(combined_position_data, idcol = "Identification_count", fill = TRUE)[qname %in% combined_assigned_qnames[[1]] & qname %in% combined_assigned_qnames[[2]]]
        plus_dual_assigned_qnames <- unique(plus_qname_identification_data$qname)
        plus_strand_read_sets[[1]] <- subset(position_data, qname %in% plus_dual_assigned_qnames)
        plus_strand_read_sets[[2]] <- subset(position_data, qname %in% combined_assigned_qnames[[1]] & !(qname %in% plus_dual_assigned_qnames))
        plus_strand_read_sets[[3]] <- subset(position_data, qname %in% combined_assigned_qnames[[2]] & !(qname %in% plus_dual_assigned_qnames))
        plus_strand_read_sets[[4]] <- subset(position_data, !(qname %in% combined_assigned_qnames[[1]]) & !(qname %in% combined_assigned_qnames[[2]]) & dna_strand == "+")
    }
    else if(i == 4){
        minus_qname_identification_data <- rbindlist(combined_position_data, idcol = "Identification_count", fill = TRUE)[qname %in% combined_assigned_qnames[[3]] & qname %in% combined_assigned_qnames[[4]]]
        minus_dual_assigned_qnames <- unique(minus_qname_identification_data$qname)
        minus_strand_read_sets[[1]] <- subset(position_data, qname %in% minus_dual_assigned_qnames)
        minus_strand_read_sets[[2]] <- subset(position_data, qname %in% combined_assigned_qnames[[3]] & !(qname %in% minus_dual_assigned_qnames))
        minus_strand_read_sets[[3]] <- subset(position_data, qname %in% combined_assigned_qnames[[4]] & !(qname %in% minus_dual_assigned_qnames))
        minus_strand_read_sets[[4]] <- subset(position_data, !(qname %in% combined_assigned_qnames[[3]]) & !(qname %in% combined_assigned_qnames[[4]]) & dna_strand == "-")
    }
}


### Isolate the read sets, outputting each read set as a SAM file ###
combined_dual_assigned_reads <- rbind(plus_strand_read_sets[[1]], minus_strand_read_sets[[1]])
combined_start_only_assigned_reads <- rbind(plus_strand_read_sets[[2]], minus_strand_read_sets[[2]])
combined_end_only_assigned_reads <- rbind(plus_strand_read_sets[[3]], minus_strand_read_sets[[3]])
combined_unassigned_reads <- rbind(plus_strand_read_sets[[4]], minus_strand_read_sets[[4]])
combined_reads <- list(combined_dual_assigned_reads, combined_start_only_assigned_reads, combined_end_only_assigned_reads, combined_unassigned_reads)
sam_set_names <- c("full_length_transcripts", "5p_intact_transcripts", "3p_intact_transcripts", "unassigned")
for(i in 1:length(combined_reads)){
    read_set_pad <- subset(sam_pad, sam_pad$qname %in% combined_reads[[i]]$qname)
    read_set_pad <- read_set_pad[order(as.numeric(read_set_pad$pos)), ]
    sam_set_pad <- rbind(sam_header, read_set_pad)
    sam_set <- list()
    for(j in 1:nrow(sam_set_pad)){
        row <- paste(sam_set_pad[j, ][!is.na(sam_set_pad[j, ])], collapse = "\t")
        sam_set[[j]] <- row
    }
    sam_set <- unlist(sam_set)
    write.table(sam_set, file = paste0(sam_set_names[i], ".sam"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE) 
}
