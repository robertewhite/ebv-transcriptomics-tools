#!/usr/bin/Rscript --vanilla

suppressMessages(suppressWarnings(library(plyr)))
suppressMessages(suppressWarnings(library(stringi)))


### Initiating user input ###
args <- commandArgs(trailingOnly = TRUE)
sam_path <- args[1]
exons_path <- args[2]
start_end_window <- as.numeric(args[3])
splice_window <- as.numeric(args[4])


### Open and store the SAM file ###
sam_conn <- file(sam_path, open = "r")
sam_file <- list()
for(j in readLines(sam_conn, warn = FALSE)){
  split_j <- unlist(strsplit(j, split = "\t"))
  sam_file <- c(sam_file, list(split_j))
}
close(sam_conn)


### Calculate the maximum column number of the SAM file ###
col_max <- 0
for(j in 1:length(sam_file)){
  if(length(sam_file[[j]]) > col_max){
    col_max <- length(sam_file[[j]])
  }
}


### Pad the SAM file rows with NAs based on the maximum column number ###
sam_pad <- data.frame(matrix(nrow = 0, ncol = col_max))
for(j in 1:length(sam_file)){
  if(length(sam_file[[j]]) < col_max){
    sam_file[[j]] <- c(sam_file[[j]], rep(NA, length.out = col_max - length(sam_file[[j]])))
  }
  sam_pad <- rbind(sam_pad, sam_file[[j]])
}
colnames(sam_pad) <- c("qname", "flag", "rname", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq", "qual", paste0("V", 12:(col_max)))
sam_header <- subset(sam_pad, grepl(sam_pad$qname, pattern = "^@"))
sam_reads <- subset(sam_pad, !grepl(sam_pad$qname, pattern = "^@"))


### Calculate the size of each BED file exon in each read ###
exons <- read.table(exons_path, sep = "\t", header = TRUE)
sam_read_exon_lengths <- data.frame(matrix(ncol = nrow(exons) + 2, nrow = 0))
for(j in 1:nrow(sam_reads)){
    sam_read_exons <- vector(mode = "character", length = nrow(exons))
    sam_read_exons[1:length(sam_read_exons)] <- "0"
    cigar <- unlist(strsplit(sam_reads$cigar[j], split = ""))
    end_of_cigar <- FALSE
    exon_start <- as.numeric(sam_reads$pos[j])
    exon_end <- exon_start - 1
    exon_length <- 0
    first_exon_found <- FALSE
    first_exon_is_W0 <- FALSE
    letter_other_than_S_found <- FALSE
    letter_other_than_H_found <- FALSE
    for(k in 1:length(cigar)){
        reverse_event <- ""
        event_length <- 0
        if(cigar[k] == "M" | cigar[k] == "D" | cigar[k] == "I" | cigar[k] == "N" | cigar[k] == "S" | cigar[k] == "H"){
            if(k == length(cigar)){
                end_of_cigar <- TRUE
            }
            for(l in (k - 1):1){
                if(cigar[l] != "M" & cigar[l] != "I" & cigar[l] != "D" & cigar[l] != "N" & cigar[l] != "S" & cigar[l] != "H"){ 
                    reverse_event <- paste(reverse_event, cigar[l], sep = "")
                }
                else{
                    break
                }
            }
            event_length <- as.numeric(stri_reverse(reverse_event))
            if(cigar[k] == "M" | cigar[k] == "D"){
                letter_other_than_S_found <- TRUE
                letter_other_than_H_found <- TRUE
                exon_end <- exon_end + event_length
            }
            if(cigar[k] == "M" | cigar[k] == "I" | cigar[k] == "D"){
                letter_other_than_S_found <- TRUE
                letter_other_than_H_found <- TRUE
                if(cigar[k] == "M"){
                    exon_length <- exon_length + event_length
                }
                else if((cigar[k] == "I" & event_length >= 2)){
                    exon_length <- exon_length + event_length
                }
                else if(cigar[k] == "D" & event_length <= 2){
                    exon_length <- exon_length + event_length
                }

            }
            if(cigar[k] == "N" & first_exon_found == FALSE){
                for(l in 1:nrow(exons)){
                    if((abs(diff(c(exon_start, as.numeric(exons$start[l])))) <= start_end_window) & (abs(diff(c(exon_end, as.numeric(exons$end[l])))) <= splice_window)){
                        if(exons$exon[l] == "W1"){
                            if(abs(diff(c(exon_start, as.numeric(exons$start[l])))) < abs(diff(c(exon_start, as.numeric(exons$start[l + 1]))))){
                                sam_read_exons[l] <- exon_length
                                first_exon_found <- TRUE
                                break
                            }
                            else if(abs(diff(c(exon_start, as.numeric(exons$start[l])))) > abs(diff(c(exon_start, as.numeric(exons$start[l + 1]))))){
                                sam_read_exons[l + 1] <- exon_length
                                first_exon_found <- TRUE
                                break
                            }
                        }
                        sam_read_exons[l] <- exon_length
                        if(exons$exon[l] == "W0"){ 
                            first_exon_is_W0 <- TRUE
                        }
                        break
                    }
                }
                exon_start <- exon_end + event_length + 1
                exon_end <- exon_start - 1
                exon_length <- 0
                if(first_exon_is_W0 == TRUE){ 
                    next
                }
                first_exon_found <- TRUE
                next
            }
            else if((cigar[k] == "N" & first_exon_found == TRUE)){
                for(l in 1:nrow(exons)){
                    if(((abs(diff(c(exon_start, as.numeric(exons$start[l])))) <= splice_window) & (abs(diff(c(exon_end, as.numeric(exons$end[l])))) <= splice_window))){
                        if(exons$exon[l] == "W1"){
                            if(abs(diff(c(exon_start, as.numeric(exons$start[l])))) < abs(diff(c(exon_start, as.numeric(exons$start[l + 1]))))){
                                sam_read_exons[l] <- exon_length
                                break
                            }
                            else if(abs(diff(c(exon_start, as.numeric(exons$start[l])))) > abs(diff(c(exon_start, as.numeric(exons$start[l + 1]))))){
                                sam_read_exons[l + 1] <- exon_length
                                break
                            }
                        }
                        sam_read_exons[l] <- exon_length
                        break
                    }
                }
                exon_start <- exon_end + event_length + 1
                exon_end <- exon_start - 1
                exon_length <- 0
            }
            else if(end_of_cigar == TRUE){
                for(l in 1:nrow(exons)){
                    if(first_exon_found == TRUE){
                        if(abs(diff(c(exon_start, as.numeric(exons$start[l])))) <= splice_window & abs(diff(c(exon_end, as.numeric(exons$end[l])))) <= start_end_window){
                            sam_read_exons[l] <- exon_length
                            break
                        }
                    }
                    else if(first_exon_found == FALSE){
                        if(abs(diff(c(exon_start, as.numeric(exons$start[l])))) <= start_end_window & abs(diff(c(exon_end, as.numeric(exons$end[l])))) <= start_end_window){
                            sam_read_exons[l] <- exon_length
                            break
                        }
                    }
                }
            }
        }
    }
    sam_read_exon_lengths <- rbind(sam_read_exon_lengths, c(sam_reads$qname[j], sam_read_exons))
}
colnames(sam_read_exon_lengths) <- c("qname", exons$exon)


### Calculate the number of IR1 repeats (according to the alignment tool used) in each read ###
aligner_IR1_repeat_counts <- data.frame(matrix(ncol = 4, nrow = 0))
sam_read_W_exon_lengths <- sam_read_exon_lengths[, grepl("^W", names(sam_read_exon_lengths))]
sam_read_W_exon_lengths <- cbind(sam_read_W_exon_lengths, sam_read_exon_lengths$qname)
colnames(sam_read_W_exon_lengths)[ncol(sam_read_W_exon_lengths)] <- "qname"
for(j in 1:nrow(sam_read_W_exon_lengths)){
    promoter <- ""
    repeat_count <- 0
    W0_count <- 0
    W0_used <- FALSE
    W1_or_W1_prime_used <- "NA"
    for(k in seq(1, ncol(sam_read_W_exon_lengths) - 1, by = 5)){
        if(sam_read_W_exon_lengths[j, k] == "0" & W0_used == FALSE){
            W0_count <- W0_count + 1
        }
        else if(sam_read_W_exon_lengths[j, k] != "0" & W0_used == FALSE){
            W0_count <- W0_count + 1
            W0_used <- TRUE
        }
        if(sam_read_W_exon_lengths[j, k + 1] != "0" | sam_read_W_exon_lengths[j, k + 2] != "0"){
            repeat_count <- repeat_count + 0.5
            if(W1_or_W1_prime_used  == "NA" & sam_read_W_exon_lengths[j, k + 1] != "0"){
                W1_or_W1_prime_used <- "W1"
            }
            else if(W1_or_W1_prime_used  == "NA" & sam_read_W_exon_lengths[j, k + 2] != "0"){
                W1_or_W1_prime_used <- "W1_prime"
            }
        }
        if(sam_read_W_exon_lengths[j, k + 3] != "0" | sam_read_W_exon_lengths[j, k + 4] != "0"){
            repeat_count <- repeat_count + 0.5
        }
    }
    if(repeat_count == 0 & as.numeric(sam_reads$pos[sam_reads$qname == sam_read_W_exon_lengths$qname[j]]) >= as.numeric(exons$end[exons$exon == "Y1"])){
        repeat_count <- NA 
    }
    if(W0_used == FALSE & as.numeric(sam_reads$pos[sam_reads$qname == sam_read_W_exon_lengths$qname[j]]) <= as.numeric(exons$end[exons$exon == "C1"])){
        promoter <- "Cp"
    }
    else if(W0_used == FALSE & as.numeric(sam_reads$pos[sam_reads$qname == sam_read_W_exon_lengths$qname[j]]) > as.numeric(exons$end[exons$exon == "C1"])){
        promoter <- NA
    }
    else if(W0_used == TRUE & as.numeric(sam_reads$pos[sam_reads$qname == sam_read_W_exon_lengths$qname[j]]) > as.numeric(exons$end[exons$exon == "C2"])){
        promoter <- paste("Wp", W0_count, sep = "")
    }
    aligner_IR1_repeat_counts <- rbind(aligner_IR1_repeat_counts, c(sam_read_W_exon_lengths$qname[j], promoter, W1_or_W1_prime_used, repeat_count))
}
colnames(aligner_IR1_repeat_counts) <- c("qname", "promoter", "W1_vs_W1_prime", "alignment_count")


### Calculate the number of IR1 repeats (according to nucleotide distance) in each read ###
nucleotide_IR1_repeat_counts <- c()
average_W1W2_length <- mean(c(198, 185.5, 193, 180.5))
C2_end_position <- exons$end[exons$exon == "C2"]
Y1_start_position <- exons$start[exons$exon == "Y1"]
for(j in 1:nrow(sam_read_exon_lengths)){
    reference_position <- as.numeric(sam_reads$pos[sam_reads$qname == sam_read_exon_lengths$qname[j]]) - 1
    repeat_count <- 0
    nucleotide_count <- 0
    cigar <- unlist(strsplit(sam_reads$cigar[sam_reads$qname == sam_read_exon_lengths$qname[j]], split = ""))
    for(k in 1:length(cigar)){
        reverse_event <- ""
        event_length <- 0
        if(reference_position <= as.numeric(C2_end_position)){
            if(cigar[k] == "M" | cigar[k] == "D" | cigar[k] == "N"){
                for(l in (k - 1):1){
                    if(cigar[l] != "M" & cigar[l] != "I" & cigar[l] != "D" & cigar[l] != "N" & cigar[l] != "S" & cigar[l] != "H"){
                        reverse_event <- paste(reverse_event, cigar[l], sep = "")
                    }
                    else{
                        break
                    }
                }
                event_length <- as.numeric(stri_reverse(reverse_event))
                reference_position <- as.numeric(reference_position) + event_length
            }
        }
        else if(reference_position > as.numeric(C2_end_position) & reference_position < (as.numeric(Y1_start_position) - 1)){
            if(cigar[k] == "M" | cigar[k] == "D" | cigar[k] == "I" | cigar[k] == "N"){
                for(l in (k - 1):1){
                    if(cigar[l] != "M" & cigar[l] != "I" & cigar[l] != "D" & cigar[l] != "N" & cigar[l] != "S" & cigar[l] != "H"){
                        reverse_event <- paste(reverse_event, cigar[l], sep = "")
                    }
                    else{
                        break
                    }
                }
                event_length <- as.numeric(stri_reverse(reverse_event))
                if(cigar[k] == "M" | cigar[k] == "D" | cigar[k] == "N"){
                    reference_position <- as.numeric(reference_position) + event_length
                }
                if(cigar[k] == "M"){
                    if(reference_position >= (as.numeric(Y1_start_position) - 1)){
                        nucleotide_count <- nucleotide_count + (reference_position - (as.numeric(Y1_start_position) - 1))
                        break
                    }
                    nucleotide_count <- nucleotide_count + event_length
                }
                else if(cigar[k] == "I"){
                    if(event_length >= 2){
                        if(reference_position >= (as.numeric(Y1_start_position) - 1)){
                            nucleotide_count <- nucleotide_count + (reference_position - (as.numeric(Y1_start_position) - 1))
                            break
                        }
                        nucleotide_count <- nucleotide_count + event_length
                    }
                }
                else if(cigar[k] == "D"){
                    if(event_length <= 2){
                        if(reference_position >= (as.numeric(Y1_start_position) - 1)){
                            nucleotide_count <- nucleotide_count + (reference_position - (as.numeric(Y1_start_position) - 1))
                            break
                        }
                        nucleotide_count <- nucleotide_count + event_length
                    }
                }
                else if(cigar[k] == "N"){
                    if(reference_position >= (as.numeric(Y1_start_position) - 1)){
                        break
                    }
                }
            }
        }
    }
    if(!is.na(aligner_IR1_repeat_counts$promoter[j]) & grepl("^Wp", aligner_IR1_repeat_counts$promoter[j])){
        W0_exon_lengths <- sam_read_exon_lengths[j, grepl("^W0", names(sam_read_exon_lengths))]
        for(k in 1:length(W0_exon_lengths)){
            if(W0_exon_lengths[k] != "0"){
                nucleotide_count <- nucleotide_count - as.numeric(W0_exon_lengths[k])   
            }
        }
    }
    if(nucleotide_count > 0){
        repeat_count <- nucleotide_count / average_W1W2_length
        repeat_count <- round_any(repeat_count, 0.1)
    }
    else if(nucleotide_count == 0){
        repeat_count <- NA
    }
    nucleotide_IR1_repeat_counts <- c(nucleotide_IR1_repeat_counts, repeat_count)
}
combined_IR1_repeat_counts <- cbind(aligner_IR1_repeat_counts, nucleotide_IR1_repeat_counts)
colnames(combined_IR1_repeat_counts) <- c(colnames(aligner_IR1_repeat_counts), "distance_count")
combined_IR1_repeat_counts <- na.omit(combined_IR1_repeat_counts)
sam_read_exon_lengths <- sam_read_exon_lengths[sam_read_exon_lengths$qname %in% combined_IR1_repeat_counts$qname, ]


### Output the exons contents and the IR1 repeat counts for Cp/Wp reads ###
write.table(sam_read_exon_lengths, "exon_contents.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(combined_IR1_repeat_counts, file = "IR1_repeat_counts.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
