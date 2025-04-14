# EBV transcriptomics tools

This repository serves strictly as a tutorial for the bioinformatics methodology behind [Mamane-Logsdon *et al*., 2025](), and as a place to store the scripts developed during the study. For data used or generated during the study, please follow links from [the front page of ebv.org.uk](http://ebv/org/uk).  

The methodology was developed as part of a pipeline for the characterisation of the Epstein-Barr virus (EBV) latency transcriptome based from nanopore direct RNA sequencing. This tutorial begins at the basecalling raw nanopore singal-level data in FAST5 format, and, therefore, assumes that the user has a FAST5 file. Also, for conciseness, this tutorial only focuses on the full-length transcript output of the methodology.

## Software dependencies

- [Guppy](https://nanoporetech.com/software/other/guppy) (v6.1.5)
- [minimap2](https://github.com/lh3/minimap2) (v2.24)
- [samtools](https://github.com/samtools/samtools) (v1.17)
- [R](https://www.r-project.org/) (v4.4.2), including the following libraries:
  - [data.table](https://cran.r-project.org/web/packages/data.table/index.html) (v1.16.4)
  - [plyr](https://cran.r-project.org/web/packages/plyr/index.html) (v1.8.9)
  - [stringi](https://cran.r-project.org/web/packages/stringi/index.html) (v1.8.4)
  - [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html) (v2.0.0)
- [Mamane-Logsdon *et al*., 2025 scripts]() (initial release)
- [nanopolish](https://github.com/jts/nanopolish) (v0.14.0)
- [perl](https://www.perl.org/) (v5.34.1)
- [Donovan-Banfield *et al*., 2020 scripts](https://zenodo.org/records/3610249) (v3)

## Tutorial

### Part 1: Basecalling raw nanopore signal-level data

**1.1. Basecall your FAST5 file using Guppy**:

```bash
guppy_basecaller --input_path fast5 --save_path fastq --flowcell flowcell --kit kit
```

To replace:
- `guppy_basecaller`: Path to your guppy_basecaller executable file.
- `fast5`: Path to the directory that stores your FAST5 file.
- `flowcell`: Flowcell used in your sequencing experiment, e.g., FLO-MIN106.
- `kit`: Kit used in your sequencing experiment, e.g., SQK-RBK004.
- `fastq`: Path to the directory that will store your basecalled output, i.e., the FASTQ file.

<br>

### Part 2: Generating an EBV genome alignment

**2.1. Perform annotation-guided splice-aware alignment of your FASTQ file to your EBV genome of choice using minimap2**:

```bash
minimap2 -ax splice -k14 --junc-bed genome_annotation.bed genome.fasta reads.fastq > alignment.sam
```

To replace:
- `minimap2`: Path to your minimap2 executable file.
- `reads.fastq`: Path to your FASTQ file.
- `genome.fasta`: Path to your FASTA file of your EBV genome of choice.
- `genome_annotation.bed`: Path to your BED file of EBV genome annotation.
- `alignment.sam`: Path to your output SAM file.

**2.2. Remove unaligned reads and secondary alignments from your SAM file using samtools**:

```bash
samtools view --excl-flags 260 --with-header alignment.sam > filtered_alignment.sam
```

To replace:
- `samtools`: Path to your samtools executable file.
- `alignment.sam`: Path to your SAM file.
- `filtered_alignment.sam`: Path to your output filtered SAM file.

<br>

### Part 3: Characterising and quantifying the full-length EBV transcripts

**3.1. Add back the W0 exon to reads in your filtered SAM file from which it has been removed by 5' soft-clipping using W0_correcting.r**:

```bash
Rscript W0_correcting.r filtered_alignment.sam exons.bed
```

To replace:
- `W0_correcting.r`: Path to your W0_correcting.r script from Mamane-Logsdon *et al*., 2025.
- `filtered_alignment.sam`: Path to your filtered SAM file.
- `exons.bed`: Path to your BED file containing the start and end coordinates of all your chosen EBV genome's W0, W1, and W1' exons. See [here]() for the expected format.

Output (not shown):
- `candidates_W0_id.sam`: SAM file containing candidate reads for W0 correction (by virtue of starting within 5 bp of the start of a W1 or W1' exon) that contain at least one of the W0-associated motifs ACAAT, ACAAAT, AGGAGT, or CAGGAGT within or adjacent to their 5′ soft clip, but not the W2-associated motif CAGGG.
-  `candidates_W2_id.sam`: SAM file containing candidate reads for W0 correction (by virtue of starting within 5 bp of the start of a W1 or W1' exon) that contain the aforementioned W2-associated motif within or adjacent to their 5′ soft clip, but not any of the aforementioned W0-associated motifs.
- `candidates_mixed_id.sam`: SAM file containing candidate reads for W0 correction (by virtue of starting within 5 bp of the start of a W1 or W1' exon) that contain both the aforementioned W0- and W2-associated motifs within or adjacent to their 5′ soft clip.
- `candidates_unassigned_id.sam`: SAM file containing candidate reads for W0 correction (by virtue of starting within 5 bp of the start of a W1 or W1' exon) that do not contain any of the aforementioned W0- and W2-associated motifs within or adjacent to their 5′ soft clip.
- `updated_filtered_alignment.sam`: SAM file containing all the reads from the input SAM file but with the W0 exon added to the start of reads that were candidates for W0 correction and also shown to have only W0 exon identification in their 5' region. This is produced by unclipping the 5′ soft-clipped bases of all candidates reads. Then, for candidate reads shown to have only W0 exon identification within or adjacent to their 5′ soft clip, it is calculated where the W0 exon would end within the read sequence based on where the W0-associated motif ends. Finally, it inserts a gap in the read between the W0 exon end and the rest of the read to represent the missing splice boundary so that the W0 exon is effectively added back into the read alignment.

**3.2. Classify reads in your W0 exon corrected SAM file into full-length transcripts, broken transcripts, and unassigned reads using transcript_classifying.r**:

```bash
Rscript transcript_classifying.r updated_filtered_alignment.sam significance_window mininum_significant_read_count noise_window clip_window
```

To replace:
- `transcript_classifying.r`: Path to your transcript_classifying.r script from Mamane-Logsdon *et al*., 2025.
- `updated_filtered_alignment.sam`: Path to your W0 corrected SAM file.
- `significance_window`: Your window size in bp for grouping read ends to define a significant 5′ or 3′ end, e.g., 2. This is the first window applied.
- `mininum_significant_read_count`: Your minimum number of reads within your significance window required to define a significant 5' or 3' end, e.g., 5.
- `noise_window`: Your window size in bp for capturing unassigned read ends near significant 5' or 3' ends to assign them to said significant ends, e.g., 20. This is done without unclipping clipped bases of unassigned reads. This is the second window applied.
- `clip_window`: Your window size in bp for capturing unassigned read ends near significant 5' or 3' ends to assign them to said significant ends, e.g., 20. This is done with unclipping clipped bases of unassigned reads. This is the third and final window applied.

Output (not shown):
- `full_length_transcripts.sam`: SAM file containing the reads assigned to have both significant 5' and 3' ends. These reads are considered full-length transcripts.
- `5p_intact_transcripts.sam`: SAM file containing the reads assigned to have only a significant 5' end. These reads are considered broken transcripts. 
- `3p_intact_transcripts.sam`: SAM file containing the reads assigned to have only a significant 3' end. These reads are considered broken transcripts. 
- `unassigned.sam`: SAM file containing the reads that weren't assigned any significant ends. These reads have unclear identities.

**3.3. Count the number of IR1 repeat units in your full-length EBV transcripts using IR1_repeat_counting.r**:

```bash
Rscript IR1_repeat_counting.r full_length_transcripts.sam exons.bed start_end_window splice_window
```

To replace:
- `IR1_repeat_counting.r`: Path to your IR1_repeat_counting.r script from Mamane-Logsdon *et al*., 2025.
- `full_length_transcripts.sam`: Path to your full-length transcripts SAM file.
- `exons.bed`: Path to your BED file containing the start and end coordinates of all your chosen EBV genome's exons of interest. At minimum, exons C1, C2, W0, W1, W1', W2, and Y1 must be given. See [here]() for the expected format.
- `start_end_window`: Your window size in bp for how many bases outside of an exon’s start or end coordinates can still be considered part of that exon when annotating transcript structure.
- `splice_window`: Your window size in bp for how many bases is close enough to an exon junction to be treated as a proper splice site when annotating transcript structure.

Output (not shown):
- `exon_contents.txt`: TXT file listing each full-length transcript alongside the lengths of each aligned exon that the transcript contains. It is produced by looping through each transcript, identifying its aligned exons based on exon start-end coordinates and the user-specified windows, and then stores the aligned exon lengths in bp.
- `IR1_repeat_counts.txt`: TXT file that details each full-length transcript's promoter type (Cp or Wp), whether an initial W1 or W1’ exon is used, and two methods of IR1 repeat unit counting, one based on alignment that counts the number of aligned IR1 repeats by mimimap2 and another based on distance that calculates the number of exon nucleotides between exons C2 and Y2 and divides it by the average IR1 repeat unit length. The latter count method is to account for misalignment by minimap2 and is considered the ground truth for IR1 repeat unit count.

<br>

### Part 4: Determining the polyA tail status of full-length EBV transcripts

**4.1. Convert your full-length transcripts SAM file to a BAM file using samtools**:

```bash
samtools view --bam full_length_transcripts.sam > full_length_transcripts.bam
```

To replace:
- `samtools`: Path to your samtools executable file.
- `full_length_transcripts.sam`: Path to your full-length transcripts SAM file.
- `full_length_transcripts.bam`: Path to your output full-length transcripts BAM file.

**4.2. Sort your full-length transcripts BAM file using samtools**:

```bash
samtools sort full_length_transcripts.bam > full_length_transcripts_sort.bam
```

To replace:
- `samtools`: Path to your samtools executable file.
- `full_length_transcripts.bam`: Path to your full-length transcripts BAM file.
- `full_length_transcripts_sort.bam`: Path to your output sorted full-length transcripts BAM file.

**4.3. Index your sorted full-length transcripts BAM file using samtools**:

```bash
samtools index full_length_transcripts_sort.bam
```

To replace:
- `samtools`: Path to your samtools executable file.
- `full_length_transcripts_sort.bam`: Path to your sorted full-length transcripts BAM file.

Output (not shown):
- `full_length_transcripts_sort.bam.bai`: Index of your sorted full-length transcripts BAM file

**4.4. Index your FASTQ file from before using nanopolish**:

```bash
nanopolish index --directory fast5 reads.fastq
```

To replace:
- `nanopolish`: Path to your nanopolish executable file.
- `fast5`: Path to the directory that stores your FAST5 file.
- `reads.fastq`: Path to your FASTQ file.

Output (not shown):
- `reads.fastq.index`: gzipped FASTA file of your FASTQ file.
- `reads.fastq.index.fai`, `reads.fastq.index.gzi`, `reads.fastq.index.readdb`: Index files of your FASTQ file.

**4.5. Determine the polyA tail status of your full-length transcripts using nanopolish**:

```bash
nanopolish polya --threads thread_number --reads reads.fastq --bam full_length_transcripts_sort.bam --genome genome.fasta > polya_results.tsv
```

To replace:
- `nanopolish`: Path to your nanopolish executable file.
- `thread_number`: Number of threads you want to use for computation, e.g., 8.
- `reads.fastq`: Path to your FASTQ file.
- `full_length_transcripts_sort.bam`: Path to your sorted full-length transcripts BAM file.
- `genome.fasta`: Path to your FASTA file of your EBV genome of choice.
- `polya_results.tsv`: Path to your output TSV file of polyA tail status of the full-length transcripts.

<br>

### Part 5: Grouping full-length EBV transcripts and correcting their 5' ends

**5.1. Convert your sorted full-length transcripts BAM file to a sorted SAM file using samtools**:

```bash
samtools view --with-header full_length_transcripts_sort.bam > full_length_transcripts_sort.sam
```

To replace:
- `samtools`: Path to your samtools executable file.
- `full_length_transcripts_sort.bam`: Path to your sorted full-length transcripts BAM file.
- `full_length_transcripts_sort.sam`: Path to your output sorted full-length transcripts SAM file.

**5.2. Group your full-length transcripts using classify_transcripts_and_polya_segmented_V2.pl**:

```bash
perl classify_transcripts_and_polya_segmented_V2.pl prefix polya/TSS_window splice_window polya_min polya_results.tsv full_length_transcripts_sort.sam genome.fasta min_copy max_entries
```

To replace:
- `classify_transcripts_and_polya_segmented_V2.pl`: Path to your classify_transcripts_and_polya_segmented_V2.pl script from Donovan-Banfield *et al*., 2020. 
- `prefix`: Your prefix for output files, e.g., ebv.
- `polya/TSS_window`: Your window size for grouping transcription start or stop sites together, e.g., 20.
- `splice_window`: Your window size for grouping splice sites together, e.g., 2.
- `polya_min`: Your minimum length for a polyA tail before a transcript will be analysed, e.g., 8.
- `polya_results.tsv`: Path to your TSV file of polyA tail status of the full-length transcripts.
- `full_length_transcripts_sort.sam`: Path to your sorted full-length transcripts SAM file.
- `genome.fasta`: Path to your FASTA file of your EBV genome of choice.
- `min_copy`: Your minimum copy number in a transcript group to consider, e.g., 1.
- `max_entries`: Your maximum number of entries in the list of most abundant transcript groups, e.g., 80.

Output (not shown):
- `ebv.canonical_transcripts_by_abundance.fasta`: FASTA file listing pseudo transcripts generated by using the start, end and splice patterns of the transcript groups found, one transcript per transcript group along with information on how many individual transcripts belong to that transcript group, the average polyA tail length, and the mapping co-ordinates.
- `ebv.GFF_all_found.gff3`: GFF3 file describing each transcript group.
- `ebv.GFF_most_popular_list.gff3`: GFF3 file only describing the most abundant transcript groups.
- `ebv.raw_soft_clip_locations.txt`: TXT file of the locations on the genome where there is soft clipping on either the plus or minus strand.
- `ebv.start_sad_stop_pattern_count.txt`: TXT file describing each transcript group's start location, end location, strand, splice acceptor/donor locations (referred to as the sad location), the average and standard deviation for the polyA tail length, and how many transcripts belong to that group.
- `ebv.raw_processed_start_sad_polya.txt`: TXT file describing for each nucleotide position how often that is the location of a transcript start, end, or a splice acceptor/donor.

**5.3. Correct the 5’ ends of your grouped transcripts**:

```bash
perl name_transcripts_and_track_ssc_V2.pl prefix ssc_nt features_table.txt genome.fasta ebv.start_sad_stop_pattern_count.txt
```

To replace:

- `name_transcripts_and_track_ssc_V2.pl`: Path to your name_transcripts_and_track_ssc_V2.pl script from Donovan-Banfield *et al*., 2020.
- `prefix`: Your prefix for output files, e.g., ebv.
- `ssc_nt`: Your number of upstream nucleotides for each pseudo transcript to include in order to account for the loss of 5’ nucleotides in nanopore sequencing, e.g., 10.
- `features_table.txt`: Path to your TXT file containing a list of features on the genome being analysed. See [here]() for the expected format.
- `genome.fasta`: Path to your FASTA file of your EBV genome of choice.
- `ebv.start_sad_stop_pattern_count.txt`: Path to your TXT file describing each transcript group's start location, end location, strand, splice acceptor/donor locations (referred to as the sad location), the average and standard deviation for the polyA tail length, and how many transcripts belong to that group.

Output (not shown):
- `ebv.TSS_start_sad_polyA_count_list.txt`: TXT file listing all the transcript groups analysed.
- `ebv.proteins_known.fasta`: FASTA file of the known proteins found after analysing the canonical features provided.
- `ebv.proteins_not_known.fasta`: FASTA file of the unknown proteins found after analysing the canonical features provided.
- `ebv.most_abundant_trans_per_orf.gff`: GFF file describing the most abundant transcript group that will code for each one of the canonical ORFs.
- `ebv.count_of_translated_features.txt`: TXT file detailing, for each ORF identified, how many transcripts in total would code for that ORF and the average polyA tail length for the dominant transcript group that codes for the indicated ORF.
- `ebv.combined_counts_and_feature_names.txt`: TXT file containing a summary for each transcript group.
- `ebv.splice_acceptor_donor_usage.txt`: TXT file listing splice acceptor donor pair usage (in message sense only) and how often each one is used.
- `ebv`: Folder named after the prefix storing a series of GFF3 files, one for each ORF found.
