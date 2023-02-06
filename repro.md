## Project overview
 - Mosquitoes fed on blood or ZIKV-infected mice
 - Some mosquitoes were "infected" with Wmel (*Wolbachia sp.*); others were given Tetracycline ("Tet")
 - We see reduced ZIKV dissemination (movement from midgut to legs to saliva) in the Wmel group
 - We are interested in how the ZIKV population evolves during Wolbachia infection

## Goal
 - Calculate Pi, PiN/PiS, mutation frequency, Gini-Simpson, and Shannon Entropy for all samples
 - Compare metrics between Tet and Wmel over time
 - Investigate bottleneck sizes through barcode analyses

# 230201
## Objectives:
 - Calculate barcode frequencies
 - Compare Tet to Wmel over time

### Process Illumina WGS fastq.gz files
Adapted file `ZIKV_OverlappingAmplicon_WGS_PE150_v3.1_LOOP2.command` as `/code/1-process.sh`
 - removed multi-lane processing
Adapted file `ZIKV_WGS_pipeline_Diversity_dinuc_v1_1.R` as `/code/RStudio/1-process.R`
 - removed multi-lane processing
 - changed data output to various CSVs

### Extract and count barcodes
Adapted file `ZIKV_WGS_Barcode_analyses.command` as `/code/2-count_barcodes.sh`
 - Mainly changed file names and paths
 - Also made it a for-loop
 - moved files into a similarly-named folder
 - removed multi-lane processing

### R scripts
Adapted the following files into `2-barcode_freqs.R`:
 - `Mean_bc_frequencies.R`
 - `Euclidean_distance_meanbcfreqs.R`
 - `Barcodes_overtime_COMPOSITE.R`
 - 