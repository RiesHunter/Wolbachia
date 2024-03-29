## Project overview
 - Mosquitoes fed on blood or ZIKV-infected mice
 - Some mosquitoes were infected with Wmel (*Wolbachia sp.*); others were given Tetracycline ("Tet")
 - We see reduced ZIKV dissemination (movement from midgut to legs to saliva) in the Wmel group
 - We are interested in how the ZIKV population evolves during Wolbachia co-infection

## Goals
 - Calculate barcode frequencies
 - Compare Tet to Wmel over time
   - Nucleotide diversity
   - Barcode frequency

### Process Illumina WGS fastq.gz files
Adapted file `ZIKV_OverlappingAmplicon_WGS_PE150_v3.1_LOOP2.command` as `/code/1-process.sh`
 - removed multi-lane processing
Adapted file `ZIKV_WGS_pipeline_Diversity_dinuc_v1_1.R` as `~/code/RStudio/1-process.R`
 - removed multi-lane processing
 - changed data output to various CSVs

### Extract and count barcodes
Adapted file `ZIKV_WGS_Barcode_analyses.command` as `~/code/2-count_barcodes.sh`
 - Mainly changed file names and paths
 - Also made it a for-loop
 - moved files into a similarly-named folder
 - removed multi-lane processing

### Remove barcode region from VCFs
Created `rm_BC_region_iSNVs.sh`
 - Removes iSNVs if they fall on 4007
 - Removes iSNVs if they fall between 4007 and 4030
 - Removes iSNVs if they fall on 4030
Ran this script on fn_ann VCF files

### R scripts
Adapted the following files into `2-barcode_freqs.R`:
 - `Mean_bc_frequencies.R`
 - `Euclidean_distance_meanbcfreqs.R`
 - `Barcodes_overtime_COMPOSITE.R`