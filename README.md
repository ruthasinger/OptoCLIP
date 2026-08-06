# Opto-CLIP Analysis Pipeline



## Overview

This repository includes all scripts used to generate the bioinformatic components main and supplementary figures for the following manuscript:

**Neuronal activation triggers rapid biphasic FMRP-mediated translational control coupling synaptic and nuclear responses**

This study implements a multi-modal analysis pipeline integrating FMRP CLIP and RiboTag RNA-seq to investigate activity-dependent RNA regulation in CA1 neurons following optogenetic stimulation.



## Repository Structure
```
├── Figure2_FigureS2_Ephys/
  └── OptoCLIP_Figure2_FigureS2_Rcode_final.R
  └── Data
      └── 5trains5Hz
      └── CurrentClamp
      └── DifferentHz
├── Figure3_FigureS3_FMRP_CLIP/
  └── OptoCLIP_Figure3_FigureS3_Rcode_final.R
  └── Data
      └── Data_FMRP_CLIP_Tx_STAR
├── Figure4_FigureS4_RiboTag/
  └── OptoCLIP_Figure4_FigureS4_Rcode_final.R
  └── Data
      └── Data_FMRP_CLIP_Tx_STAR
      └── Data_RiboTag
├── Figure5_FigureS5_Figure6_FigureS7/
  └── OptoCLIP_Figure5_FigureS5_Figure6_FigureS7_Rcode_final.R
├── R_functions/
├── RDS_files/
```



## **Figures**

### Figure 2 and Figure S2 — Electrophysiology

#### **Script:**

- OptoCLIP_Figure2_FigureS2_Rcode_final.R

#### **Input data:**

- **CurrentClamp**

  This folder contains 3 raw abf files corresponding to 3 different patched cells that underwent stepwise injection of -50 to 200 pA current in 25 pA increments.

- **DifferentHz**

  This folder contains 6 raw abf files. The files correspond to either ChR2:mCherry+ cells or 3 Control (ChR2 negative) cells that were patched during exposure to 3 different blue LED frequencies (1, 5, and 10 Hz)

- **5trains5Hz**

  This folder contains 5 raw abf files corresponding to 3 ChR2:mCherry+ cells or 2 Control (ChR2 negative) cells that were patched during exposure to 5 trains of 30" 5 Hz blue LED light followed by 30" light off

#### **Output (11 plots):**

- **Figures 2c-e:** Whole cell patch clamp traces of ChR2+ and ChR2- cells exposed to 10 seconds of 1Hz, 5Hz, or 10Hz blue LED light pulses

- **Figures 2c'-e':** Zoomed in 1 second traces of 2C-2E

- **Figures 2f, g:** Whole cell patch clamp traces of ChR2+ and ChR2- cells exposed to trains of 30" 5Hz blue light followed by 30" dark repeated for a total of 5 trains. 

- **Supplementary Figures 2a, b:** Current clamp traces of ChR2+ cells that underwent stepwise injection of -50 to 200 pA current in 25 pA increments. Rheobase current (+50 pA) is highlighted in green

- **Supplementary Figure 2c:** Quantification of number of induced action potentials by pA of injected current. Shading indicates +/- SEM from 3 different patched cells. 



### Figure 3 and Figure S3 — Opto FMRP-CLIP

#### **Script:**

- OptoCLIP_Figure3_FigureS3_Rcode_final.R

#### **Input data:**

- **Data_FMRP_CLIP_Tx_STAR**

  This folder contains 15 bed files (*.transcriptome.unique.bed). 3 Cre negative samples, 4 Control samples, 4 Opto_5minute samples, and 4 Opto_30minute samples. Code for generating these files from raw, demultiplexed fastq files (found on GEO: GSE286379) can be found in this data folder entitled "RAS_bashscript_FMRPCLIP_fastq2transcriptome"

#### **Output (5 plots, 2 source data csv files):** 

- **Figure 3d:** Quantification of uniquely mapped FMRP CLIP tags by condition. 

- **Figure 3e:** Bar plots showing genomic annotations of uniquely mapped FMRP CLIP tags separated by transcripts and regions found in increasing biological replicates.

- **Figure 3f:** Meta transcript coverage of 1000 transcripts with most CLIP tags separated by condition.

- **Figure 3g:** Boxplots depicting the distribution of transcript-level CLIP tag proportions in 5′UTR, CDS, and 3′UTR regions separated by condition.

- **Supplementary Figure 3a:** Heatmap of pairwise pearson correlation analysis of each sample compared to every other sample.

- **Source Data for Figure 3d**

- **Source Data for Figure 3g**

- **Source Data for Supplementary Figure 3a**

### Figure 4 and Figure S4 — Opto RiboTag

#### **Script:**

- OptoCLIP_Figure4_FigureS4_Rcode_final.R

#### **Input data:**

- **Data_FMRP_CLIP_Tx_STAR**

  This folder contains 16 bed files (*.transcriptome.unique.bed). 4 Control_5minute samples, 4 Control_30minute samples, 4 Opto_5minute samples, and 4 Opto_30minute samples. 

- **Data_RiboTag**

  - **Data_InputandIP**

    This folder contains 29 folders with nested quant.sf files output from salmon transcriptome mapping. Each sample has an Input and an IP folder.

  - **Data_5min**

    This folder contains 8 folders with nested quant.sf files output from salmon transcriptome mapping of just the IP samples for 5 minute timepoint.

  - **Data_30min**

    This folder contains 16 folders with nested quant.sf files output from salmon transcriptome mapping of just the IP samples for all timepoints. 

  - **Data_Nolight**

    This folder contains 11 folders with nested quant.sf files output from salmon transcriptome mapping of RiboTag IP samples for comparing 3 groups, ChR2 + Light, ChR2 + No Light, and Control + Light.
    
    Code for generating these quant.sf files from raw, paired-end fastq files (found on GEO: GSE286381) can be found in Data_RiboTag/RAS_fastq2salmon.R.
    
    The docker containing all software details used in this fastq2salmon pipeline can be found here https://hub.docker.com/r/rubrc/brcpipeline_r4.1.2_base. This platform was developed and executed by The Rockefeller University's Bioinformatics Resource Center (https://www.rockefeller.edu/bioinformatics/)

- **RDS files**

  There are 3 RDS files will need to be read in as detailed in the beginning of the OptoCLIP_Figure4_FigureS4_Rcode_final.R script. These files are located in "RDS_files"

#### **Output (12 plots, 11 source data csv files, 1 RDS file):**

- **Figure 4e, f:** Volcano plots from DESeq comparing control to opto_5min and opto_30min to opto_5min

- **Figure 4g:** Stacked bar plots of DEGs mapping to synaptic, nuclear, or other Gene Ontology (GO) terms 

- **Figure 4h:** Bubble plot showing GO analysis on Opto RiboTag DEGs  

- **Figure 4i:** Violin and box plots showing subcellular localization of Opto RiboTag DEGs 

- **Figure 4j (2 plots, 1 for each timepoint comparison):** Stacked bar plot showing the percentage of FMRP target and non-target transcripts among subsets of Opto RiboTag transcripts 

- **Supplementary Figure 4a:** PCA on OptoRiboTag samples

- **Supplementary Figure 4b:** Heatmap of enrichment of CA1 excitatory marker genes and de-enrichment of other cell type markers in RiboTag IP compared to Input samples

- **Supplementary Figure 4c:** Heatmap of Euclidean distances from 3 groups of RiboTag samples: ChR2 + Light, ChR2+ No Light, and Control + Light

- **Supplementary Figure 4d:** Quantification of Euclidean distances from 3 groups of RiboTag samples: ChR2 + Light, ChR2+ No Light, and Control + Light

- **Supplementary Figure 4e:** Upset plot showing overlap of RiboTag DEGs across different comparisons

- **Supplementary Figure 4f:** Comparison of Opto RiboTag datasets across both time points

- **Source Data for Figure 4e,f**

- **Source Data for Figure 4g**

- **Source Data for Figure 4h**

- **Source Data for Figure 4i and 4i stats**

- **Source Data for Figure 4j**

- **Source Data for Supplementary Figure 4b**

- **Source Data for Supplementary Figure 4d and 4d stats**

- **Source Data for Supplementary Figure 4e**

- **Source Data for Supplementary Figure 4f**

- **RDS object FMRP_RiboTag_merged:** This RDS object will be used in Figures 5, S5, 6, and S7


### Figures 5, S5, 6, and S7- FMRP CLIP score and comparison to RiboTag data

#### **Script:**

- OptoCLIP_Figure5_FigureS5_Figure6_FigureS7_Rcode_final.R

#### **Input data:**

- **RDS files**

  There are 2 RDS files that will need to be read in as detailed in the beginning of thi script. These files are located in "RDS_files"

#### **Output (36 plots, 14 source data csv files):**

- **Figures 5b-d:** Scatter plots of FMRP CLIP binding vs RiboTag for each condition

- **Figures 5e:** Density plots showing the distribution of average CLIP scores for transcripts each condition. 

- **Figures 5f:** Subset of the density plot showing only transcripts that are not stringent in Opto 30 minutes but are stringent in other conditions

- **Figures 5g, h:** Scatter plots of per transcript FMRP CLIP scores comparing different conditions

- **Figures 5i, j:** Volcano plots comparing FMRP CLIP scores across different conditions

- **Figures 5k:** UpSet plot showing overlap between transcripts differentially bound by FMRP across different conditions. 

- **Figures 5l:** Functional classification of transcripts within each intersection category shown in (5k). 

- **Figures 6a, b:** Scatter plots differentially bound FMRP CLIP transcripts versus differentially ribosome associated transcripts 

- **Figures 6c-f_Overlap:** Stacked bar plots showing overlap of transcripts differentially bound by FMRP and differentially ribosome associated. Output graphs were split in Adobe.Illustrator and aligned according to category.

- **Figures 6c-f_GO:** Stacked bar plots showing GO classification of trascripts differentially bound by FMRP and differentially ribosome associated. Output graphs were split in Adobe.Illustrator and aligned according to category.

- **Supplementary Figures 5a-c:** Scatter plots of FMRP CLIP binding from individual samples compared to average RiboTag across the matching. These plots were used to calculate the linear relationship between those variables and calculate a CLIP score per transcript per sample.

- **Supplementary Figure 5d:** Bubble plot showing Gene Ontology (GO) analysis on transcripts differentially bound by FMRP in different conditions.

- **Supplementary Figure 5e:** Comparison of Opto-CLIP datasets across both time points

- **Supplementary Figure 5f:** Enrichment of dendritic FMRP targets (Hale et al 2021) in differentially bound Opto CLIP transcripts 

- **Supplementary Figures 7a, b:** Density ridge plots showing the distribution of RiboTag log2 fold-change subset by transcripts differentially bound by FMRP 

- **Supplementary Figure 7c:** Bubble plot showing Gene Ontology (GO) analysis on transcripts differentially bound by FMRP in different conditions.

- **Supplementary Figures 7d, e:** Violin and box plots showing subcellular localization of transcripts differentially bound by FMRP in different conditions

- **Supplementary Figure 7f:** Upset plot showing overlap of transcripts differentially bound by FMRP and differentially ribosome associated 

- **Source Data for Figures 5b-h**

- **Source Data for Figures 5i,j**

- **Source Data for Figure 5k**

- **Source Data for Figure 5l**

- **Source Data for Figures 6a, b**

- **Source Data for Figures 6c-f**

- **Source Data for Supplementary Figure 5d**

- **Source Data for Supplementary Figure 5f**

- **Source Data for Supplementary Figures 7a, b**

- **Source Data for Supplementary Figures 7a, b stats**

- **Source Data for Supplementary Figure 7c**

- **Source Data for Supplementary Figures 7d, e**

- **Source Data for Supplementary Figures 7d, e stats**

- **Source Data for Supplementary Figure 7f**


## **R_functions**

Reusable analysis functions:

- **run_deseq_contrasts.R**
    - Flexible DESeq2 pipeline with:
        - filtering strategies
        - GO enrichment
        - volcano plots
        - Venn/UpSet analysis


## **RDS_files**

Precomputed reference objects and intermediate datasets:

- mm10_Tx_final.rds and mm10_txdb_lengths.rds
  - mm10 transcript annotations
  - Script showing how this object was made: mm10_gtf_RDS_code.R
  - Input data to make this object: RDS_files/Data_mm10gtf 
  - Used in:
    - OptoCLIP_Figure3_FigureS3_Rcode_final.R
    - OptoCLIP_Figure4_FigureS4_Rcode_final.R
- CH_RT_resdata.rds
  - External dataset from subcellular RiboTag on CA1 microdissected neuropils and cell bodies (Hale et al 2021; PMID: 34939924)
  - Script showing how this object was made: CH_RiboTag_RDS_code.R
  - Input data to make this object: RDS_files/Data_CHale_et_al_2021_RiboTag_salmon
  - Used in:
    - OptoCLIP_Figure4_FigureS4_Rcode_final.R
    - OptoCLIP_Figure5_FigureS5_Figure6_FigureS7_Rcode_final.R
- FMRP_RiboTag_merged.rds
  - Intermediate dataframe created in OptoCLIP_Figure4_FigureS4_Rcode_final.R for simplicity of downstream analyses
  - Used in OptoCLIP_Figure5_FigureS5_Figure6_FigureS7_Rcode_final.RUsed



## **Contact**

Ruth A. Singer, PhD

Postdoctoral Research Associate

Laboratory of Robert B. Darnell

The Rockefeller University

email: <rsinger01@rockefeller.edu>

