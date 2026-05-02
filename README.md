# Galectin-1 Invasive Subpopulations

## Data Availability
The raw sequencing data and annotated matricies used in this project are availible at GEO under GSE328755

## Data Acquisition
Patient-derived GBM cell lines, HK177, HK408, and GS54 were used in these studies. Cells were encapulated in HA-based hydrogels mimicking the stiffness of tumor and peritumoral brain tissue. 
### Sequencing 
Cells were extracted from hydrogels sorted using FACS. Samples were prepared using 10X Genomics’ Chromium Single Cell 3ʹ v3 kit. 1000-1500 cells per condition were sequenced on an Illumina NovaSeq SR100 chip at a depth of 50,000 reads/cell. Annotated count matrices were generated from the raw fastq files using Cell Ranger. Genes were aligned to the reference genome, NIH GRCh38. 
### Migration
Migration images were aquired every 4 hours using Incucyte S3 with the 4X objective using phase contrast imaging. 
### Mouse Section Imaging
Mice with orthotropic xenografts of the HK408 line were used to validate protein expression. Images were aquired using Leica SP8 laser confocal microscope with the 20× objective lense.

## Environments
Environments and dependencies were handled with conda. For MAST analysis which require R implementation, environments/environment_r.yml defines a seperate environment containing rpy packages. For image analysis environments/image.yml was used. For scRNAseq analysis environments/scrna.yml was used.

## Workflow
### Sequening 
QC for HK177 and GS54 was conducted first to remove low quality reads. To analzye the differentially expressed genes, analysis through MAST was conducted. MAST datafiles were uploaded to QIAGEN Ingenuity Pathway Analysis (IPA) to analyze differentially expressed pathways. 
### Migration
Following migration image acquisition, images were aligned using crop_align_incucyte.py.
