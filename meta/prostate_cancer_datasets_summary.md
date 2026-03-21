| dataset_name | source_database | data_modality | accession_id | sample_size | matched_multiomics | data_level | notes |
|-------------|----------------|---------------|--------------|-------------|--------------------|------------|-------|
| TCGA-PRAD RNA-Seq (gene-level) | TCGA/GDC | RNA-Seq (bulk) | TCGA-PRAD | 497 | True | Processed counts + raw (Level 3) | Gene expression counts for 497 prostate adenocarcinoma cases; part of TCGA multi-omics with matched clinical, DNA, CNV, methylation and RPPA data |
| TCGA-PRAD miRNA-Seq | TCGA/GDC | miRNA-Seq | TCGA-PRAD | 494 | True | Processed (RPM counts) | miRNA expression for 494 cases; complements gene expression and other omics data |
| TCGA-PRAD DNA methylation | TCGA/GDC | DNA methylation (Illumina 450K) | TCGA-PRAD | 498 | True | Beta values (Level 3) | CpG-site and gene-level methylation for 498 cases |
| TCGA-PRAD Copy Number Variation | TCGA/GDC | Somatic Copy Number Alteration (SCNA) | TCGA-PRAD | 499 | True | Segment-level CNV (Level 3) | High-level and gene-level copy number data for 499 cases |
| TCGA-PRAD Whole-Exome Sequencing | TCGA/GDC | Somatic mutation (WXS) | TCGA-PRAD | 498 | True | MAF/VCF (Level 3) | Somatic SNV and INDEL calls for 498 cases using exome sequencing; matched normal tissue |
| TCGA-PRAD Whole-Genome Sequencing | TCGA/GDC | Somatic mutation (WGS) | TCGA-PRAD | 414 | True | Raw and processed WGS (Level 2/3) | WGS data for 414 prostate cancers providing structural variants and coverage for copy-number analysis |
| TCGA-PRAD Proteomics (RPPA) | TCGA/GDC | Proteomics (RPPA) | TCGA-PRAD | 352 | True | Protein abundance (relative) | RPPA data for 352 cases covering ~200 proteins; complements gene expression and genomic data |
| TCGA-PRAD ATAC-Seq | TCGA/GDC | Chromatin accessibility (ATAC-Seq) | TCGA-PRAD | 26 | True | Peak calls and raw counts | ATAC-Seq on 26 prostate tumors for chromatin accessibility profiles |
| TCGA-PRAD Clinical & Phenotype | TCGA/GDC | Clinical/Phenotype | TCGA-PRAD | 500 | True | Clinical tables and survival data | Clinical characteristics, treatments, and outcomes for 500 prostate cancer patients; used for survival analyses |
| TCGA-PRAD Imaging (TCIA) | TCIA | Imaging (Radiology & Histopathology) | TCIA-TCGA-PRAD | 469 | True | DICOM radiology files & digitized slides | Imaging data for 469 cases including CT, MRI and PET scans; 1,172 SVS slides and 16,790 imaging files total |
| Single-patient prostate spatial transcriptomics | EGA | Spatial transcriptomics | EGAD00001007921 | 2 | False | BAM & counts | Visium spatial transcriptomics from two tissue sections of a single untreated prostate cancer patient |
| CRUK-ICGC Prostatectomy batches 4-6 WGS | ICGC/EGA | Whole Genome Sequencing | EGAD00001003225 | 221 | True | BAM/FASTQ | Batch 4-6 WGS from 111 men (221 samples) exploring multifocal prostate tumors |
| ICGC WGS 51 tumors (10 patients) | ICGC/EGA | Whole Genome Sequencing | EGAD00001000891 | 62 | True | BAM/FASTQ | WGS sequencing of 51 tumors and matched normals from 10 patients |
| ICGC WGS 20 men primary tumors | ICGC/EGA | Whole Genome Sequencing | EGAD00001000892 | 40 | True | BAM/FASTQ | WGS of primary tumors and matched blood from 20 patients |
| ICGC WGS pelvic lymph node metastases | ICGC/EGA | Whole Genome Sequencing | EGAD00001002002 | 20 | True | BAM/FASTQ | WGS of pelvic lymph node metastases with matched blood from 10 patients |
| ICGC Normal prostatectomy WGS | ICGC/EGA | Whole Genome Sequencing | EGAD00001004125 | 71 | True | BAM/FASTQ | Normal prostate tissue WGS from cancer and non-cancer individuals |
| ICGC multi-region WGS 3 patients | ICGC/EGA | Whole Genome Sequencing | EGAD00001000689 | 18 | True | BAM/FASTQ | Deep WGS of tumor and normal samples from three men |
| ICGC SNP6 genotype data | ICGC/EGA | Genotyping / SNP array | EGAD00010000498 | 18 | True | Raw SNP calls | Affymetrix SNP6 genotyping data for CNV estimation |
| CRUK-ICGC DNA methylation sequencing | ICGC/EGA | DNA methylation (capture sequencing) | EGAD00001010184 | 376 | True | BAM & coverage tables | DNA methylation sequencing for tumor-control pairs |
| SU2C mCRPC RNA-Seq | SU2C/PCF via EGA | RNA-Seq (bulk) | EGAD00001008991 | 118 | False | FASTQ/BAM | RNA-Seq of metastatic castration-resistant prostate cancer |
| SU2C metastatic WES (Cell 2015) | SU2C/PCF via cBioPortal | Whole-exome sequencing | prad_su2c_2015 | 300 | True | MAF/VCF | WES of metastatic prostate cancer tumor/normal pairs |
| mCRPC integrative genomics (PNAS 2019) | SU2C/PCF via cBioPortal | Genomics & Transcriptomics | prad_su2c_2019 | 429 | True | MAF/VCF + RNA counts | Integrative genomic and transcriptomic analysis |
| GSE147250 RNA-Seq (TP53/RB1 loss) | GEO | RNA-Seq (bulk) | GSE147250 | 162 | False | Counts and FPKM | RNA-Seq exploring TP53/RB1 loss |
| GSE80609 RNA-Seq | GEO | RNA-Seq (bulk) | GSE80609 | 45 | False | Counts and FPKM | RNA-Seq across prostate cancer progression |
| GSE158593 RNA-Seq (NE vs AdCa) | GEO | RNA-Seq (bulk) | GSE158593 | 48 | False | Counts | RNA-Seq of LuCaP xenografts |
| PXD013422 tissue proteomics | PRIDE | Proteomics (LC-MS/MS) | PXD013422 | nan | False | Raw mass spectra & quantification | Proteomics identifying ~1,904 proteins |
| PXD022005 Lysine succinylation | PRIDE | Proteomics (succinylome) | PXD022005 | nan | False | Raw mass spectra | Succinylome analysis |
| PXD028651 Serum proteomics | PRIDE | Proteomics (SWATH-MS) | PXD028651 | nan | False | Raw SWATH-MS spectra | Serum proteomics |
| LuCaP PDX proteome & phosphoproteome | MassIVE | Proteomics & phosphoproteomics | MSV000092139 | 48 | False | Raw MS/MS and quantified peptides | Proteome + phosphoproteome profiling |
| DepMap CRISPR gene dependency | DepMap | Functional dependency (CRISPR) | DepMap 24Q2 | 1320 | False | Gene effect scores | CRISPR knockout screens |
| DepMap RNAi gene dependency | DepMap | Functional dependency (RNAi) | DEMETER2 | 712 | False | Gene effect scores | RNAi knockdown screens |
| Project SCORE CRISPR dependency | Sanger / Project SCORE | Functional dependency (CRISPR) | ProjectSCORE | 324 | False | Gene dependency scores | CRISPR screens |
| GDSC drug sensitivity | GDSC | Drug response (viability) | GDSC1/GDSC2 | 1000 | False | IC50/area-under-curve | Drug screening |
| CCLE mRNA expression for prostate lines | CCLE | RNA-Seq (cell lines) | CCLE_Expression | 8 | True | TPM/FPKM counts | Expression for prostate cell lines |
| CCLE miRNA expression | CCLE | miRNA expression | CCLE_miRNA | 7 | True | Counts | miRNA expression |
| CCLE copy number | CCLE | Copy number | CCLE_CNV | 8 | True | Log2 ratios | Copy number profiles |
| CCLE DNA methylation | CCLE | DNA methylation | CCLE_Methylation | 6 | True | Beta values | Methylation profiles |
| CCLE RPPA proteomics | CCLE | Proteomics (RPPA) | CCLE_RPPA | 7 | True | Relative abundance | Protein profiling |
| ENCODE LNCaP ATAC-Seq and ChIP-Seq | ENCODE | ATAC-Seq & ChIP-Seq | LNCaP datasets | 5 | False | Aligned reads & peak calls | Chromatin accessibility and TF binding |
| Roadmap Epigenomics prostate tissue | Roadmap Epigenomics | Histone modification & DNA methylation | Roadmap E091/E092 | 2 | False | Signal tracks & peak calls | Epigenomic maps |