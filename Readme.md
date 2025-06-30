

# scDIFF: Automatic cell type annotation using scATAC-seq data by incorporating bulk-level genomic and epigenomic information in diffusive transformers 

## Project Overview
We proposed **scDIFF**, which integrates bulk-level genomic and epigenomic information in a diffusive transformer to annotate cell types from scATATC-seq data. The performance of scDIFF is comparable or better than state-of-the-art methods in all comparisons. Our scDIFF can eliminate batch effects between different datasets, which enabled it to integrate datasets from different sequencing platforms together. With optimized epigenomic profiles, scDIFF present superior annotation performances, indicating it has a great potential to facilitate analysis of scATAC-seq data. These results demonstrate that scDIFF provides both accurate cell type annotations and meaningful insights into the gene regulatory mechanisms at single cell level.



## Project Directory Structure

The directory structure of the project is as follows:

```
├── data/                       # Not included in the repository (structure described below)
│   ├── BoneMarrow/             # Example dataset for Bone Marrow
│   │   ├── BoneMarrowA_BoneMarrowB.h5ad    # Single-cell data in h5ad format
│   │   ├── H3k4me1_mm9.bigWig              # Epigenetic feature (ChIP-seq data)
│   │   ├── H3k4me3_mm9.bigWig
│   │   ├── H3k27ac_mm9.bigWig
│   └── Cortex/                 # Example dataset for Cortex
│       ├── MosP1_MouseBrain(10x).h5ad      # Single-cell data in h5ad format
│       ├── H3k4me1_mm10.bigWig             # Epigenetic feature (ChIP-seq data)
│       ├── H3k4me3_mm10.bigWig
│       └── H3k27ac_mm10.bigWig
│
├── dataPreprocessing/          # Data preprocessing scripts and utilities
│   ├── create_h5ad.ipynb       # Create h5ad files
│   ├── create_pairs.ipynb      # Creating dataset pairs
│   ├── intersect_union.sh      # Shell script for BED intersection and union operations
│   ├── mm10bed_to_hg19bed.sh   # From mm10 BED to hg19 BED format
│   ├── mm10bed_to_mm10bigwig.sh # From mm10 BED to mm10 BigWig format
│   └── mm10bed_to_mm9bigwig.sh  # From mm10 BED to mm9 BigWig format
│
├── drawFigures/                # Figure generation and visualization scripts
│   ├── Bar.ipynb               
│   ├── Coverage_Plot.R         
│   ├── Line_chart.ipynb        
│   ├── Signac.R                
│   ├── comparison_plot.ipynb   
│   ├── confusion_matrix.ipynb  
│   └── evaluateCrossMethods.ipynb 
│
├── figures/                    
│   └── model.svg               
│
├── output/                     # Example output directory  
│   ├── BoneMarrowA_BoneMarrowB/
│   │   ├── args.txt                 
│   │   ├── CACNN_best_model.pt      # Best model weights
│   │   ├── CACNN_output.h5ad        # CACNN_output data
│   │   ├── CACNN_train.log          
│   │   ├── embedding.h5ad           # Embedding 
│   │   ├── result.csv               # performance
│   │   ├── BoneMarrowA_BoneMarrowB_ref0_query1.csv # Reference-query labels
│   │   ├── edge.txt                 
│   │   └── model.pkl                
│   └── MosP1_MouseBrain(10x)/       
│       ├── CACNN_best_model.pt      
│       ├── CACNN_output.h5ad        
│       ├── CACNN_train.log          
│       ├── MosP1_MouseBrain(10x)_ref0_query1.csv 
│       ├── args.txt                 
│       ├── edge.txt                 
│       ├── embedding.h5ad           
│       ├── model.pkl                
│       └── result.csv               
│
├── scDIFF/                     # Core code for the project
│   ├── CACNN/                   
│   │   ├── genome/             # Not included in the repository (structure described below)
│   │   │   ├── mm9.fa.h5       
│   │   │   └── mm10.fa.h5      
│   │   ├── dataset.py          
│   │   ├── model.py            # CACNN model architecture
│   │   ├── train.py            
│   │   ├── utils.py            
│   │   └── ECA_layer.py        
│   └── DIFFormer/              
│       ├── data_utils.py       
│       ├── dataset.py          
│       ├── model.py            # DIFFomer model architecture
│       ├── train.py            
│       ├── eval.py             
│       ├── metrics.py          
│       ├── logger.py           
│       └── parse.py            
│
├── .gitignore                  
├── BoneMarrowA_BoneMarrowB.ipynb  
├── MosP1_MouseBrain(10x).ipynb   
├── Readme.md                   
└── requirements.txt            
```



## Installation
To reproduce scDIFF, we suggest first create a conda environment by:

```
conda create -n scDIFF python=3.9
conda activate scDIFF
```

and then run the following code to install the required package:

```
pip install -r requirements.txt
```

## Usage
### data preprocessing
In order to run scDIFF, we need to first create an anndata object from the raw data. The .h5ad file should have cells as obs and peaks as var. The var must include at least three columns: chr, start, and end, which indicate the genomic region of each peak. Additionally, the obs should contain two columns: Batch and CellType (reference data), where Batch distinguishes between reference and query data, and CellType indicates the true label of the cell.

For the sequence information, we use the corresponding reference genome files provided by the UCSC Genome Browser as input. To set this up, create a folder named genome in the ./scDIFF/CACNN/ directory and download the reference genome file (e.g., mm9.fa.h5) into this folder. This file will be used to extract sequence features for the peaks.

To obtain tissue-specific histone modification information, histone ChIP-seq peak files (in BED format) are aligned to the genome using LiftOver, with '1' indicating the presence and '0' indicating the absence of modification signals for each corresponding protein. . In cases where multi-source reference sets derive from different tissues, we consider two strategies for handling tissue-specific histone modifications: Union strategy: if a nucleotide locates within any tissue's ChIP-seq peak, it is marked as '1', otherwise marked as '0'; Intersection strategy: if a nucleotide locates within all tissues' ChIP-seq peaks, it is marked as '1', otherwise marked as '0'. All ChIP-seq peaks are directly downloaded from the ensemble database.These processed histone modification signals are stored in .bigWig format.

For optimal performance, we recommend filtering out peaks that are accessible in less than 1% of cells. Processed data should be organized in the data directory according to the tissue type, as described in the directory structure above.

### Stage 1: Feature extraction
```
# Stage 1: Feature extraction
cd scDIFF/CACNN

python train.py -i ../../data/BoneMarrow/BoneMarrowA_BoneMarrowB.h5ad \ 
            -g mm9 \  
            -o ../../output/BoneMarrowA_BoneMarrowB \ 
            --tissue BoneMarrow \ 
            --bw_list H3k4me1_mm9.bigWig H3k4me3_mm9.bigWig H3k27ac_mm9.bigWig \ 
            --epifeature_dim 3 
```
Running the above command will generate three output files in the output path:
+ CACNN_train.log: recording logs during training
+ CACNN_best_model.pt: storing the model weights with the best AUC score during training
+ CACNN_output.h5ad: an anndata file storing the embedding extracted by CACNN.

### Stage 2: Classifier via DIFFormer
```
# Stage 2: Classifier via DIFFormer
%cd ../DIFFormer
python train.py --data_dir ../../output/BoneMarrowA_BoneMarrowB/CACNN_output.h5ad \
                --train_name_list BoneMarrowA --test_name BoneMarrowB \
                --save_path ../../output \
                --save_name BoneMarrowA_BoneMarrowB \
                --patience 50 \
                --dropout 0.0 
```
Running the above command will generate three output files in the output path:
+ model.pkl: storing the model weights with the best valid loss during training.
+ embedding.h5ad: an anndata file storing the embedding extracted by GraphTransformer. And .obs['Pred'] saves the results of the prediction.

## Tutorial
### Tutorial 1: Cell annotations within samples (BoneMarrowA_BoneMarrowB)
1. Install the required environment according to [Installation](#installation).
2. Create a data folder in the same directory as the 'scDIFF' folder and create subfolder data/BoneMarrow.
3. After downloading the datasets and processing them, place them in data/BoneMarrow.
4. Create a folder genome in the ./scDIFF/CACNN/ directory and download mm9.fa.h5.
5. For more detailed information, run the tutorial [BoneMarrowA_BoneMarrowB.ipynb](./BoneMarrowA_BoneMarrowB.ipynb) for how to do data preprocessing and training.
### Tutorial 2: Cell annotations on datasets cross platforms (MosP1_MouseBrain(10x))
1. Install the required environment according to [Installation](#installation).
2. Create a data folder in the same directory as the 'scDIFF' folder and create subfolder data/Cortex.
3. After downloading the datasets and processing them, place them in data/Cortex.
4. Create a folder genome in the ./scDIFF/CACNN/ directory and download mm10.fa.h5.
5. For more detailed information, run the tutorial [MosP1_MouseBrain(10x).ipynb](./MosP1_MouseBrain(10x).ipynb) for how to do data preprocessing and training.

# Datasets
### scATAC-seq Datasets Summary

| Dataset | Technology | Reference Genome | GEO Accession |
|---------|------------|------------------|---------------|
| BoneMarrowA/B | sciATAC-seq | mm9 | GSE111586 |
| LungA/B | sciATAC-seq | mm9 | GSE111586 |
| LargeIntestineA/B | sciATAC-seq | mm9 | GSE111586 |
| WholeBrainA/B | sciATAC-seq | mm9 | GSE111586 |
| SmallIntestine | sciATAC-seq | mm9 | GSE111586 |
| Cerebellum | sciATAC-seq | mm9 | GSE111586 |
| PreFrontalCortex | sciATAC-seq | mm9 | GSE111586 |
| MosA1/A2 | snATAC-seq | GRCm38 | GSE126724 |
| MosM1/M2 | snATAC-seq | GRCm38 | GSE126724 |
| MosP1/P2 | snATAC-seq | GRCm38 | GSE126724 |
| MouseBrain(10x) | 10x sequencing | mm10 | https://support.10xgenomics.com/single-cell-atac/datasets/1.1.0/atac_v1_adult_brain_fresh_5k |
| NormalCortex | 10x sequencing | mm10 | https://www.10xgenomics.com/datasets/fresh-cortex-from-adult-mouse-brain-p-50-1-standard-1-2-0 |
| KidneyA-E | snATAC-seq | hg38 | GSE185948 |

### ChIP-seq Datasets Overview

| Tissue | Type | Accession ID | Developmental Stage |
|--------|------|--------------|-------------------|
| Bone marrow | H3K4me3 | ENCFF366WUW | adult (2 months) |
| Bone marrow | H3K4me1 | ENCFF675QTV | adult (2 months) |
| Bone marrow | H3K27ac | ENCFF009GVZ | adult (2 months) |
| Cerebellum | H3K4me3 | ENCFF026XUO | adult (2 months) |
| Cerebellum | H3K4me1 | ENCFF934THC | adult (2 months) |
| Cerebellum | H3K27ac | ENCFF763RNO | adult (2 months) |
| Forebrain | H3K4me3 | ENCFF160SCR | postnatal (0 days) |
| Forebrain | H3K4me1 | ENCFF937JHP | postnatal (0 days) |
| Forebrain | H3K27ac | ENCFF442GIT | postnatal (0 days) |
| Intestine | H3K4me3 | ENCFF609BZI | postnatal (0 days) |
| Intestine | H3K4me1 | ENCFF523POG | postnatal (0 days) |
| Intestine | H3K27ac | ENCFF182RQN | postnatal (0 days) |
| Lung | H3K4me3 | ENCFF082TBM | postnatal (0 days) |
| Lung | H3K4me1 | ENCFF453ZKR | postnatal (0 days) |
| Lung | H3K27ac | ENCFF535ZRH | postnatal (0 days) |
| Large intestine | H3K4me3 | ENCFF609BZI | postnatal (0 days) |
| Large intestine | H3K4me1 | ENCFF523POG | postnatal (0 days) |
| Large intestine | H3K27ac | ENCFF182RQN | postnatal (0 days) |
| Cortex | H3K4me3 | ENCFF218SDV | adult (2 months) |
| Cortex | H3K4me1 | ENCFF407TGJ | adult (2 months) |
| Cortex | H3K27ac | ENCFF958PXM | adult (2 months) |
| Kidney | H3K4me3 | ENCFF349XHZ | adult (50 years) |
| Kidney | H3K4me1 | ENCFF955PYL | adult (50 years) |
| Kidney | H3K27ac | ENCFF418XDA | adult (50 years) |

**Note**: Accession ID refers to unique identifiers assigned by the ENCODE database for each dataset.

### Data pre-processing
#### scATAC-seq was pre-processed using the following steps:
- `create_h5ad.ipynb`
- `Create_pairs.ipynb`, with the generated files then placed into the ./data/corresponding tissue directories

####  ChIP-seq was pre-processed using the following steps:

- `intersect_union.sh`
- `mm10bed_to_hg19bed.sh`
- `mm10bed_to_mm9bigwig.sh`
- `mm10bed_to_mm10bigwig.sh`

Scripts were selected based on the genome version of the modification information used, with the generated BigWig files named according to modification information and then placed into the ./data/corresponding tissue directories.

# Figure Reproduction Instructions

## Main Figures

### Figure 2
- **Figure 2(A)**: Generated using `./drawFigures/evaluateCrossMethods.ipynb`; `./drawFigures/comparison_plot.ipynb`
- **Figure 2(B)(C)**: Generated using `./drawFigures/confusion_matrix.ipynb`
- **Figure 2(D)(E)**: Generated using corresponding code in Visualization sections of Tutorial 1 and Tutorial 2

### Figure 3
Generated using `./drawFigures/Line_chart.ipynb`

### Figure 4
- **Figure 4(A)**: Generated using `./drawFigures/evaluateCrossMethods.ipynb`; `./drawFigures/comparison_plot.ipynb`
- **Figure 4(B)**: Generated using corresponding code in Visualization sections of Tutorial 1 and Tutorial 2
- **Figure 4(C)**: Generated using corresponding code in Visualization sections of Tutorial 1 and Tutorial 2

### Figure 5
Generated using `./drawFigures/evaluateCrossMethods.ipynb`; `./drawFigures/Bar.ipynb`

### Figure 6
- **Figure 6(A)**: Generated using `./drawFigures/evaluateCrossMethods.ipynb`; `./drawFigures/Bar.ipynb`
- **Figure 6(B-F)**: Generated using corresponding code in Visualization sections of Tutorial 1 and Tutorial 2

### Figure 7
- **Figure 7(A)**: Generated using `./drawFigures/evaluateCrossMethods.ipynb`; `./drawFigures/Bar.ipynb`
- **Figure 7(B-H)**: Generated using corresponding code in Visualization sections of Tutorial 1 and Tutorial 2

### Figure 8
- **Figure 8(A)**: Generated using `./drawFigures/Coverage_Plot.R`
- **Figure 8(B)**: Generated using SNPsea software package
- **Figure 8(C)(D)**: Generated using `./drawFigures/Signac.R`

---

## Supplementary Figures

### Figure S1-S2
**Figure S1(A-C), S2(A-C)**: Generated using `./drawFigures/evaluateCrossMethods.ipynb`; `./drawFigures/comparison_plot.ipynb`

### Figure S3-S5
- **Figure S3(A-C)**: Generated using `./drawFigures/evaluateCrossMethods.ipynb`; `./drawFigures/Bar.ipynb`
- **Figure S4(A-C)**: Generated using `./drawFigures/evaluateCrossMethods.ipynb`; `./drawFigures/Bar.ipynb`
- **Figure S5(A-C)**: Generated using `./drawFigures/evaluateCrossMethods.ipynb`; `./drawFigures/Bar.ipynb`

### Figure S6-S7
- **Figure S6(A-D)**: Generated using `./drawFigures/Coverage_Plot.R`
- **Figure S7(A-E)**: Generated using `./drawFigures/Signac.R`