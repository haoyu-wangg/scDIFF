

# scDIFF: Automatic Cell Type Annotation of scATAC-seq Data Using a Diffusion-Based Transformer Integrating Histone Modification Information

## Project Overview
**scDIFF**  is a novel method for supervised cell type annotation in single-cell ATAC-seq (scATAC-seq) data, addressing challenges of high dimensionality and sparsity. By integrating genomic sequences and histone modification information around accessibility peaks, scDIFF generates optimized cell embeddings and performs annotation through a diffusion-based Transformer network. It consistently outperforms existing methods on diverse scATAC-seq datasets, providing functional insights through expression and motif enrichment analyses.


## Project Directory Structure
The directory structure of the project is as follows:
```
├── data/                       # Not included in the repository (structure described below)
│   ├── BoneMarrow/             # Example dataset for Bone Marrow
│   │   ├── BoneMarrowA_BoneMarrowB.h5ad    # Single-cell data in h5ad format
│   │   ├── H3k4me1_mm9.bigWig              # Epigenetic feature (ChIP-seq data)
│   │   ├── H3k4me3_mm9.bigWig
│   │   ├── H3k27ac_mm9.bigWig
│   ├── Cortex/                 # Example dataset for Cortex
│       ├── MosP1_MouseBrain(10x).h5ad      # Single-cell data in h5ad format
│       ├── H3k4me1_mm10.bigWig             # Epigenetic feature (ChIP-seq data)
│       ├── H3k4me3_mm10.bigWig
│       ├── H3k27ac_mm10.bigWig
│
├── figures/                    # Directory for architecture and visualization outputs
│   ├── model.svg               # Example model architecture
│
├── output/                     # Example output directory 
│   ├── BoneMarrowA_BoneMarrowB/
│   │   ├── args.txt                 # Configuration and runtime arguments
│   │   ├── CACNN_best_model.pt      # Best model weights
│   │   ├── CACNN_output.h5ad        # CACNN_output data
│   │   ├── CACNN_train.log          # Training logs
│   │   ├── embedding.h5ad           # Embedding results
│   │   ├── result.csv               # Final results in tabular format
│   ├── MosP1_MouseBrain(10x)/       # Another example output (similar structure)
│
├── scDIFF/                     # Core code for the project
│   ├── CACNN/                  # CACNN implementation
│   │   ├── genome/             # Not included in the repository (structure described below)
│   │   │   ├── mm9.fa.h5       # Reference genome for mm9 (mouse genome version 9)
│   │   │   ├── mm10.fa.h5      # Reference genome for mm10 (mouse genome version 10)
│   │   ├── dataset.py          # Data loading and preprocessing
│   │   ├── model.py            # CACNN model architecture
│   │   ├── train.py            # Training script for CACNN
│   │   ├── utils.py            # Utility functions
│   │   ├── ECA_layer.py        # Efficient Channel Attention (ECA) mechanism
│   ├── DIFFomer/               # DIFFomer implementation
│       ├── data_utils.py       # Data handling and preprocessing
│       ├── dataset.py          # Dataset preparation for DIFFomer
│       ├── model.py            # DIFFomer model architecture
│       ├── train.py            # Training script for DIFFomer
│       ├── eval.py             # Evaluation script for DIFFomer
│       ├── metrics.py          # Evaluation metrics
│       ├── logger.py           # Logging utilities
│       ├── parse.py            # Argument parser for training and evaluation
│
├── .gitignore                  # Git ignore rules
├── BoneMarrowA_BoneMarrowB.ipynb  # Example Jupyter Notebook for Bone Marrow analysis
├── MosP1_MouseBrain(10x).ipynb   # Example Jupyter Notebook for Cortex analysis
├── Readme.md                   # Project documentation
├── requirements.txt            # Python dependencies

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

To obtain tissue-specific histone modification information, histone ChIP-seq peak files (in BED format) are aligned to the genome using LiftOver, with '1' indicating the presence and '0' indicating the absence of modification signals for each corresponding protein. These processed histone modification signals are stored in .bigWig format.

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
### Tutorial 2: Cell annotations on datasets cross platforms (MosP1_Cerebellum)
1. Install the required environment according to [Installation](#installation).
2. Create a data folder in the same directory as the 'scDIFF' folder and create subfolder data/Cortex.
3. After downloading the datasets and processing them, place them in data/Cortex.
4. Create a folder genome in the ./scDIFF/CACNN/ directory and download mm10.fa.h5.
5. For more detailed information, run the tutorial [MosP1_Cerebellum.ipynb](./MosP1_Cerebellum.ipynb) for how to do data preprocessing and training.