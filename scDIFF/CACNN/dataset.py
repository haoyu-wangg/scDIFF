import os
import numpy as np
import h5py
import pyBigWig
import scanpy as sc
from anndata import AnnData
from scipy.sparse import csr_matrix
from torch.utils.data import Dataset, DataLoader
import torch

def load_adata(data) -> AnnData:
    '''
    Load single-cell data as AnnData
    data: path to the input h5ad file
    return: AnnData
    '''
    adata = sc.read_h5ad(data)
    if adata.X.max() > 1:
        print("Binarizing data")
        adata.X.data = (adata.X.data > 0).astype(np.float32)
    return adata


class SingleCellDataset(Dataset):
    '''
    Preprocess data and make dataset
    data: AnnData
    genome: reference genome
    bw_path: path to the directory containing BigWig files
    bw_list: list of BigWig filenames
    seq_len: length to extend/trim sequences to
    '''

    def __init__(self, data: AnnData, genome, bw_path, bw_list, seq_len=1344):
        # Filter genes that need to be accessible in at least 1% of cells
        sc.pp.filter_genes(data, min_cells=int(round(0.01 * data.shape[0])))
        self.data = data
        self.seq_len = seq_len
        self.bw_path = bw_path  # Save BigWig path for worker reinitialization
        self.bw_list = bw_list  # Save BigWig filenames for worker reinitialization

        # Load genome
        self.genome = h5py.File(genome, 'r')
        self.bwls = []  # This will be initialized in each worker

        # Save metadata
        self.obs = self.data.obs.copy()
        del self.data.obs
        self.var = self.data.var.copy()
        del self.data.var
        self.X = csr_matrix(self.data.X.T)
        del self.data.X

        if "chr" in self.var.keys():
            self.chroms = self.var["chr"]

    def _initialize_bwls(self, bw_path, bw_list):
        '''
        Initialize BigWig files for the current worker process
        '''
        self.bwls = []
        for bw_file in bw_list:
            bw_file_path = os.path.join(bw_path, bw_file)
            if os.path.exists(bw_file_path):
                try:
                    bw = pyBigWig.open(bw_file_path)
                    self.bwls.append(bw)  # Append the BigWig object to the list
                except Exception as e:
                    print(f"Failed to load BigWig file: {bw_file_path}. Error: {e}")
            else:
                print(f"BigWig file not found: {bw_file_path}")

    def __len__(self):
        return self.X.shape[0]

    def __getitem__(self, index):
        # Reinitialize BigWig files if not initialized (e.g., in each worker process)
        if not self.bwls:
            self._initialize_bwls(self.bw_path, self.bw_list)

        # Retrieve sequence with a center length of seq_len from peak
        chrom, start, end = self.var["chr"].iloc[index], self.var["start"].iloc[index], self.var["end"].iloc[index]
        mid = (int(start) + int(end)) // 2
        left, right = mid - self.seq_len // 2, mid + self.seq_len // 2
        left_pad, right_pad = 0, 0

        if left < 0:
            left_pad = -left
            left = 0
        if right > self.genome[chrom].shape[0]:
            right_pad = right - self.genome[chrom].shape[0]
            right = self.genome[chrom].shape[0]

        seq = self.genome[chrom][left:right]

        # Imputation: Pad the sequence if it is shorter than seq_len
        if len(seq) < self.seq_len:
            seq = np.concatenate((
                np.full(left_pad, -1, dtype=seq.dtype),
                seq,
                np.full(right_pad, -1, dtype=seq.dtype),
            ))

        # Extract epifeatures for the region
        epifeatures = []

        # Define thresholds for different BigWig files

        # thresholds = {
        #     "H3k4me3": 5.0,  # -log10(p-value) > 5
        #     "H3k4me1": 2.0,  
        #     "H3k27ac": 2.0   
        # }

        for bw_file, bw in zip(self.bw_list, self.bwls):
            try:
                # Extract features for the region `left:right`
                feature = bw.values(chrom, left, right, numpy=True)
                # Handle NaN values in BigWig output (replace with 0)
                feature = np.nan_to_num(feature, nan=0.0)

                # Pad features to align with sequence length, if necessary
                if len(feature) < self.seq_len:
                    feature = np.concatenate((
                        np.full(left_pad, 0.0, dtype=feature.dtype),
                        feature,
                        np.full(right_pad, 0.0, dtype=feature.dtype),
                    ))

                epifeatures.append(feature)
            except Exception as e:
                print(f"Error fetching epifeatures for region {chrom}:{left}-{right} in file {bw_file}: {e}")
                epifeatures.append(np.full(self.seq_len, 0.0))  # Add zero-padded feature as a fallback

        # Concatenate all epifeatures into a single NumPy array
        epifeatures = np.stack(epifeatures, axis=0)  # Shape: (num_features, seq_len)

        return seq, epifeatures, self.X[index].toarray().flatten()

