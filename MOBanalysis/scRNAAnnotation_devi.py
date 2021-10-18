import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
import umap

import scvi
from scvi.model import CondSCVI
from scvi.model import DestVI
import torch

import sys
filename1 = sys.argv[1]
filename2 = sys.argv[2]

#matplotlib inline
#pip install --quiet git+https://github.com/yoseflab/scvi-tools@master#egg=scvi-tools


sc_adata = sc.read_h5ad('~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/scAnnotation/OB_sc.h5ad')
G = 3000
sc.pp.filter_genes(sc_adata, min_counts=10)

sc_adata.layers["counts"] = sc_adata.X.copy()

sc.pp.highly_variable_genes(
    sc_adata,
    n_top_genes=G,
    subset=True,
    layer="counts",
    flavor="seurat_v3"
)

sc.pp.normalize_total(sc_adata, target_sum=10e4)
sc.pp.log1p(sc_adata)
sc_adata.raw = sc_adata

st_adata = sc.read_h5ad(filename1)
st_adata.layers["counts"] = st_adata.X.copy()
sc.pp.normalize_total(st_adata, target_sum=10e4)
sc.pp.log1p(st_adata)
st_adata.raw = st_adata

# filter genes to be the same on the spatial data
intersect = np.intersect1d(sc_adata.var_names, st_adata.var_names)
st_adata = st_adata[:, intersect].copy()
sc_adata = sc_adata[:, intersect].copy()
G = len(intersect)


scvi.data.setup_anndata(sc_adata, layer="counts", labels_key="CellType")
sc_model = CondSCVI(sc_adata, weight_obs=True)
sc_model.train(max_epochs=100)
scvi.data.setup_anndata(st_adata, layer="counts")
st_model = DestVI.from_rna_model(st_adata, sc_model)
st_model.train(max_epochs=250)
st_model.get_proportions()
st_adata.obsm["proportions"] = st_model.get_proportions()
prop = st_adata.obsm["proportions"]
prop.to_csv(filename2)

