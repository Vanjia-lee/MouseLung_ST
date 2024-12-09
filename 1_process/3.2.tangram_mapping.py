################################################################################
# mapping single cell to spatial data by Tangram with the mode of `cells` and
# 1544 marker genes of 89 cell-type-related clusters from single-cell lung atlas 
# as training genes.

################################################################################

import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import torch
import diopy
sys.path.append('./')  # uncomment for local import
import tangram as tg

##################################################################################
# 0.Set Directory
projectDir = '~/project/lung_10xST'
refDir = os.path.join(projectDir, "data/reference_sc/GSE165063_Development_2021/processing/")
outdir = os.path.join(projectDir, 'results/Tangram_mapping/tangram_origQC_training_genes_mode_cells')
isExists = os.path.exists(outdir)
if not isExists:
    os.mkdir(outdir)

###################################################################################
# 1.load sc data
path = os.path.join(refDir, 'SC_ref_E12_P0_new_subcelltype.h5')
ad_sc =  diopy.input.read_h5(path, assay_name = "RAW_COUNTS")

###################################################################################
# 2. Mapping with 89 cell type from SC lung atlas
### choose 1544 marker genes of the single-cell lung atlas
df_genes = pd.read_csv(os.path.join(refDir, "SC_clusters89_DEGs.csv"), index_col = 0)
markers = np.reshape(df_genes.values, (-1, ))
markers = list(markers)
len(markers)
#1544

stages = ["P0_S1","E185_S3","E175_S2","E165_S2","E155_S4","E145_S4","E135_S3","E125_S2"]
for sta in stages:
    print("Now is mapping analysis for " + sta)
    path = os.path.join(projectDir, 'results/Tangram_mapping/ST_h5', sta + '_afterQC_rmHbMtRp.h5')
    ad_sp = diopy.input.read_h5(path, assay_name='Spatial')
    tg.pp_adatas(ad_sc, ad_sp, genes = markers)
    ad_map = tg.map_cells_to_space(adata_sc=ad_sc, adata_sp=ad_sp,
                                   mode='cells', device='cpu')

    cell_spot = pd.DataFrame(ad_map.X)
    cell_spot.index = ad_sc.obs.index
    cell_spot.columns = ad_sp.obs.index
    cell_spot.to_csv(os.path.join(outdir,sta + "_tangram_mode_cells_training_genes_ct89.csv"))

    # map with celltype89
    tg.project_cell_annotations(ad_map, ad_sp, annotation="new_subcelltype")
    df = ad_sp.obsm["tangram_ct_pred"]
    df.to_csv(os.path.join(outdir, sta + "_tangram_ct89_train_mode_cells_pred.csv"))
    perc=0
    df = df.clip(df.quantile(perc), df.quantile(1 - perc), axis=1)
    df = (df - df.min()) / (df.max() - df.min())
    df.to_csv(os.path.join(outdir, sta + "_tangram_ct89_train_mode_cells_pred_scale.csv"))

    print("Mapping analys for " + sta + "is completed.")
