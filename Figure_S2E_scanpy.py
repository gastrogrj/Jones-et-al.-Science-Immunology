"""Fig S2E Scanpy rerun, condensed from the recovered S2_python_mouse notebook.

Run Export_Figure_S2E_input.R first. Input is the same post-QC colon cell set as
Figure 2. This compares downstream pipelines, not independent biological samples.
The source notebook plotted scaled/regressed values in its dot plot. Here the
dot plot uses the preserved log-normalised RNA matrix; the correction can change
that plot. QC regression covariates are the original count-derived Seurat values,
not proportions calculated from already log-transformed expression. Rerun and
check the five clusters, marker overlap and correlations before replacing S2E.
"""
from pathlib import Path
import importlib.metadata as metadata
import json
import numpy as np
import pandas as pd
import scipy.io
import scanpy as sc
import anndata as ad
import harmonypy as hm
import matplotlib.pyplot as plt

ROOT = Path("output/Figure_S2E")
INPUT = ROOT / "input"
ROOT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 7, "svg.fonttype": "none", "pdf.fonttype": 42})
genes = pd.read_csv(INPUT / "genes.csv")["gene"].astype(str)
cells = pd.read_csv(INPUT / "cells.csv", index_col="cell")
logcounts = scipy.io.mmread(INPUT / "logcounts.mtx").T.tocsr()
counts = scipy.io.mmread(INPUT / "counts.mtx").T.tocsr()
if logcounts.shape != (len(cells), len(genes)):
    raise ValueError("Exported matrix dimensions and labels differ")
rna = ad.AnnData(X=logcounts, obs=cells.copy(), var=pd.DataFrame(index=genes))
rna.layers["counts"] = counts
rna.uns["log1p"] = {"base": None}
required = ["percent.mt", "percent.ribo", "percent.Xist", "sample"]
if not set(required).issubset(rna.obs.columns):
    raise ValueError("Missing original covariates: " + str(set(required) - set(rna.obs.columns)))
sc.pp.highly_variable_genes(rna, n_top_genes=4000, flavor="seurat")
adata = rna[:, rna.var.highly_variable].copy()
sc.pp.regress_out(adata, ["percent.mt", "percent.ribo", "percent.Xist"])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, random_state=0)
harmony = hm.run_harmony(adata.obsm["X_pca"].astype(np.float64), adata.obs, ["sample"], random_state=0)
coords = np.asarray(harmony.Z_corr)
if coords.shape[0] != adata.n_obs:
    coords = coords.T
if coords.shape[0] != adata.n_obs:
    raise ValueError("Harmony coordinates do not match the cell count")
adata.obsm["X_pca_harmony"] = coords
sc.pp.neighbors(adata, n_pcs=30, use_rep="X_pca_harmony", random_state=0)
sc.tl.umap(adata, random_state=0)
# Retain the original selected Louvain resolution; omit the exploratory scan.
sc.tl.louvain(adata, resolution=0.48, key_added="cluster", random_state=0)
if adata.obs["cluster"].nunique() != 5:
    raise ValueError("Expected five clusters; check package versions and source inputs before relabelling")
rna.obs["cluster"] = adata.obs["cluster"]
rna.obsm["X_umap"] = adata.obsm["X_umap"]
sc.pl.umap(rna, color="cluster", size=8, frameon=False, show=False)
for ext in ("svg", "pdf"):
    plt.savefig(ROOT / f"FigS2E_scanpy_UMAP.{ext}", bbox_inches="tight")
plt.close("all")
marker_genes = ["Cxcl10", "Cxcl9", "Acod1", "Gbp5", "Gbp2", "Slpi", "Cers6", "Ifi203", "Irf1", "Inhba"]
dot = sc.pl.dotplot(rna, marker_genes, groupby="cluster", use_raw=False, return_fig=True,
                    colorbar_title="Mean RNA expression\nlog1p(counts per 10,000)")
dot.savefig(ROOT / "FigS2E_scanpy_marker_dotplot.svg")
dot.savefig(ROOT / "FigS2E_scanpy_marker_dotplot.pdf")
# Cluster 4 was the IAM population in the original Scanpy notebook. Verify its
# marker expression after rerunning; cluster indices are not biological labels.
sc.tl.rank_genes_groups(rna, groupby="cluster", groups=["4"], reference="rest",
                        method="wilcoxon", use_raw=False)
markers = sc.get.rank_genes_groups_df(rna, group="4").rename(columns={
    "names": "gene", "logfoldchanges": "avg_log2FC", "pvals": "p_val", "pvals_adj": "p_val_adj"})
markers.to_csv(ROOT / "Scanpy_IAM_all_markers.csv", index=False)
top = markers.query("avg_log2FC > 0 and p_val_adj < 0.05").sort_values(
    ["avg_log2FC", "gene"], ascending=[False, True]).head(200)
top.to_csv(ROOT / "Scanpy_IAM_top200.csv", index=False)
rna.obs.to_csv(ROOT / "Scanpy_cell_annotations.csv")
pd.DataFrame(rna.obsm["X_umap"], index=rna.obs_names, columns=["UMAP_1", "UMAP_2"]).to_csv(
    ROOT / "Scanpy_UMAP_coordinates.csv")
versions = {}
for package in ("scanpy", "anndata", "numpy", "pandas", "scipy", "harmonypy", "matplotlib", "louvain", "igraph"):
    try:
        versions[package] = metadata.version(package)
    except metadata.PackageNotFoundError:
        versions[package] = "not reported"
(ROOT / "python_versions.json").write_text(json.dumps(versions, indent=2))
print("Run Figure_S2E_concordance.R next; compare the computed overlap and correlations with the panel.")
