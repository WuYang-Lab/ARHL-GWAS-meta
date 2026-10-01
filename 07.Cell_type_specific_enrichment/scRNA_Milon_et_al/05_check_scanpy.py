# check in scanpy
import scanpy as sc, numpy as np

ad = sc.read_h5ad("Milon_et_al.h5ad")

if "counts" not in ad.layers.keys():
    if ad.raw is not None:
        ad.layers["counts"] = ad.raw.X.copy()
    else:
        ad.layers["counts"] = ad.X.copy()

if "n_counts" not in ad.obs:
    ad.obs["n_counts"] = np.asarray(ad.layers["counts"].sum(axis=1)).ravel()
if "n_genes" not in ad.obs:
    ad.obs["n_genes"] = np.asarray((ad.layers["counts"]>0).sum(axis=1)).ravel()

ad.write_h5ad("Milon_et_al_ready.h5ad", compression="gzip")