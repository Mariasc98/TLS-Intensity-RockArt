# -------------------------------------------------------
# NMF Analysis on 4-Band Raster Dataset (R,G,B,intensity)
# Script developed by Maria Sotomayor Chicote 
# -------------------------------------------------------
import numpy as np
import rasterio
import matplotlib.pyplot as plt
from sklearn.decomposition import NMF

# Helper: percentile stretch (for visualization)
def stretch01(arr, p_low=2, p_high=98):
    lo, hi = np.nanpercentile(arr, [p_low, p_high])
    return np.clip((arr - lo) / (hi - lo + 1e-12), 0, 1)

# Helper: read single-band raster
def read_band(path):
        with rasterio.open(path) as src:
            band = src.read(1).astype(float)
            profile = src.profile
            nodata = src.nodata
        return band, profile, nodata

# NMF function
def nmf_from_geotiff(
    r_path, 
    g_path, 
    b_path, 
    i_path,
    out_path, 
    n_components,
    init,
    max_iter,
    random_state,
    alpha_W,
    alpha_H,
    l1_ratio
):
    """
    Perform Non-negative Matrix Factorization (NMF)
    on four aligned raster bands (R,G,B,intensity)
    Returns:
        nmf_img (n_components, H, W)
        nmf
        W 
        H 
    """

    # Read raster bands
    R, profile, ndR = read_band(r_path)
    G, _, ndG = read_band(g_path)
    B, _, ndB = read_band(b_path)
    I, _, ndI = read_band(i_path)

    # Build a valid-data pixel mask
    def valid(arr, nod):
        m = np.isfinite(arr)
        if nod is not None:
            m &= (arr != nod)
        return m

    mask = valid(R, ndR) & valid(G, ndG) & valid(B, ndB) & valid(I, ndI)

    # Reshape into 3D array [H,W,4] 
    X = np.stack([R, G, B, I], axis=-1)

    # Flatten valid pixels into PCA matrix [N, 4]
    Xv = X[mask]

    # Ensure non-negativity (required for NMF)
    min_val = Xv.min()
    if min_val < 0:
        Xv = Xv - min_val

    # Compute NMF
    nmf = NMF(
        n_components=n_components,
        init=init,          
        max_iter=max_iter,
        random_state=random_state,
        alpha_W=alpha_W,
        alpha_H=alpha_H,
        l1_ratio=l1_ratio
    )

    W = nmf.fit_transform(Xv)   
    H = nmf.components_         

    # Rebuild component rasters 
    H_img, W_img = R.shape
    nmf_img = np.full((n_components, H_img, W_img), np.nan)

    for k in range(n_components):
        comp = np.full((H_img, W_img), np.nan)
        comp[mask] = W[:, k]    
        nmf_img[k] = comp

    # Save NMF to GeoTIFF 
    if out_path is not None:
        profile_out = profile.copy()
        profile_out.update(
            dtype="float32",
            count=n_components,
            nodata=np.nan
        )
        with rasterio.open(out_path, "w", **profile_out) as dst:
            for k in range(n_components):
                dst.write(nmf_img[k].astype("float32"), k+1)

    return nmf_img, nmf, W, H

# Visualization: NMF component images
def plot_nmf_components(nmf_img):
    n_components = nmf_img.shape[0]
    fig, axs = plt.subplots(1, n_components, figsize=(5*n_components, 5))

    for k in range(n_components):
        axs[k].imshow(stretch01(nmf_img[k]), cmap="gray")
        axs[k].set_title(f"NMF Component {k+1}")
        axs[k].axis("off")

    plt.show()

# -----------------------------------------------------------
# User settings and Execution
# -----------------------------------------------------------
if __name__ == "__main__":

    # Edit you file paths here
    R_PATH = "your_path\\raster_R.tif"
    G_PATH = "your_path\\raster_G.tif"
    B_PATH = "your_path\\raster_B.tif"
    I_PATH = "your_path\\raster_Intensity.tif"

    # Select a name file for the output multiband raster. If you don't want the output, write None
    OUTPUT_PATH = None

    # Select parameters. Those selected for the analysis of the case studies within the MA Thesis appear as default
    N_COMPONENTS = 4
    INIT = "nndsvd"
    MAX_ITER = 2000
    RANDOM_STATE = 0
    ALPHA_W = 0.0
    ALPHA_H = 0.0
    L1_RATIO = 0.0


    nmf_img, nmf_model, W, H = nmf_from_geotiff(
        r_path=R_PATH,
        g_path=G_PATH,
        b_path=B_PATH,
        i_path=I_PATH,
        out_path=OUTPUT_PATH,
        n_components=N_COMPONENTS,
        max_iter=MAX_ITER,
        init=INIT,
        random_state=RANDOM_STATE,
        alpha_W=ALPHA_W,
        alpha_H=ALPHA_H,
        l1_ratio=L1_RATIO
    )

    plot_nmf_components(nmf_img)


