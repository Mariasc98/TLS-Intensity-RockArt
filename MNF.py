# -------------------------------------------------------
# MNF Analysis on 4-Band Raster Dataset (R,G,B,intensity)
# Script developed by Maria Sotomayor Chicote 
# -------------------------------------------------------
import numpy as np
import rasterio
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

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

# MNF function
def mnf_from_geotiff(
    r_path, 
    g_path, 
    b_path, 
    i_path,
    out_path,
    n_components
):
    """
    Performs Minimum Noise Fraction (MNF) transform on 4 aligned rasters (R,G,B,intensity)
    Returns:
        mnf_img  = MNF components as rasters (n_components, H, W)
        eigenvalues = MNF eigenvalues (signal-to-noise)
        noise_whitening_matrix
        pca
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

    # Flatten valid pixels into MNF matrix [N, 4] 
    Xv = X[mask]   

    # Compute MNF

    # -----------------------------------------------------------
    # Step 1: Estimate noise — local difference between pixels
    # -----------------------------------------------------------
    # noise = X - local_mean

    noise = Xv[1:] - Xv[:-1]                    # differences
    noise_cov = np.cov(noise, rowvar=False)     # noise covariance

    # -----------------------------------------------------------
    # Step 2: Noise whitening
    # -----------------------------------------------------------
    eigvals, eigvecs = np.linalg.eigh(noise_cov)

    # Avoid division by zero
    eigvals[eigvals < 1e-12] = 1e-12

    D_noise = np.diag(1.0 / np.sqrt(eigvals))
    W_noise = eigvecs @ D_noise @ eigvecs.T     # whitening matrix

    # Whiten the data
    X_white = Xv @ W_noise.T

    # -----------------------------------------------------------
    # Step 3: PCA on noise-whitened data
    # -----------------------------------------------------------
    pca = PCA(n_components=n_components)
    Mv = pca.fit_transform(X_white)

    mnf_eigenvalues = pca.explained_variance_

    print("\n=== MNF Eigenvalues (signal-to-noise) ===")
    for i, ev in enumerate(mnf_eigenvalues, 1):
        print(f"MNF Component {i}: {ev:.6f}")

    # -----------------------------------------------------------
    # Step 4: Rebuild MNF components 
    # -----------------------------------------------------------
    H, W = R.shape
    mnf_img = np.full((n_components, H, W), np.nan)
    for k in range(n_components):
        comp = np.full((H, W), np.nan)
        comp[mask] = Mv[:, k]
        mnf_img[k] = comp

    # Save MNF to GeoTIFF
    if out_path is not None:
        profile_out = profile.copy()
        profile_out.update(
            dtype="float32",
            count=n_components,
            nodata=np.nan
        )
        with rasterio.open(out_path, "w", **profile_out) as dst:
            for k in range(n_components):
                dst.write(mnf_img[k].astype("float32"), k+1)

    return mnf_img, mnf_eigenvalues, W_noise, pca

# Visualization: MNF component images
def plot_mnf_components(mnf_img):
    n_components = mnf_img.shape[0]
    fig, axs = plt.subplots(1, n_components, figsize=(5*n_components, 5))

    for k in range(n_components):
        axs[k].imshow(stretch01(mnf_img[k]), cmap="gray")
        axs[k].set_title(f"MNF Component {k+1}")
        axs[k].axis("off")

    plt.show()

# -----------------------------------------------------------
# User settings and Execution
# -----------------------------------------------------------
if __name__ == "__main__":

    # Edit your file paths here
    R_PATH = "your_path\\raster_R.tif"
    G_PATH = "your_path\\raster_G.tif"
    B_PATH = "your_path\\raster_B.tif"
    I_PATH = "your_path\\raster_Intensity.tif"

    # Select a name file for the output multiband raster. If you don't want the output, write None
    OUTPUT_PATH = None

    # Select parameters. Those selected for the analysis of the case studies within the MA Thesis appear as default
    N_COMPONENTS = 4

    
    mnf_img, eigenvalues, W_noise, pca_model = mnf_from_geotiff(
        r_path=R_PATH,
        g_path=G_PATH,
        b_path=B_PATH,
        i_path=I_PATH,
        out_path=OUTPUT_PATH,
        n_components=N_COMPONENTS
    )

    plot_mnf_components(mnf_img)










