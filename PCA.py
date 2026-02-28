# -------------------------------------------------------
# PCA Analysis on 4-Band Raster Dataset (R,G,B,intensity)
# Script developed by Maria Sotomayor Chicote 
# -------------------------------------------------------
import numpy as np
import rasterio
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# Helper: percentile stretch (for visualization)
def stretch01(arr, p_low=2, p_high=98):
    lo, hi = np.nanpercentile(arr, [p_low, p_high])
    return np.clip((arr - lo) / (hi - lo + 1e-12), 0, 1)

# Helper: read single-band raster
def read_band(path):
    with rasterio.open(path) as src:
        band = src.read(1)
        profile = src.profile
        nodata = src.nodata
    return band, profile, nodata

# PCA function
def pca_from_geotiff(
    r_path, 
    g_path, 
    b_path, 
    i_path,
    n_components,
    whiten,
    svd_solver,
    random_state,
    out_path
):
    """
    Performs Principal Component Analysis (PCA) on four raster bands (R, G, B, Intensity).
    Returns:
        pca_img: PCA components as 2D rasters
        pca: scikit-learn PCA model
        eigenvalues: array of eigenvalues (explained variance)
        explained_variance_ratio: ratio of total variance explained
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
    

    # Compute PCA 
    pca = PCA(
        n_components=n_components, 
        whiten=whiten, 
        svd_solver=svd_solver,
        random_state=random_state)
    
    Pv = pca.fit_transform(Xv) 

    # Extract eigenvalues
    eigenvalues = pca.explained_variance_
    explained_ratio = pca.explained_variance_ratio_

    print("\n=== PCA Eigenvalues (explained variance) ===")
    for i, ev in enumerate(eigenvalues, 1):
        print(f"Component {i}: {ev:.6f}")

    print("\n=== Explained Variance Ratio ===")
    for i, r in enumerate(explained_ratio, 1):
        print(f"Component {i}: {r*100:.2f}%")

    print(f"\nTotal explained variance: {explained_ratio.sum()*100:.2f}%\n")

    # Rebuild components to raster format
    H, W = R.shape
    pca_img = np.full((n_components, H, W), np.nan, dtype=np.float32)
    for k in range(n_components):
        comp = np.full((H, W), np.nan)
        comp[mask] = Pv[:, k]
        pca_img[k] = comp

    # Save PCA to GeoTIFF
    if out_path is not None:
        profile_out = profile.copy()
        profile_out.update(
            dtype="float32",
            count=n_components,
            nodata=np.nan
        )
        with rasterio.open(out_path, "w", **profile_out) as dst:
            for k in range(n_components):
                dst.write(pca_img[k].astype("float32"), k+1)

    return pca_img, pca, eigenvalues, explained_ratio

# Visualization: PCA component images
def plot_pca_components(pca_img):
    n_components = pca_img.shape[0]
    fig, axs = plt.subplots(1, n_components, figsize=(5*n_components, 5))

    for k in range(n_components):
        axs[k].imshow(stretch01(pca_img[k]), cmap="gray")
        axs[k].set_title(f"PCA Component {k+1}")
        axs[k].axis("off")

    plt.show()

# Visualization: Scree plot (Eigenvalues)
def plot_pca_eigenvalues(eigenvalues, explained_ratio):
    plt.figure(figsize=(6, 4))
    plt.plot(range(1, len(eigenvalues)+1), eigenvalues, marker='o')
    plt.title("PCA Scree Plot (Eigenvalues)")
    plt.xlabel("Component")
    plt.ylabel("Eigenvalue (Explained Variance)")
    plt.grid(True)
    plt.show()

    plt.figure(figsize=(6, 4))
    plt.bar(range(1, len(explained_ratio)+1), explained_ratio*100)
    plt.title("Explained Variance Ratio (%)")
    plt.xlabel("Component")
    plt.ylabel("% of Total Variance")
    plt.show()


# -----------------------------------------------------------
# User Settings and Execution
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
    WHITEN = False
    SVD_SOLVER = "auto"
    RANDOM_STATE = None



    pca_img, pca, eigenvalues, var_ratio = pca_from_geotiff(
        r_path=R_PATH,
        g_path=G_PATH,
        b_path=B_PATH,
        i_path=I_PATH,
        n_components=N_COMPONENTS,
        whiten=WHITEN,
        svd_solver=SVD_SOLVER,
        random_state=RANDOM_STATE,
        out_path=OUTPUT_PATH
        )
    
plot_pca_components(pca_img)
plot_pca_eigenvalues(eigenvalues, var_ratio)
