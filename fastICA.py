# -------------------------------------------------------
# ICA Analysis on 4-Band Raster Dataset (R,G,B,intensity)
# Script developed by Maria Sotomayor Chicote 
# -------------------------------------------------------
import numpy as np
import rasterio
from sklearn.decomposition import FastICA
import matplotlib.pyplot as plt

# Helper: percentile stretch (for visualization)
def stretch01(arr, p_low=2, p_high=98):
    lo, hi = np.nanpercentile(arr, [p_low, p_high])
    out = (arr - lo) / (hi - lo + 1e-12)
    return np.clip(out, 0, 1)

# Helper: read single-band raster
def read_band(path):
        with rasterio.open(path) as src:
            band = src.read(1)  # shape (H, W)
            profile = src.profile
            nodata = src.nodata
        return band, profile, nodata

# ICA function
def fastica_from_geotiff(
    r_path, 
    g_path, 
    b_path, 
    i_path,
    algorithm,
    fun,
    tol,
    whiten,
    max_iter, 
    n_components, 
    random_state,
    out_path
):
    """
    Perform FastICA on four aligned raster bands (R,G,B,intensity).

    Returns:
        components_img (n_components, H, W)
        mixing_matrix
        unmixing_matrix
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

    # Flatten valid pixels into ICA matrix [N, 4] 
    Xv = X[mask]

    # Compute ICA 
    ica = FastICA(
        algorithm=algorithm,
        fun=fun,
        tol=tol,
        n_components=n_components, 
        whiten=whiten,
        random_state=random_state, 
        max_iter=max_iter)
    
    Iv = ica.fit_transform(Xv)  

    # Rebuild components to raster format
    H, W = R.shape
    components_img = np.full((n_components, H, W), np.nan, dtype=np.float32)
    for k in range(n_components):
        comp = np.full((H, W), np.nan, dtype=np.float32)
        comp[mask] = Iv[:, k]
        components_img[k] = comp

    # Save ICA to GeoTIFF 
    if out_path is not None:
        profile_out = profile.copy()
        profile_out.update(
            dtype='float32',
            count=n_components,
            nodata=np.nan
        )
        
        with rasterio.open(out_path, 'w', **profile_out) as dst:
            for k in range(n_components):
                dst.write(components_img[k].astype('float32'), indexes=k+1)

    return components_img, ica.mixing_, ica.components_

# Visualization: ICA components images
def plot_ica_components(components_img, cols=2):
    n_components = components_img.shape[0]
    rows = int(np.ceil(n_components / cols))

    fig, axs = plt.subplots(rows, cols, figsize=(6 * cols, 5 * rows))
    axs = np.asarray(axs).reshape(rows, cols)

    for k in range(n_components):
        ax = axs[k // cols, k % cols]
        ax.imshow(stretch01(components_img[k]), cmap="gray")
        ax.set_title(f"ICA Component {k+1}", fontsize=14)
        ax.axis("off")

    # hide unused axes
    for i in range(n_components, rows * cols):
        axs[i // cols, i % cols].axis("off")

    plt.tight_layout()
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
    ALGORITHM = "parallel"
    FUN = "cube"
    TOL = 1e-5
    N_COMPONENTS = 4
    WHITEN = "unit-variance"
    RANDOM_STATE = 0   
    MAX_ITER=2000



    components_img, mixing, unmixing = fastica_from_geotiff(
        r_path=R_PATH,
        g_path=G_PATH,
        b_path=B_PATH,
        i_path=I_PATH,
        algorithm=ALGORITHM,
        fun=FUN,
        tol=TOL,
        whiten=WHITEN,
        n_components=N_COMPONENTS,
        random_state=RANDOM_STATE,
        max_iter=MAX_ITER,
        out_path=OUTPUT_PATH,
    )

    plot_ica_components(components_img)







