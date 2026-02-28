# -------------------------------------------------------
# Spectral Index Analysis on 3D Point Clouds 
# Script developed by Maria Sotomayor Chicote 
# -------------------------------------------------------
import numpy as np
import pandas as pd
import open3d as o3d
import pye57
import matplotlib.pyplot as plt

# Load Point Cloud from E57
def load_e57(file_path):
    # Read E57 file
    e57 = pye57.E57(file_path)
    data = e57.read_scan(0,intensity=True, colors=True, ignore_missing_fields=True)  

    # Extract XYZ, Intensity, and Color
    x, y, z = data["cartesianX"], data["cartesianY"], data["cartesianZ"]
    intensity = data.get("intensity", np.zeros_like(x))  # Handle missing intensity
    red, green, blue = data["colorRed"], data["colorGreen"], data["colorBlue"]

    # Normalize colors 
    red, green, blue = np.array(red), np.array(green), np.array(blue)
    

    # Create DataFrame
    df = pd.DataFrame({"X": x, "Y": y, "Z": z, "Red": red, "Green": green, "Blue": blue, "Intensity":intensity})

    return df

# Compute Rock Art Indices
def compute_indices(df):
    df["IRR"] = df["Intensity"] / (df["Red"] + 1e-6)   
    df["IRD"] = df["Intensity"] - df["Red"]
    df["RNI"] = df["Red"] / (df["Intensity"] + 1e-6)  
    df["NRAI"] = ((df["Intensity"] - df["Red"]) - (df["Green"] - df["Blue"])) / \
                 ((df["Intensity"] + df["Red"]) + (df["Green"] + df["Blue"]) + 1e-6)

    # Normalize indices between 0-1
    for col in ["IRR", "IRD", "RNI", "NRAI"]:
        df[col] = (df[col] - df[col].min()) / (df[col].max() - df[col].min())

    return df

# Visualise resulting 3D Point Cloud with Open3D
def visualize_point_cloud(df, index_column):
    # Convert DataFrame to Open3D PointCloud
    pcd = o3d.geometry.PointCloud()
    points = np.vstack((df["X"], df["Y"], df["Z"])).T
    pcd.points = o3d.utility.Vector3dVector(points)

    # Use selected index as color map
    index_values = df[index_column].values
    colors = plt.get_cmap("gray")(index_values)[:, :3]  # Use grayscale for the output point cloud
    pcd.colors = o3d.utility.Vector3dVector(colors)

    # Visualize
    o3d.visualization.draw_geometries([pcd], window_name=f"Point Cloud Visualization - {index_column}")


# -----------------------------------------------------------
# User Settings and Execution
# -----------------------------------------------------------
if __name__ == "__main__":

    # Edit you file path here (e57 format)
    e57_file = "your_path\\pointcloud.e57" 

    # Select the index to run ("IRR","IRD","RNI" or "NRAI")
    index_column = "IRR"

    df = load_e57(e57_file)
    df = compute_indices(df)
    visualize_point_cloud(df, index_column)  
