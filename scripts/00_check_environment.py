import os
import sys
import h5py
import torch
import requests
import scipy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import torch_geometric
import matplotlib.pyplot as plt

print(f"Python version: {sys.version}")
print(f"Torch version: {torch.__version__}")
print(f"Torch Geometric version: {torch_geometric.__version__}")
print(f"Scipy version: {sc.__version__}")
print(f"Numpy version: {np.__version__}")
print(f"Pandas version: {pd.__version__}")
print(f"Seaborn version: {sns.__version__}")  # type: ignore
print(f"Matplotlib version: {plt.matplotlib.__version__}")  # type: ignore
print(f"h5py version: {h5py.__version__}")

print(f"data/raw exists: {os.path.isdir('data/raw')}")
print(f"data/processed exists: {os.path.isdir('data/processed')}")

array = np.array([1, 2, 3], dtype=np.float32)

with h5py.File("test.hdf5", "w") as f:
    f.create_dataset("default", data=array)

with h5py.File("test.hdf5", "r") as f:
    dataset = f["default"][:]  # type: ignore
    print(dataset)

kegg_url = "https://rest.kegg.jp/info/kegg"

response = requests.get(kegg_url)
if response.status_code == 200:
    print("Request successful!")
    print(f"KEGG API response preview: {response.text[:100]}")
else:
    print(f"Request failed with status code: {response.status_code}")
