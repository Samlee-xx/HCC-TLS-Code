
###### This script used to split into small image to perform segmentation

import os
import tifffile
import numpy as np

input_dir = r'...\ICI_img_normalized'
output_dir = r'...\ICI_img_split'
os.makedirs(output_dir, exist_ok=True)

MAX_SIZE = 2048 

def save_patch(patch, base_name, suffix):
    out_path = os.path.join(output_dir, f"{base_name}_{suffix}.tiff")
    tifffile.imwrite(out_path, patch.astype(np.float32))

for filename in os.listdir(input_dir):
    if not filename.lower().endswith(('.tif', '.tiff')):
        continue

    input_path = os.path.join(input_dir, filename)
    base_name = os.path.splitext(filename)[0]
    img = tifffile.imread(input_path)  # shape: (C, H, W)

    if img.ndim != 3:
        print(f"Ignore:{filename}")
        continue

    C, H, W = img.shape
    print(f"Processing: {filename}, shape: {img.shape}")

    if H <= MAX_SIZE and W <= MAX_SIZE:
        save_patch(img, base_name, "full")
        continue

    if H > MAX_SIZE and W > MAX_SIZE:
        img_q1 = img[:, :H//2, :W//2]
        img_q2 = img[:, :H//2, W//2:]
        img_q3 = img[:, H//2:, :W//2]
        img_q4 = img[:, H//2:, W//2:]
        save_patch(img_q1, base_name, "q1")
        save_patch(img_q2, base_name, "q2")
        save_patch(img_q3, base_name, "q3")
        save_patch(img_q4, base_name, "q4")
    elif H > MAX_SIZE:
        img_top = img[:, :H//2, :]
        img_bottom = img[:, H//2:, :]
        save_patch(img_top, base_name, "top")
        save_patch(img_bottom, base_name, "bottom")
    elif W > MAX_SIZE:
        img_left = img[:, :, :W//2]
        img_right = img[:, :, W//2:]
        save_patch(img_left, base_name, "left")
        save_patch(img_right, base_name, "right")

    print(f"Completion:{base_name}")
    
###### Split channels into single channel


import tifffile
import numpy as np
import os

input_dir = r'...\ICI_img_split'
output_root = r'...\Split_channels'
os.makedirs(output_root, exist_ok=True)

channel_names = [
    "CD45", "VEGF", "HLA-DR", "SMA", "CD15", "LYVE1", "CD69", "CD3", "Ki-67", "CD163", "PanKeratin",
    "CXCR5", "PDL1", "CXCR6", "TCF1", "TOX", "Tim3", "FoxP3", "CD4", "CX3CR1", "Ecadherin", "CD68", 
    "Tbet", "CD20", "CD8a", "Eomes", "Arginase1", "PD1", "CD204", "GranzymeB", "CD39", "Collagen", 
    "CD103", "CD56", "CD38", "CD45RO", "CD33", "CD34", "CK7", "DNA1", "DNA2", "HH3", "ICSK1", "ICSK2", "ICSK3"
]

for filename in os.listdir(input_dir):
    if filename.lower().endswith('.tif') or filename.lower().endswith('.tiff'):
        input_path = os.path.join(input_dir, filename)
        print(f"Processing: {filename}")
        img = tifffile.imread(input_path)
        if img.ndim != 3 or img.shape[0] != len(channel_names):
            print(f"Ignore: {filename})
            continue
        basename = os.path.splitext(filename)[0]
        output_dir = os.path.join(output_root, basename)
        os.makedirs(output_dir, exist_ok=True)
        for i in range(img.shape[0]):
            channel_img = img[i]
            channel_name = channel_names[i]
            out_path = os.path.join(output_dir, f"{channel_name}.tiff")
            tifffile.imwrite(out_path, channel_img)

        print(f"Completion: {filename} into {output_dir}")

###### Segmentation 

import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import skimage.io as io
import xarray as xr
from alpineer import io_utils

from ark.segmentation import marker_quantification, segmentation_utils
from ark.utils import (deepcell_service_utils, example_dataset,
                       plot_utils)
                       
base_dir = r"...\Split_channels"
result_path = r'...\Segmentation_Out'

cell_table_dir = os.path.join(result_path, "cell_table")
deepcell_input_dir = os.path.join(result_path, "deepcell_input")
deepcell_output_dir = os.path.join(result_path, "deepcell_output")
deepcell_visualization_dir = os.path.join(result_path, "deepcell_visualization")

# create directories if do not exist
for directory in [cell_table_dir, deepcell_input_dir, deepcell_output_dir, deepcell_visualization_dir]:
    if not os.path.exists(directory):
        os.makedirs(directory)

# input data file path
tiff_dir = base_dir
fovs = io_utils.list_folders(tiff_dir)
# nuclear channel name(s) (or nucs = None)
nucs = ['HH3','DNA2']
# membrane channel name(s) (or mems = None)
mems = ['CD3','CD4','CD8a','CD15','CD20','CD45','CD45RO']
deepcell_service_utils.generate_deepcell_input(
    deepcell_input_dir,
    tiff_dir,
    nucs,
    mems,
    fovs,
    img_sub_folder=None
)

rescale_factor = 2.0
deepcell_service_utils.create_deepcell_output(deepcell_input_dir, 
                                              deepcell_output_dir, 
                                              fovs=fovs, 
                                              scale=rescale_factor)

segmentation_utils.save_segmentation_labels(
    segmentation_dir=deepcell_output_dir,
    data_dir=deepcell_input_dir,
    output_dir=deepcell_visualization_dir,
    fovs=io_utils.remove_file_extensions(fovs),
    channels=['nuclear_channel', 'membrane_channel']
)

import numpy as np
np.float = float
# set to True to add nuclear cell properties to the expression matrix
nuclear_counts = False

# set to True to bypass expensive cell property calculations
# only cell label, size, and centroid will be extracted if True
fast_extraction = False

# now extract the segmented imaging data to create normalized and transformed expression matrices
# note that if you're loading your own dataset, please make sure all the imaging data is in the same folder
# with each fov given its own folder and all fovs having the same channels
#The bactch size means run rov numbers simutaneously
cell_table_size_normalized, cell_table_arcsinh_transformed = \
    marker_quantification.generate_cell_table(segmentation_dir=deepcell_output_dir,
                                              tiff_dir=tiff_dir,
                                              img_sub_folder=None,
                                              fovs=fovs,
                                              batch_size=5,
                                              nuclear_counts=nuclear_counts,
                                              fast_extraction=fast_extraction)








