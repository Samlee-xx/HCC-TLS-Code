###### CODEX cell segmentation with ARK and Mesmer ######

import os
import subprocess

from alpineer import io_utils
from ark.segmentation import marker_quantification, segmentation_utils
from ark.utils import deepcell_service_utils


# TLS-containing regions were delineated using H&E staining and the
# fluorescence signals of CD20 and CD3. Individual marker channels were
# exported from QuPath v0.4.3 as TIFF files before this step.
# Before run Deepcell for segmentation in Docker, you have to acquire the token in website https://www.deepcell.org/

qupath_root = r"F:\TLS DATA\01-CODEX\02-Segmentation\01 Qupath out"
segmentation_root = r"F:\TLS DATA\01-CODEX\02-Segmentation\02 Segmentation"

docker_image = "vanvalenlab/deepcell-applications:latest"
image_mpp = "0.5"

slide_names = sorted([
    x for x in os.listdir(qupath_root)
    if os.path.isdir(os.path.join(qupath_root, x))
])


###### Mesmer segmentation and cell-table generation ######

for slide_name in slide_names:
    print("Processing:", slide_name)

    tiff_dir = os.path.join(qupath_root, slide_name)
    result_dir = os.path.join(segmentation_root, slide_name)

    cell_table_dir = os.path.join(result_dir, "cell_table")
    deepcell_input_dir = os.path.join(result_dir, "deepcell_input")
    deepcell_output_dir = os.path.join(result_dir, "deepcell_output")
    deepcell_visualization_dir = os.path.join(result_dir, "deepcell_visualization")

    for one_dir in [cell_table_dir,
                    deepcell_input_dir,
                    deepcell_output_dir,
                    deepcell_visualization_dir]:
        os.makedirs(one_dir, exist_ok=True)

    fovs = io_utils.list_folders(tiff_dir)

    nuclear_channels = ["DAPI"]

    membrane_channels = [
        "CD45", "CD20", "CD3e", "CD8", "CD4",
        "CD66b", "Pan-Cytokeratin", "SMA", "CD31", "CD11c"
    ]

    deepcell_service_utils.generate_deepcell_input(
        deepcell_input_dir,
        tiff_dir,
        nuclear_channels,
        membrane_channels,
        fovs,
        img_sub_folder=None
    )

    fov_names = io_utils.remove_file_extensions(fovs)

    # Both CODEX batches were run with Mesmer in Docker.
    for fov_name in fov_names:
        for compartment in ["whole-cell", "nuclear"]:
            if compartment == "whole-cell":
                output_name = fov_name + "_whole_cell.tiff"
            else:
                output_name = fov_name + "_nuclear.tiff"

            docker_command = [
                "docker", "run", "--rm", "-i",
                "-e", "DEEPCELL_ACCESS_TOKEN",
                "-v", result_dir + ":/data",
                docker_image,
                "mesmer",
                "--nuclear-image", "/data/deepcell_input/" + fov_name + ".tiff",
                "--nuclear-channel", "0",
                "--membrane-image", "/data/deepcell_input/" + fov_name + ".tiff",
                "--membrane-channel", "1",
                "--output-directory", "/data/deepcell_output",
                "--output-name", output_name,
                "--compartment", compartment,
                "--image-mpp", image_mpp,
                "--batch-size", "1"
            ]

            subprocess.run(docker_command, check=True)

    segmentation_utils.save_segmentation_labels(
        segmentation_dir=deepcell_output_dir,
        data_dir=deepcell_input_dir,
        output_dir=deepcell_visualization_dir,
        fovs=fov_names,
        channels=["nuclear_channel", "membrane_channel"]
    )

    cell_table_size_normalized, cell_table_arcsinh_transformed = \
        marker_quantification.generate_cell_table(
            segmentation_dir=deepcell_output_dir,
            tiff_dir=tiff_dir,
            img_sub_folder=None,
            fovs=fovs,
            batch_size=5,
            nuclear_counts=False,
            fast_extraction=False
        )

    cell_table_size_normalized.to_csv(
        os.path.join(cell_table_dir, "cell_table_size_normalized.csv"),
        index=False
    )

    cell_table_arcsinh_transformed.to_csv(
        os.path.join(cell_table_dir, "cell_table_arcsinh_transformed.csv"),
        index=False
    )
