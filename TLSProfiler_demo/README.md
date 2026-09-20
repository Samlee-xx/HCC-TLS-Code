# TLSProfiler inference demo

TLSProfiler identifies and classifies tertiary lymphoid structures (TLSs) in five-channel mIHC images. This repository provides the inference code used in the accompanying manuscript and an example image set for testing the workflow.

The analysis includes:

1. validation of DAPI, CD45, CD20, CD3, and CD66B TIFF images;
2. whole-cell segmentation with DeepCell Mesmer;
3. marker-based assignment of B cells, T cells, neutrophils, other immune cells, and other cells;
4. DBSCAN-based TLS detection and polygon reconstruction; and
5. BF-TLS/TNC-TLS classification with the final logistic-regression model.

This release covers inference. Model-training scripts and training data are not included.

## Contact

Shuai Li  
PhD, Fudan University  
Email: 21110880009@m.fudan.edu.cn

## Repository contents

- `app.py`: Streamlit interface.
- `run_tlsprofiler.py`: command-line inference workflow.
- `segment_and_quantify.py`: image validation, DeepCell segmentation, and marker quantification.
- `call_mihc_cell_types.py`: marker thresholding and cell-type assignment.
- `detect_and_classify_tls.py`: TLS detection, polygon reconstruction, and subtype prediction.
- `visualize_results.py`: input, threshold, cell-type, and TLS-region figures.
- `model/tlsprofiler_log1p_BTN.json`: final logistic-regression coefficients and model settings.
- `demo/input/`: five-channel 8000 x 8000 example.

## Requirements

- Python 3.10
- Docker Desktop or Docker Engine
- A valid DeepCell access token for Mesmer segmentation

The Python package versions used for testing are listed in `requirements.txt`.

## Setup

Create and activate a Python 3.10 environment, install the dependencies, and download the DeepCell image:

```powershell
conda create -n tlsprofiler python=3.10 -y
conda activate tlsprofiler
python -m pip install -r requirements.txt
docker pull vanvalenlab/deepcell-applications@sha256:fbfba01a02827ef785cda974e56407d331a35daa2038601c15a4ea680ef951ce
```

Set the DeepCell token in the current PowerShell session:

```powershell
$env:DEEPCELL_ACCESS_TOKEN = '<your-token>'
```

## Browser interface

With the Python environment activated, run:

```powershell
.\run_streamlit.ps1
```

Open the local URL printed by Streamlit and upload the five TIFF files from `demo/input/` to the corresponding channels. The default image resolution is `0.5` microns per pixel, and the default minimum TLS size is `300` cells.

The interface shows the uploaded channels, marker thresholds, cell-type assignments, retained TLS regions, subtype probabilities, and result tables. Only TLSs meeting the selected minimum cell count are shown in the final TLS figure. The complete detection table is retained in the output for inspection.

The 8000 x 8000 example may take several minutes when DeepCell runs on CPU. A warning about slower peak finding for images larger than 5000 x 5000 pixels is expected and does not indicate a failed run.

## Command line

```powershell
python .\run_tlsprofiler.py `
  --sample-id 768241_demo `
  --dapi .\demo\input\DAPI.tif `
  --cd45 .\demo\input\CD45.tif `
  --cd20 .\demo\input\CD20.tif `
  --cd3 .\demo\input\CD3.tif `
  --cd66b .\demo\input\CD66B.tif `
  --minimum-tls-cells 300 `
  --output-dir .\runs\768241_demo
```

The command creates:

- `01-segmentation/deepcell/`: DeepCell input and segmentation output;
- `01-segmentation/cell_table_arcsinh.csv.gz`: quantified cells;
- `02-cell_types/`: marker thresholds and cell-type assignments;
- `03-tls/tls_candidates_and_subtypes.csv`: all detected TLS regions;
- `03-tls/retained_tls_and_subtypes.csv`: TLSs retained at the selected minimum cell count;
- `03-tls/tls_polygons.csv`: reconstructed TLS polygons in WKT format;
- `04-visualizations/`: input, threshold, cell-type, and TLS-region figures; and
- `TLSProfiler_summary.json`: summary of the run.

## Default analysis settings

The example was tested with `image_mpp = 0.5`, an arcsinh cofactor of `0.01`, DBSCAN `eps = 100`, DBSCAN `min_samples = 50`, an alpha-shape radius of `150`, a minimum B/T percentage of `10%`, and a minimum polygon size of `100` cells. The final display threshold defaults to `300` cells and can be changed in the interface or with `--minimum-tls-cells`.

The model JSON records the feature order, transformation, coefficients, and probability cutoff used for BF-TLS/TNC-TLS prediction.

## Using a precomputed segmentation mask

To start from an existing whole-cell label mask, provide it with `--mask`. This skips DeepCell segmentation while retaining the same cell-type and TLS-classification steps.
