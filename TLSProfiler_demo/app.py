from __future__ import annotations

import io
import os
import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st

from run_tlsprofiler import run_pipeline


st.set_page_config(page_title="TLSProfiler demo", layout="wide")
st.title("TLSProfiler")
st.caption("Five-channel mIHC inference demo")

sample_id = st.text_input("Sample ID", value="sample_001")
columns = st.columns(5)
uploads = {}
for column, marker in zip(columns, ["DAPI", "CD45", "CD20", "CD3", "CD66B"]):
    with column:
        uploads[marker] = st.file_uploader(marker, type=["tif", "tiff"], key=marker)

with st.expander("Advanced settings"):
    setting_columns = st.columns(4)
    with setting_columns[0]:
        mask_upload = st.file_uploader("Optional whole-cell label mask", type=["tif", "tiff"], key="mask")
    with setting_columns[1]:
        cofactor = st.number_input(
            "Arcsinh cofactor", min_value=0.0001, value=0.01, step=0.005, format="%.4f"
        )
    with setting_columns[2]:
        image_mpp = st.number_input(
            "Image resolution (microns per pixel)", min_value=0.05, value=0.5, step=0.05
        )
    with setting_columns[3]:
        minimum_tls_cells = st.number_input(
            "Minimum cells per TLS", min_value=100, value=300, step=50
        )

if not os.environ.get("DEEPCELL_ACCESS_TOKEN", "").strip() and mask_upload is None:
    st.warning("DEEPCELL_ACCESS_TOKEN is not set. Add it to the PowerShell session before running without a mask.")

run_clicked = st.button("Run TLSProfiler", type="primary")
if run_clicked:
    st.session_state.pop("tlsprofiler_result", None)
    if any(upload is None for upload in uploads.values()):
        st.error("Upload all five marker channels.")
    else:
        try:
            with st.spinner("Running segmentation, cell calling, TLS detection, and visualization..."):
                with tempfile.TemporaryDirectory(prefix="tlsprofiler_") as tmp:
                    work = Path(tmp)
                    channel_paths = {}
                    for marker, upload in uploads.items():
                        path = work / f"{marker}.tif"
                        path.write_bytes(upload.getvalue())
                        channel_paths[marker] = path

                    mask_path = None
                    if mask_upload is not None:
                        mask_path = work / "whole_cell_mask.tif"
                        mask_path.write_bytes(mask_upload.getvalue())

                    output = work / "output"
                    summary = run_pipeline(
                        sample_id=sample_id,
                        channel_paths=channel_paths,
                        output_dir=output,
                        mask_path=mask_path,
                        cofactor=cofactor,
                        image_mpp=image_mpp,
                        minimum_tls_cells=int(minimum_tls_cells),
                        segmentation_mode="docker",
                    )
                    tls_path = Path(summary["tls_table"])
                    thresholds_path = output / "02-cell_types" / "marker_thresholds.csv"
                    st.session_state["tlsprofiler_result"] = {
                        "summary": summary,
                        "tls_bytes": tls_path.read_bytes(),
                        "thresholds_bytes": thresholds_path.read_bytes(),
                        "visualizations": {
                            name: Path(path).read_bytes() for name, path in summary["visualizations"].items()
                        },
                    }
        except Exception as error:
            st.error("TLSProfiler did not complete.")
            st.exception(error)

result = st.session_state.get("tlsprofiler_result")
if result is not None:
    summary = result["summary"]
    tls = pd.read_csv(io.BytesIO(result["tls_bytes"]))
    thresholds = pd.read_csv(io.BytesIO(result["thresholds_bytes"]))
    metric_columns = st.columns(3)
    metric_columns[0].metric("Retained TLS", summary["retained_tls_count"])
    metric_columns[1].metric("BF-TLS", summary["BF_TLS_count"])
    metric_columns[2].metric("TNC-TLS", summary["TNC_TLS_count"])

    tab_input, tab_thresholds, tab_cells, tab_tls, tab_data = st.tabs(
        ["Input images", "Marker thresholds", "Cell types", "TLS regions", "Results"]
    )
    with tab_input:
        st.image(result["visualizations"]["input_channels"], use_container_width=True)
        st.download_button(
            "Download input-channel figure",
            result["visualizations"]["input_channels"],
            file_name=f"{sample_id}_01_input_channels.png",
            mime="image/png",
            key="download_input_channels",
        )

    with tab_thresholds:
        st.image(result["visualizations"]["marker_thresholds"], use_container_width=True)
        shown_columns = [
            "marker",
            "threshold",
            "threshold_method",
            "threshold_pool",
            "target_positive_pct_all_retained_cells",
        ]
        st.dataframe(thresholds[[column for column in shown_columns if column in thresholds]], use_container_width=True)
        st.download_button(
            "Download threshold figure",
            result["visualizations"]["marker_thresholds"],
            file_name=f"{sample_id}_02_marker_thresholds.png",
            mime="image/png",
            key="download_threshold_figure",
        )

    with tab_cells:
        st.image(result["visualizations"]["cell_types"], use_container_width=True)
        st.download_button(
            "Download cell-type figure",
            result["visualizations"]["cell_types"],
            file_name=f"{sample_id}_03_cell_types.png",
            mime="image/png",
            key="download_celltype_figure",
        )

    with tab_tls:
        st.image(result["visualizations"]["tls_regions"], use_container_width=True)
        if tls.empty:
            st.info(
                f"No TLS reached the selected minimum of {summary['minimum_tls_cells']:,} cells."
            )
        st.download_button(
            "Download TLS-region figure",
            result["visualizations"]["tls_regions"],
            file_name=f"{sample_id}_04_tls_regions.png",
            mime="image/png",
            key="download_tls_figure",
        )

    with tab_data:
        if tls.empty:
            st.info(
                f"No TLS reached the selected minimum of {summary['minimum_tls_cells']:,} cells."
            )
        else:
            st.dataframe(tls, use_container_width=True)
        download_columns = st.columns(2)
        with download_columns[0]:
            st.download_button(
                "Download TLS results",
                result["tls_bytes"],
                file_name=f"{sample_id}_TLSProfiler_results.csv",
                mime="text/csv",
                key="download_tls_results",
            )
        with download_columns[1]:
            st.download_button(
                "Download marker thresholds",
                result["thresholds_bytes"],
                file_name=f"{sample_id}_marker_thresholds.csv",
                mime="text/csv",
                key="download_threshold_table",
            )
