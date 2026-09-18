# Code for TLS spatial analysis and subtype characterization in HCC

This repository contains the analysis scripts used to identify, characterize, and subtype tertiary lymphoid structures (TLSs) in hepatocellular carcinoma (HCC) across spatial transcriptomic, CODEX, Visium HD, CosMx, and immunotherapy-related datasets.

The scripts cover:

- TLS detection, quality control, spatial clustering, and subtype classification;
- CODEX-based DBSCAN TLS identification and clustering;
- spatial organization analyses, including Delaunay-based TLS layering, cellular neighborhood analysis, cell–cell proximity, Ripley’s cross-K analysis, and spatial co-occurrence analysis;
- Visium and Visium HD mini-bulk analysis, batch correction, TPM conversion, differential expression, functional enrichment, and cell-type deconvolution using CIBERSORTx and cell2location;
- integration of CODEX and Visium HD TLS profiles;
- tumor-proximity axis analysis and associated gene enrichment analysis;
- prognostic and immunotherapy-related analyses, including TLS signature scoring, survival analysis, logistic regression, and neutrophil composition analysis;
- TLSProfiler feature selection and model construction, as well as HookNet-TLS analysis.

For code-related questions, you can contact the first author, Shuai Li, at 21110880009@m.fudan.edu.cn.
