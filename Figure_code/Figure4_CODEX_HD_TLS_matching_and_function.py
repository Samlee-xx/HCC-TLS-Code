###### Integration of CODEX and Visium HD TLSs ######

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr


###### files ######

analysis_dir = Path(
    r"F:\TLS DATA\03-CODEX-ref_to_HD\03-analysis-final"
    r"\04-Correlation analysis"
)
result_dir = Path(
    r"F:\TLS Material\06-Code\result\Figure4_CODEX_HD_integration"
)




codex_file = (
    result_dir
    / "00-CODEX_154_TLS_clustering"
    / "CODEX_154_TLS_assignment.csv"
)
hd_file = (
    result_dir
    / "01-HD_TLS_cell2location"
    / "HD_TLS_cell2location_mean_proportion.csv"
)
count_file = (
    analysis_dir
    / "01-HD_TLS_region_count_sum"
    / "HD_TLS_region_gene_count_sum.csv"
)
gene_set_file = (
    analysis_dir
    / "04-Correlation"
    / "gene signature.xlsx"
)


shared_cell_types = [
    "B cells",
    "CD4 T cells",
    "CD8 T cells",
    "Cholangiocyte-like cells",
    "Endothelial/lymphatic cells",
    "Macrophages",
    "Neutrophils",
    "Fibroblasts"
]

functional_programs = [
    "CD8 activation",
    "CD40 pathway",
    "CTL pathway",
    "GC B cell",
    "B activation",
    "B cell receptor signaling",
    "Neutrophil chemokines"
]



codex = pd.read_csv(codex_file).set_index("tls_id")
codex.shape[0] == 154

hd = pd.read_csv(hd_file, index_col="tls_id")
candidate_ids = [tls_id for tls_id in hd.index if tls_id in codex.index]

codex_composition = pd.DataFrame({
    "B cells": codex.loc[candidate_ids, "dpct_B"],
    "CD4 T cells": codex.loc[candidate_ids, "dpct_CD4T"],
    "CD8 T cells": codex.loc[candidate_ids, "dpct_CD8T"],
    "Cholangiocyte-like cells": codex.loc[candidate_ids, "dpct_Cholangiocyte"],
    "Endothelial/lymphatic cells": codex.loc[
        candidate_ids, "dpct_Endothelial_Lymphatic"
    ],
    "Macrophages": codex.loc[candidate_ids, "dpct_Macrophage"],
    "Neutrophils": codex.loc[candidate_ids, "dpct_Neutrophil"],
    "Fibroblasts": codex.loc[candidate_ids, "dpct_Fibroblast"]
}).astype(float)

hd_composition = pd.DataFrame({
    "B cells": hd.loc[candidate_ids, "B Cell"],
    "CD4 T cells": hd.loc[candidate_ids, "CD4 T Cell"],
    "CD8 T cells": hd.loc[candidate_ids, "CD8 T Cell"],
    "Cholangiocyte-like cells": hd.loc[candidate_ids, "Bile_Duct"],
    "Endothelial/lymphatic cells": hd.loc[candidate_ids, "Endothelial"],
    "Macrophages": hd.loc[candidate_ids, "Macrophage"],
    "Neutrophils": hd.loc[candidate_ids, "Neutrophil"],
    "Fibroblasts": hd.loc[candidate_ids, "Fibroblast"]
}).astype(float)

hd_composition = hd_composition * 100

complete = ~(
    codex_composition.isna().any(axis=1)
    | hd_composition.isna().any(axis=1)
)
candidate_ids = complete.index[complete].tolist()
codex_composition = codex_composition.loc[candidate_ids]
hd_composition = hd_composition.loc[candidate_ids]


###### match CODEX and Visium HD TLSs ######

matching_rows = []

for tls_id in candidate_ids:
    correlation, p_value = pearsonr(
        codex_composition.loc[tls_id, shared_cell_types],
        hd_composition.loc[tls_id, shared_cell_types]
    )

    matching_rows.append({
        "tls_id": tls_id,
        "tls_subtype": codex.loc[tls_id, "TLS_subtype"],
        "pearson_r": correlation,
        "p_value": p_value,
        "matched": correlation > 0 and p_value < 0.05
    })

matching = pd.DataFrame(matching_rows)
matched = matching.loc[matching["matched"]].copy()


matching.to_csv(
    result_dir / "CODEX_HD_candidate_TLS_pairs.csv",
    index=False
)
matched.to_csv(
    result_dir / "CODEX_HD_matched_TLS_pairs.csv",
    index=False
)

harmonized = pd.concat(
    [
        codex_composition.add_prefix("CODEX_"),
        hd_composition.add_prefix("Visium_HD_")
    ],
    axis=1
)
harmonized.to_csv(
    result_dir / "CODEX_HD_eight_shared_components.csv",
    index_label="tls_id"
)


###### log1p(CPM) expression profiles ######

count_table = pd.read_csv(count_file)
count_table = count_table.dropna(subset=["gene_name"])

gene_count = count_table.set_index("gene_name")[candidate_ids]
gene_count.index = gene_count.index.astype(str).str.upper()
gene_count = gene_count.groupby(level=0).sum().astype(float)

library_size = gene_count.sum(axis=0)
cpm = gene_count.div(library_size, axis=1) * 1_000_000
log_cpm = np.log1p(cpm)


###### functional  programs ######

gene_set_table = pd.read_excel(gene_set_file, sheet_name=0)


gene_sets = {}
for column, program in zip(gene_set_table.columns, functional_programs):
    genes = (
        gene_set_table[column]
        .dropna()
        .astype(str)
        .str.upper()
        .drop_duplicates()
        .tolist()
    )
    gene_sets[program] = [gene for gene in genes if gene in log_cpm.index]

gene_set_summary = pd.DataFrame({
    "functional_program": functional_programs,
    "n_genes_in_Table_S14": [
        gene_set_table[column].dropna().astype(str).str.upper().nunique()
        for column in gene_set_table.columns
    ],
    "n_genes_in_HD_data": [
        len(gene_sets[program]) for program in functional_programs
    ]
})
gene_set_summary.to_csv(
    result_dir / "functional_program_gene_matching.csv",
    index=False
)


###### ssGSEA ######

def calculate_ssgsea(expression, gene_sets, alpha=0.25):
    genes = expression.index.to_numpy()
    values = expression.to_numpy(float)
    ranks = expression.rank(axis=0, method="average").to_numpy(float)

    gene_masks = {
        name: np.isin(genes, one_gene_set)
        for name, one_gene_set in gene_sets.items()
    }

    score = pd.DataFrame(
        index=expression.columns,
        columns=gene_sets.keys(),
        dtype=float
    )

    for sample_number, tls_id in enumerate(expression.columns):
        order = np.argsort(-values[:, sample_number], kind="mergesort")
        rank_weight = ranks[order, sample_number] ** alpha

        for program, gene_mask in gene_masks.items():
            hits = gene_mask[order]
            n_hits = int(hits.sum())
            hit_weight = rank_weight * hits

            running_hit = np.cumsum(hit_weight) / hit_weight.sum()
            running_miss = np.cumsum(~hits) / (len(hits) - n_hits)
            score.loc[tls_id, program] = np.sum(
                running_hit - running_miss
            )

    return score


ssgsea_score = calculate_ssgsea(log_cpm, gene_sets, alpha=0.25)
ssgsea_score.to_csv(
    result_dir / "Visium_HD_TLS_ssGSEA_scores.csv",
    index_label="tls_id"
)


###### Pearson correlations in the 30 matched TNC-TLSs ######

tnc_ids = matched.loc[
    matched["tls_subtype"].eq("TNC_TLS"),
    "tls_id"
].tolist()


correlation_rows = []

for cell_type in shared_cell_types:
    for program in functional_programs:
        correlation, p_value = pearsonr(
            codex_composition.loc[tnc_ids, cell_type],
            ssgsea_score.loc[tnc_ids, program]
        )

        correlation_rows.append({
            "cell_type": cell_type,
            "functional_program": program,
            "n_TNC_TLS": len(tnc_ids),
            "pearson_r": correlation,
            "p_value": p_value
        })

correlation_result = pd.DataFrame(correlation_rows)


###### Benjamini-Hochberg correction across all 56 correlations ######

def bh_adjust(p_values):
    p_values = np.asarray(p_values, dtype=float)
    order = np.argsort(p_values)
    ranked_p = p_values[order]
    adjusted = ranked_p * len(ranked_p) / np.arange(1, len(ranked_p) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]

    answer = np.empty_like(adjusted)
    answer[order] = np.minimum(adjusted, 1)
    return answer


correlation_result["BH_FDR"] = bh_adjust(
    correlation_result["p_value"]
)
correlation_result.to_csv(
    result_dir / "TNC_TLS_CODEX_celltype_vs_HD_function_Pearson.csv",
    index=False
)

correlation_result.pivot(
    index="cell_type",
    columns="functional_program",
    values="pearson_r"
).loc[shared_cell_types, functional_programs].to_csv(
    result_dir / "TNC_TLS_CODEX_HD_Pearson_r_matrix.csv"
)

correlation_result.pivot(
    index="cell_type",
    columns="functional_program",
    values="BH_FDR"
).loc[shared_cell_types, functional_programs].to_csv(
    result_dir / "TNC_TLS_CODEX_HD_BH_FDR_matrix.csv"
)
