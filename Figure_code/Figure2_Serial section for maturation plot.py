###### Preliminary maturation-stage assignment in serial TLS sections ######


from pathlib import Path

import numpy as np
import pandas as pd


###### paths and parameters ######

input_path = Path(
    r"E:\01-TLS\01-Omics\05-Serial\02-analysis"
    r"\07-After adjust\07-clean\02-input"
    r"\serial_TLS_composition_final_annotation.csv"
)
output_dir = Path(r"F:\TLS Material\06-Code\Serial_TLS_source_data")
output_dir.mkdir(parents=True, exist_ok=True)

tls_names = [f"TLS{number}" for number in range(1, 8)]
tls_number_map = dict(zip([1, 3, 4, 5, 10, 14, 12], range(1, 8)))
expected_sections = list(range(1, 18))

sfl_cutoff = 15.0
pfl_cutoff = 10.0

# These margins only flag sections for review; they do not change the stage.
sfl_borderline_width = 3.0
pfl_borderline_width = 2.0
minimum_B_lineage_cells = 50


###### read the serial TLS composition ######

data = pd.read_csv(input_path, dtype={"slide": str, "fov": str})

parsed_fov = data["fov"].str.extract(r"^(.+)-Section(\d+)$")
if parsed_fov.isna().any().any():
    raise ValueError("Some serial FOV names do not contain a valid section number.")

source_tls_number = parsed_fov[0].str.extract(r"TLS(\d+)$")[0].astype(int)
data["TLS"] = source_tls_number.map(tls_number_map)
data["TLS"] = data["TLS"].apply(
    lambda value: f"TLS{int(value)}" if pd.notna(value) else np.nan
)
data["section_index"] = parsed_fov[1].astype(int)
data["section_number"] = 1 + (data["section_index"] - 1) * 3
data["z_position_um"] = (data["section_number"] - 1) * 5

data = data[data["TLS"].notna()].copy()
data["fov"] = data["TLS"] + "-Section" + data["section_index"].astype(str)

for tls_name in tls_names:
    one_tls = data[data["TLS"] == tls_name]
    observed_sections = sorted(one_tls["section_index"].unique().tolist())
    if one_tls.shape[0] != 17 or observed_sections != expected_sections:
        raise ValueError(f"{tls_name} does not contain all 17 sections.")


###### calculate B-lineage composition ######

b_lineage_columns = ["B", "B-CD21", "B-GC", "PC"]
for column in b_lineage_columns:
    if column not in data.columns:
        data[column] = 0.0

data["B_lineage_fraction_of_TLS"] = data[b_lineage_columns].sum(axis=1)
data["estimated_B_lineage_cells"] = (
    data["B_lineage_fraction_of_TLS"] * data["final_tls_cell_count"]
)

positive_B_lineage = data["B_lineage_fraction_of_TLS"] > 0

data["B_GC_pct_of_B_lineage"] = 0.0
data["B_CD21_plus_GC_pct_of_B_lineage"] = 0.0

data.loc[positive_B_lineage, "B_GC_pct_of_B_lineage"] = (
    data.loc[positive_B_lineage, "B-GC"] /
    data.loc[positive_B_lineage, "B_lineage_fraction_of_TLS"] * 100
)
data.loc[positive_B_lineage, "B_CD21_plus_GC_pct_of_B_lineage"] = (
    (
        data.loc[positive_B_lineage, "B-CD21"] +
        data.loc[positive_B_lineage, "B-GC"]
    ) /
    data.loc[positive_B_lineage, "B_lineage_fraction_of_TLS"] * 100
)


###### provisional AGG/PFL/SFL assignment ######

data["preliminary_maturation_stage"] = "AGG-TLS"

pfl_rows = (
    (data["B_GC_pct_of_B_lineage"] < sfl_cutoff) &
    (data["B_CD21_plus_GC_pct_of_B_lineage"] >= pfl_cutoff)
)
sfl_rows = data["B_GC_pct_of_B_lineage"] >= sfl_cutoff

data.loc[pfl_rows, "preliminary_maturation_stage"] = "PFL-TLS"
data.loc[sfl_rows, "preliminary_maturation_stage"] = "SFL-TLS"


###### sections requiring spatial and H&E review ######

near_sfl_cutoff = (
    data["B_GC_pct_of_B_lineage"] - sfl_cutoff
).abs() <= sfl_borderline_width

near_pfl_cutoff = (
    (data["B_GC_pct_of_B_lineage"] < sfl_cutoff) &
    (
        data["B_CD21_plus_GC_pct_of_B_lineage"] - pfl_cutoff
    ).abs() <= pfl_borderline_width
)

low_B_lineage = data["estimated_B_lineage_cells"] < minimum_B_lineage_cells

data["manual_review_flag"] = (
    near_sfl_cutoff | near_pfl_cutoff | low_B_lineage
)
data["manual_review_reason"] = ""
data.loc[near_pfl_cutoff, "manual_review_reason"] = "near_AGG_PFL_cutoff"
data.loc[near_sfl_cutoff, "manual_review_reason"] = "near_PFL_SFL_cutoff"
data.loc[low_B_lineage, "manual_review_reason"] = "B_lineage_cells_below_50"

data = data.sort_values(["TLS", "section_index"])

output_columns = [
    "slide", "fov", "TLS", "section_index", "section_number",
    "z_position_um", "final_tls_cell_count", "B", "B-CD21", "B-GC", "PC",
    "B_lineage_fraction_of_TLS", "estimated_B_lineage_cells",
    "B_GC_pct_of_B_lineage", "B_CD21_plus_GC_pct_of_B_lineage",
    "preliminary_maturation_stage", "manual_review_flag",
    "manual_review_reason"
]

data[output_columns].to_csv(
    output_dir / "serial_TLS_maturation_stage_preliminary.csv",
    index=False
)


###### summarize maturation stages across the 17 sections ######

summary_rows = []

for tls_name in tls_names:
    one_tls = data[data["TLS"] == tls_name]
    stage_counts = one_tls["preliminary_maturation_stage"].value_counts()

    summary_rows.append({
        "TLS": tls_name,
        "n_sections": one_tls.shape[0],
        "n_AGG_TLS": int(
            (one_tls["preliminary_maturation_stage"] == "AGG-TLS").sum()
        ),
        "n_PFL_TLS": int(
            (one_tls["preliminary_maturation_stage"] == "PFL-TLS").sum()
        ),
        "n_SFL_TLS": int(
            (one_tls["preliminary_maturation_stage"] == "SFL-TLS").sum()
        ),
        "majority_stage": stage_counts.index[0],
        "majority_fraction": stage_counts.iloc[0] / one_tls.shape[0],
        "n_sections_for_manual_review": int(one_tls["manual_review_flag"].sum())
    })

pd.DataFrame(summary_rows).to_csv(
    output_dir / "serial_TLS_maturation_stage_summary.csv",
    index=False
)
