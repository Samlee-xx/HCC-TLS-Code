###### Serial-section CODEX and z-axis TLS subtype stability analysis ######

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


###### paths and parameters ######

input_dir = Path(
    r"E:\01-TLS\01-Omics\05-Serial\02-analysis"
    r"\07-After adjust\07-clean\02-input"
)
output_dir = Path(r"F:\TLS Material\06-Code\Serial_TLS_source_data")
output_dir.mkdir(parents=True, exist_ok=True)

training_path = input_dir / "training_TLS_composition.csv"
serial_path = input_dir / "serial_TLS_composition_final_annotation.csv"

random_seed = 20260710
positive_class = "TNC-TLS"
class_order = ["BF-TLS", "TNC-TLS"]
expected_sections = list(range(1, 18))

tls_names = [f"TLS{number}" for number in range(1, 8)]
tls_number_map = dict(zip([1, 3, 4, 5, 10, 14, 12], range(1, 8)))

features = [
    "B", "CD4T", "CD8T", "Cholangiocyte_like", "DC",
    "Endothelial", "Fibroblast", "HEV", "Lymphatic",
    "Macrophage", "Neutrophil"
]


###### model definitions ######

def make_logistic_model():
    return Pipeline([
        ("scale", StandardScaler()),
        (
            "classifier",
            LogisticRegression(
                penalty="l2",
                class_weight="balanced",
                solver="liblinear",
                max_iter=5000,
                random_state=random_seed
            )
        )
    ])


def probability_of_TNC(model, feature_table):
    class_names = list(model.classes_)
    return model.predict_proba(feature_table)[:, class_names.index(positive_class)]


###### prepare the 667-TLS training cohort ######

training = pd.read_csv(
    training_path,
    dtype={"tls_id": str, "slide": str, "fov": str}
)
training = training[training["tls_subtype"].isin(class_order)].copy()
training = training.reset_index(drop=True)


training_features = pd.DataFrame({
    "B": training["pct_B"],
    "CD4T": training["pct_CD4T"],
    "CD8T": training["pct_CD8T"],
    "Cholangiocyte_like": training["pct_Cholangiocyte"],
    "DC": training["pct_DC"],
    "Endothelial": training["pct_Endothelial"],
    "Fibroblast": training["pct_Fibroblast"],
    "HEV": training["pct_HEV"],
    "Lymphatic": training["pct_Lymphatic"],
    "Macrophage": training["pct_Macrophage"],
    "Neutrophil": training["pct_Neutrophil"]
})[features].fillna(0).astype(float)

training_labels = training["tls_subtype"].astype(str)
training_slides = training["slide"].astype(str)


###### leave-one-slide-out cross-validation ######

cv_prediction = pd.Series(index=training.index, dtype="object")
cv_probability = pd.Series(index=training.index, dtype="float64")

for held_out_slide in sorted(training_slides.unique()):
    train_rows = training_slides != held_out_slide
    test_rows = training_slides == held_out_slide

    one_model = make_logistic_model()
    one_model.fit(
        training_features.loc[train_rows],
        training_labels.loc[train_rows]
    )

    cv_prediction.loc[test_rows] = one_model.predict(
        training_features.loc[test_rows]
    )
    cv_probability.loc[test_rows] = probability_of_TNC(
        one_model,
        training_features.loc[test_rows]
    )

if cv_prediction.isna().any() or cv_probability.isna().any():
    raise ValueError("Missing LOSO predictions.")

cv_table = training[
    ["tls_id", "slide", "fov", "tls_subtype"]
].copy()
cv_table["prob_TNC_TLS"] = cv_probability
cv_table["predicted_TLS_subtype"] = cv_prediction
cv_table["correct"] = (
    cv_table["tls_subtype"] == cv_table["predicted_TLS_subtype"]
)
cv_table.to_csv(
    output_dir / "serial_TLS_classifier_LOSO_predictions.csv",
    index=False
)

binary_label = (training_labels == positive_class).astype(int)
metrics = pd.DataFrame([{
    "model": "logistic_l2",
    "validation": "leave_one_slide_out",
    "n_training_TLS": training.shape[0],
    "n_training_slides": training_slides.nunique(),
    "accuracy": accuracy_score(training_labels, cv_prediction),
    "balanced_accuracy": balanced_accuracy_score(
        training_labels,
        cv_prediction
    ),
    "roc_auc_TNC_TLS": roc_auc_score(binary_label, cv_probability)
}])
metrics.to_csv(
    output_dir / "serial_TLS_classifier_LOSO_metrics.csv",
    index=False
)


###### refit the logistic model on the complete training cohort ######

logistic_model = make_logistic_model()
logistic_model.fit(training_features, training_labels)

logistic_coefficients = pd.DataFrame({
    "feature": features,
    "coefficient_for_TNC_TLS": (
        logistic_model.named_steps["classifier"].coef_[0]
    )
})
logistic_coefficients.to_csv(
    output_dir / "serial_TLS_logistic_coefficients.csv",
    index=False
)


###### prepare the serial-section TLS composition ######

serial = pd.read_csv(serial_path, dtype={"slide": str, "fov": str})

parsed_fov = serial["fov"].str.extract(r"^(.+)-Section(\d+)$")
if parsed_fov.isna().any().any():
    raise ValueError("Some serial FOV names do not contain a valid section number.")

source_tls_number = parsed_fov[0].str.extract(r"TLS(\d+)$")[0].astype(int)
serial["TLS"] = source_tls_number.map(tls_number_map)
serial["TLS"] = serial["TLS"].apply(
    lambda value: f"TLS{int(value)}" if pd.notna(value) else np.nan
)
serial["section_index"] = parsed_fov[1].astype(int)

serial = serial[serial["TLS"].notna()].copy()
serial["fov"] = (
    serial["TLS"] + "-Section" + serial["section_index"].astype(str)
)

for tls_name in tls_names:
    one_tls = serial[serial["TLS"] == tls_name]
    observed_sections = sorted(one_tls["section_index"].unique().tolist())
    if one_tls.shape[0] != 17 or observed_sections != expected_sections:
        raise ValueError(f"{tls_name} does not contain all 17 sections.")

tls_order = {tls_name: index for index, tls_name in enumerate(tls_names)}
serial["tls_order"] = serial["TLS"].map(tls_order)
serial = serial.sort_values(["tls_order", "section_index"])

# The 17 CODEX sections correspond to original section numbers 1, 4, ..., 49.
serial["section_number"] = 1 + (serial["section_index"] - 1) * 3
serial["z_position_um"] = (serial["section_number"] - 1) * 5

serial_columns = [
    "B", "B-CD21", "B-GC", "PC", "CD4T", "CD8T", "Bile_Duct",
    "DC", "Endothelial", "Fibroblast", "HEV", "Lymphatic",
    "Macrophage", "Neutrophil"
]
for column in serial_columns:
    if column not in serial.columns:
        serial[column] = 0.0

# Serial values are fractions; multiplication by 100 matches the training scale.
serial_features = pd.DataFrame({
    "B": (
        serial["B"] + serial["B-CD21"] + serial["B-GC"] + serial["PC"]
    ) * 100,
    "CD4T": serial["CD4T"] * 100,
    "CD8T": serial["CD8T"] * 100,
    "Cholangiocyte_like": serial["Bile_Duct"] * 100,
    "DC": serial["DC"] * 100,
    "Endothelial": serial["Endothelial"] * 100,
    "Fibroblast": serial["Fibroblast"] * 100,
    "HEV": serial["HEV"] * 100,
    "Lymphatic": serial["Lymphatic"] * 100,
    "Macrophage": serial["Macrophage"] * 100,
    "Neutrophil": serial["Neutrophil"] * 100
})[features].fillna(0).astype(float)

serial["prob_TNC_logistic"] = probability_of_TNC(
    logistic_model,
    serial_features
)
serial["pred_subtype_logistic"] = np.where(
    serial["prob_TNC_logistic"] >= 0.5,
    "TNC-TLS",
    "BF-TLS"
)

serial.to_csv(
    output_dir / "serial_TLS_predictions_TLS.csv",
    index=False
)


###### z-axis stability summaries ######

stability_rows = []

for tls_name in tls_names:
    one_tls = serial[serial["TLS"] == tls_name]
    one_tls = one_tls.sort_values("section_index")
    labels = one_tls["pred_subtype_logistic"].tolist()
    counts = one_tls["pred_subtype_logistic"].value_counts()

    stability_rows.append({
        "TLS": tls_name,
        "n_sections": one_tls.shape[0],
        "n_BF_TLS": int((one_tls["pred_subtype_logistic"] == "BF-TLS").sum()),
        "n_TNC_TLS": int((one_tls["pred_subtype_logistic"] == "TNC-TLS").sum()),
        "majority_subtype": counts.index[0],
        "majority_fraction": counts.iloc[0] / one_tls.shape[0],
        "perfectly_stable": counts.shape[0] == 1,
        "switch_count": sum(
            labels[index] != labels[index - 1]
            for index in range(1, len(labels))
        ),
        "median_probability_TNC": one_tls["prob_TNC_logistic"].median()
    })

stability = pd.DataFrame(stability_rows)
stability.to_csv(
    output_dir / "serial_TLS_z_axis_stability_summary.csv",
    index=False
)

overall_stability = pd.DataFrame([{
    "n_TLS": stability["TLS"].nunique(),
    "n_perfectly_stable": int(stability["perfectly_stable"].sum()),
    "n_majority_fraction_ge_0_8": int(
        (stability["majority_fraction"] >= 0.8).sum()
    ),
    "median_majority_fraction": stability["majority_fraction"].median(),
    "median_switch_count": stability["switch_count"].median()
}])
overall_stability.to_csv(
    output_dir / "serial_TLS_z_axis_stability_overall.csv",
    index=False
)
