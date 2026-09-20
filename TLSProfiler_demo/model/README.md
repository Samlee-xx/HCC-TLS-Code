# TLSProfiler classifier

This directory contains the final logistic-regression coefficients used for TLS subtype prediction:

- `tlsprofiler_log1p_BTN.json`: unpenalized logistic-regression coefficients exported for Python inference.

The model uses the TLS-level percentages of B cells, T cells, and neutrophils. Each percentage is transformed with `log1p`. TNC-TLS is the positive class, and the probability cutoff stored in the JSON file is used for subtype assignment.

The JSON file also records the minimum polygon size (`100` cells) and the default output threshold (`300` cells). These two values describe the inference workflow and are not fitted model coefficients. The final output threshold can be changed at inference time.

This release contains the fitted inference model but does not include the training data or model-fitting scripts.
