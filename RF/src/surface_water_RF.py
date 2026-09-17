"""
ARG source prediction script.

Trains a Random Forest classifier on labeled ARG abundance profiles,
selects the top 100 important ARG features, and predicts the source
of ARGs in new environmental samples.

Input:
    - summary.xlsx (sheet: args_features_matrix): training data (samples x ARG abundances + source label)
    - profile_summary.xlsx: ARG abundance profiles of query samples

Output:
    - prediction_results.xlsx: predicted source and probabilities for each sample
"""

import os
import warnings
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split

warnings.filterwarnings("ignore")

# ===================== 1. Train model =====================
# Load labeled training data
train_path = r""
train_df = pd.read_excel(train_path, sheet_name="args_features_matrix")

X = train_df.iloc[:, 2:]           # ARG abundance matrix (samples x genes)
y = train_df["source"]             # Source label for each sample
feature_names = X.columns.tolist()

# Fit baseline Random Forest to estimate feature importance
rf_base = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1)
rf_base.fit(X, y)

# Rank features by importance and select top 100
importance_df = pd.DataFrame({
    "feature": feature_names,
    "importance": rf_base.feature_importances_
}).sort_values("importance", ascending=False)

top_n = 100
selected_features = importance_df["feature"].head(top_n).tolist()
X_selected = X[selected_features]

# Split into train/test sets (stratified to preserve class balance)
X_train, X_test, y_train, y_test = train_test_split(
    X_selected, y, test_size=0.3, random_state=42, stratify=y
)

# Train final model on selected features
final_model = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1)
final_model.fit(X_train, y_train)

# ===================== 2. Preprocess query data =====================
# Load ARG abundance profiles of query samples (long format)
predict_path = ""
predict_df = pd.read_excel(predict_path)

# Pivot to wide format: rows = samples, columns = ARG genes (sseqid)
pivot_df = predict_df.pivot_table(
    index="sample_name",
    columns="sseqid",
    values="abundance(copies/cell)",
    aggfunc="mean"
).fillna(0)

# Align query data to selected features (missing genes -> 0, preserve order)
for feature in selected_features:
    if feature not in pivot_df.columns:
        pivot_df[feature] = 0
X_predict = pivot_df[selected_features]

# ===================== 3. Predict =====================
# Predict source class and class probabilities
predictions = final_model.predict(X_predict)
prediction_probs = final_model.predict_proba(X_predict)

# Build results table
results_df = pd.DataFrame({
    "sample_name": pivot_df.index,
    "predicted_source": predictions,
    "prediction_confidence": prediction_probs.max(axis=1)
})

# Add per-class probability columns
for i, class_name in enumerate(final_model.classes_):
    results_df[f"prob_{class_name}"] = prediction_probs[:, i]

# ===================== 4. Save output =====================
output_dir = ""
os.makedirs(output_dir, exist_ok=True)
results_path = os.path.join(output_dir, "prediction_results.xlsx")
results_df.to_excel(results_path, index=False)