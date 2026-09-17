import pandas as pd
import joblib
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
import warnings
import os
warnings.filterwarnings('ignore')

# ===================== 1. Train the final model =====================
print("="*60)
print("Step 1: Train the final model (top 100 features by importance)")
print("="*60)

# 1.1 Load training data
train_path = r""
train_df = pd.read_excel(train_path, sheet_name='args_features_matrix')
print(f"Training data shape: {train_df.shape}")

# 1.2 Prepare data (using full gene names)
X_train = train_df.iloc[:, 2:]  # gene abundance data
y_train = train_df['source']    # target variable
feature_names = X_train.columns.tolist()

# 1.3 Compute feature importance
print("\nComputing feature importance...")
rf_base = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1)
rf_base.fit(X_train, y_train)

# Sort features by importance
feature_importances = rf_base.feature_importances_
importance_df = pd.DataFrame({
    'feature': feature_names,
    'importance': feature_importances
}).sort_values('importance', ascending=False)

# 1.4 Select top 100 features (full names)
top_n = 100
selected_features = importance_df['feature'].head(top_n).tolist()
X_selected = X_train[selected_features]

print(f"\nNumber of selected features: {len(selected_features)}")
print(f"Top 10 features:")
for i, (feature, imp) in enumerate(zip(selected_features[:10], importance_df['importance'].head(10).values)):
    print(f"  {i+1:2d}. {feature:30s} importance: {imp:.6f}")

# 1.5 Retrain model with selected features
X_train_final, X_test_final, y_train_final, y_test_final = train_test_split(
    X_selected, y_train, test_size=0.3, random_state=42, stratify=y_train
)

final_model = RandomForestClassifier(
    n_estimators=100,
    random_state=42,
    n_jobs=-1
)
final_model.fit(X_train_final, y_train_final)

# 1.6 Evaluate model
train_acc = final_model.score(X_train_final, y_train_final)
test_acc = final_model.score(X_test_final, y_test_final)

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
cv_scores = cross_val_score(final_model, X_selected, y_train, cv=cv, n_jobs=-1)

print(f"\nModel performance:")
print(f"Training accuracy: {train_acc:.4f}")
print(f"Test accuracy: {test_acc:.4f}")
print(f"5-fold CV mean: {cv_scores.mean():.4f} (+/- {cv_scores.std():.4f})")
print(f"Overfitting degree: {train_acc - test_acc:.4f}")

# 1.7 Save model and features
output_dir = r""
os.makedirs(output_dir, exist_ok=True)

model_path = os.path.join(output_dir, "final_model_100_features.pkl")
features_path = os.path.join(output_dir, "selected_features_100.pkl")
importance_path = os.path.join(output_dir, "feature_importance_ranking.xlsx")

joblib.dump(final_model, model_path)
joblib.dump(selected_features, features_path)
importance_df.to_excel(importance_path, index=False)

print(f"\nModel and features saved to: {output_dir}")
print(f"Model file: {model_path}")
print(f"Feature file: {features_path}")
print(f"Importance ranking: {importance_path}")

# ===================== 2. Load and preprocess prediction data =====================
print("\n" + "="*60)
print("Step 2: Load and preprocess prediction data")
print("="*60)

# 2.1 Load prediction data
predict_path = r""
predict_df = pd.read_excel(predict_path)
print(f"Prediction data shape: {predict_df.shape}")
print(f"Prediction data columns: {predict_df.columns.tolist()}")

# 2.2 Pivot data into sample x gene format (using full sseqid names)
print("\nPivoting data into sample x gene format...")
pivot_df = predict_df.pivot_table(
    index='sample_name',
    columns='sseqid',  # use full sseqid names directly
    values='abundance(copies/cell)',
    aggfunc='mean'
).fillna(0)

print(f"Pivoted data shape: {pivot_df.shape}")
print(f"Number of samples: {pivot_df.shape[0]}")
print(f"Number of genes: {pivot_df.shape[1]}")

# 2.3 Ensure prediction data contains all required features (full names)
print("\nChecking feature matching...")
missing_features = [f for f in selected_features if f not in pivot_df.columns]
print(f"Number of missing features: {len(missing_features)}")
if missing_features:
    print(f"First 10 missing features: {missing_features[:10]}")

# Add missing feature columns filled with 0
for feature in selected_features:
    if feature not in pivot_df.columns:
        pivot_df[feature] = 0

# Ensure feature order matches training
X_predict = pivot_df[selected_features]
print(f"Final prediction data shape: {X_predict.shape}")

# ===================== 3. Make predictions =====================
print("\n" + "="*60)
print("Step 3: Make predictions")
print("="*60)

# 3.1 Predict
predictions = final_model.predict(X_predict)
prediction_probs = final_model.predict_proba(X_predict)

# 3.2 Create results DataFrame
results_df = pd.DataFrame({
    'sample_name': pivot_df.index,
    'predicted_source': predictions,
    'prediction_confidence': prediction_probs.max(axis=1)
})

# Add probability for each class
class_names = final_model.classes_
for i, class_name in enumerate(class_names):
    results_df[f'prob_{class_name}'] = prediction_probs[:, i]

# 3.3 Save prediction results
results_path = os.path.join(output_dir, "prediction_results.xlsx")
results_df.to_excel(results_path, index=False)
print(f"Prediction results saved to: {results_path}")

# 3.4 Display prediction results
print(f"\nPrediction results summary:")
print(f"Number of predicted samples: {len(results_df)}")
print(f"Predicted class distribution:")
print(results_df['predicted_source'].value_counts())

print(f"\nPrediction results sample (first 10):")
print(results_df.head(10))

# ===================== 4. Generate feature matching report =====================
print("\n" + "="*60)
print("Step 4: Generate feature matching report")
print("="*60)

# Count matched features for each sample
matched_counts = []
for idx in pivot_df.index:
    sample_row = pivot_df.loc[idx, selected_features]
    matched = (sample_row > 0).sum()
    matched_counts.append(matched)

results_df['matched_features_count'] = matched_counts
results_df['matched_features_percentage'] = results_df['matched_features_count'] / top_n

# Save detailed report
detailed_path = os.path.join(output_dir, "detailed_prediction_report.xlsx")
results_df.to_excel(detailed_path, index=False)
print(f"Detailed prediction report: {detailed_path}")

# Print summary
print(f"\n" + "="*60)
print("Task Complete Summary")
print("="*60)
print(f"OK Model training done: using {top_n} important features")
print(f"OK Model performance: test accuracy {test_acc:.4f}, CV mean {cv_scores.mean():.4f}")
print(f"OK Predictions done: {len(results_df)} samples")
print(f"OK Average feature match rate: {results_df['matched_features_percentage'].mean():.4f}")
print(f"OK Results saved: {output_dir}")

# Show high-confidence predictions
high_confidence = results_df[results_df['prediction_confidence'] > 0.8]
if len(high_confidence) > 0:
    print(f"\nHigh-confidence predictions (>0.8): {len(high_confidence)} samples")
    print(high_confidence[['sample_name', 'predicted_source', 'prediction_confidence']].head())
else:
    print("\nWarning: no high-confidence (>0.8) predictions found")
