import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt
import seaborn as sns

# Set style
plt.style.use('default')
sns.set_palette("husl")
plt.rcParams['font.family'] = 'Arial'

# 1. Load data
df = pd.read_excel(r"",
                   sheet_name='')

# 2. Prepare data
X = df.iloc[:, 2:]  # Gene abundance data
y = df['source']    # Target variable
class_names = y.unique()

# 3. Split data
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.3, random_state=42, stratify=y
)

# 4. Define and train model
rf_model = RandomForestClassifier(
    n_estimators=100,
    random_state=42,
    n_jobs=-1
)

rf_model.fit(X_train, y_train)
y_pred = rf_model.predict(X_test)

# 5. 5-fold cross-validation
cv_scores = cross_val_score(rf_model, X, y, cv=5, n_jobs=-1)

# ===================== Visualization 1: Confusion Matrix =====================
plt.figure(figsize=(20, 8))

# Subplot 1: Confusion Matrix
plt.subplot(1, 2, 1)
plt.text(-0.15, 1.08, '(a)', transform=plt.gca().transAxes,
         fontsize=30, fontweight='bold', va='top')

cm = confusion_matrix(y_test, y_pred)
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues',
            xticklabels=class_names, yticklabels=class_names,
            cbar_kws={'label': 'Count'},
            annot_kws={'size': 20,})
plt.xticks(fontsize=18)   # 根据需要调整字号

# 设置y轴刻度标签字体大小
plt.yticks(fontsize=18)

plt.xlabel('Predicted Label', fontsize=18)
plt.ylabel('True Label', fontsize=18)
plt.title('Confusion Matrix', fontsize=20)
plt.tight_layout()

# ===================== Visualization 2: 5-Fold Cross Validation =====================
plt.subplot(1, 2, 2)
plt.text(-0.12, 1.06, '(b)', transform=plt.gca().transAxes,
         fontsize=30, fontweight='bold', va='top')
# Create bar plot for CV scores
fold_numbers = np.arange(1, len(cv_scores) + 1)
bars = plt.bar(fold_numbers, cv_scores, width=0.4, color='#5172a5',
               linewidth=1, alpha=0.8)

# Add value labels on top of bars
for bar, score in zip(bars, cv_scores):
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2., height + 0.005,
            f'{score:.3f}', ha='center', va='bottom', fontsize=18, fontweight='bold')

# Add mean line
mean_score = cv_scores.mean()
plt.axhline(y=mean_score, color='#ffaf94', linestyle='--', linewidth=2,
           label=f'Mean: {mean_score:.3f}')

# Add standard deviation text
std_score = cv_scores.std()
plt.text(0.5, 0.05, f'Std: {std_score:.3f}',
         transform=plt.gca().transAxes, fontsize=18,
         bbox=dict(boxstyle='round', facecolor='#8ea8fb', alpha=0.5))

plt.xlabel('Fold Number', fontsize=18)
plt.ylabel('Accuracy', fontsize=18)
plt.title('5-Fold Cross-Validation Scores', fontsize=20)
plt.xticks(fold_numbers)
plt.ylim(0, 1.1)
plt.legend(loc='lower right',fontsize=18)
plt.grid(True, alpha=0.3, axis='y')
plt.tick_params(axis='both', labelsize=18)

plt.tight_layout()



plt.show()