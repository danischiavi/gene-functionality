import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, roc_auc_score, roc_curve, confusion_matrix, matthews_corrcoef
from sklearn.impute import KNNImputer
import joblib

#### PROTEIN-CODING FEATURES MATRIX ####
# Load the protein-coding datasets
df1 = pd.read_csv('./results/prediction/Model/functional-protein-exon2-dataset-features.csv')
df2 = pd.read_csv('./results/prediction/Model/functional-protein-exon3-dataset-features.csv')
df3 = pd.read_csv('./results/prediction/Model/protein-exon2-negative-control-dataset-features.csv')
df4 = pd.read_csv('./results/prediction/Model/protein-exon3-negative-control-dataset-features.csv')

# Drop IDs columns
df1 = df1.drop(columns=['ID', 'Chromosome', 'Start', 'End', 'Sequence', 'GeneID'])
df2 = df2.drop(columns=['ID', 'Chromosome', 'Start', 'End', 'Sequence', 'GeneID'])
df3 = df3.drop(columns=['ID', 'Chromosome', 'Start', 'End', 'Sequence', 'DistanceGene'])
df4 = df4.drop(columns=['ID', 'Chromosome', 'Start', 'End', 'Sequence', 'DistanceGene'])

# Merge all df 
df_protein = pd.concat([df1, df3, df2, df4], ignore_index=True)

# Handle missing Value
X = df_protein.drop(columns=['Functional']) # exclude categorical 
y = df_protein['Functional']
X.isna().sum()

imputer = KNNImputer(n_neighbors=5)
imputed_data = imputer.fit_transform(X)
X = pd.DataFrame(imputed_data, columns=X.columns)
X.isna().sum()

df_protein = pd.concat([y, X,], axis=1)
df_protein.to_csv('./results/prediction/Model/protein-coding-features-matrix.csv', index=False)



#### RandomForest MODEL ####

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# PROTEIN-CODING-DATASET
rf_protein = RandomForestClassifier(n_estimators = 250, max_features = 3, max_depth = 10, random_state=42)
rf_protein.fit(X_train, y_train)
y_pred = rf_protein.predict(X_test)

joblib.dump(rf_protein, './results/prediction/Model/RF_protein_model.pkl')

# Metrics
mcc = matthews_corrcoef(y_test, y_pred)
auc_roc = roc_auc_score(y_test, rf_protein.predict_proba(X_test)[:, 1])
tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
sensitivity = tp / (tp + fn)
specificity = tn / (tn + fp)

print(f"Sensitivity: {sensitivity}\nSpecificity: {specificity}\nMCC: {mcc}\nAUC-ROC: {auc_roc}")
print(confusion_matrix(y_test, y_pred))

# Feature's importance: Mean Decrease in Impurity  
importances = rf_protein.feature_importances_
Importance_MDI_df = pd.DataFrame({
	'Feature': X.columns,
	'Importance_MDI': importances
}).sort_values(by='Importance_MDI', ascending=False)


#### PERFORMANCE ANALYSIS #### 
# Bootstrapping
# Initialize variables
n_bootstraps = 100
metrics = {
	'sensitivity': [],
	'specificity': [],
	'mcc': [],
	'auc_roc': []
}

feature_importances = {feature: [] for feature in X.columns}
feature_names = X_test.columns

# Bootstrap loop
for i in range(n_bootstraps):
	# Split data into new training and test sets
	X_train_boot, X_test_boot, y_train_boot, y_test_boot = train_test_split(X, y, test_size=0.2, random_state=i)
		
	# Train model
	rf_protein = RandomForestClassifier(n_estimators=250, max_features=3, max_depth=10, random_state=42)
	rf_protein.fit(X_train_boot, y_train_boot)
	
	# Predict on new test data
	y_pred = rf_protein.predict(X_test_boot)
	
	# Metrics
	mcc = matthews_corrcoef(y_test_boot, y_pred)
	auc_roc = roc_auc_score(y_test_boot, rf_protein.predict_proba(X_test_boot)[:, 1])
	tn, fp, fn, tp = confusion_matrix(y_test_boot, y_pred).ravel()
	sensitivity = tp / (tp + fn)
	specificity = tn / (tn + fp)
	
	# Store metrics
	metrics['sensitivity'].append(sensitivity)
	metrics['specificity'].append(specificity)
	metrics['mcc'].append(mcc)
	metrics['auc_roc'].append(auc_roc)
	importances = rf_protein.feature_importances_
	for feature, importance in zip(X.columns, importances):
		feature_importances[feature].append(importance)

# Convert metrics to DataFrame for analysis
metrics_df = pd.DataFrame(metrics)

summary = {
    'Metric': ['Sens', 'Spec', 'MCC', 'AUC-ROC'],
    'Mean': [
        metrics_df['sensitivity'].mean(),
        metrics_df['specificity'].mean(),
        metrics_df['mcc'].mean(),
        metrics_df['auc_roc'].mean()
    ],
    'Std': [
        metrics_df['sensitivity'].std(),
        metrics_df['specificity'].std(),
        metrics_df['mcc'].std(),
        metrics_df['auc_roc'].std()
    ]
}

summary_df = pd.DataFrame(summary) 

# Convert feature importances to DataFrame for analysis
importance_df = pd.DataFrame(feature_importances)
mean_importances = importance_df.mean()
std_importances = importance_df.std()

# Create summary DataFrame for feature importances
importance_summary = pd.DataFrame({
    'Feature' : feature_names,
	'MDI_mean': mean_importances,
    'MDI_std': std_importances
}).sort_values(by='MDI_mean', ascending=False)

# Save to TSV 
file_path = './results/prediction/Model/models_performance.tsv'
with open(file_path, 'a') as f:
    f.write("\n### PROTEIN-CODING SUMMARY ####\n")
    summary_df.round(3).to_csv(f, sep='\t', index=False)
    f.write("\nFeature Importance\n")
    importance_summary.round(3).to_csv(f, sep='\t', index=False)

print(f"Summary has been appended to '{file_path}'.")

