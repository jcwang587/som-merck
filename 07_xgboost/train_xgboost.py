import pandas as pd
import xgboost as xgb
import matplotlib.pyplot as plt
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RepeatedStratifiedKFold
from imblearn.pipeline import Pipeline
from imblearn.over_sampling import SMOTE
import shap
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import seaborn as sns
from sklearn.model_selection import train_test_split

# Load the dataset
dataset = pd.read_csv("../data/dataset/dataset_merck_bde.csv")
dataset = dataset.dropna()

print(dataset.shape)
# Define the features and target
features = [
    "bde",
    "relative_ir",
    "sasa_hydrogen_maestro",
    "heavy_atoms",
]
target = "som"

X = dataset[features]
y = dataset[target]

# Split the dataset into training and test sets
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.1, random_state=42
)

# Define the pipeline
pipeline = Pipeline(
    [
        ("smote", SMOTE(random_state=42)),
        (
            "classifier",
            xgb.XGBClassifier(objective="binary:logistic", eval_metric="auc"),
        ),
    ]
)

# Define the parameter grid for hyperparameter tuning, 40, 400, 0.05
param_grid = {
    "classifier__n_estimators": [50, 100, 150, 200],
    "classifier__max_depth": [3, 5, 7, 9],
    "classifier__learning_rate": [0.01, 0.05, 0.1],
}

# Define the evaluation procedure
cv = RepeatedStratifiedKFold(n_splits=10, n_repeats=3, random_state=1)

# Initialize GridSearchCV
grid_search = GridSearchCV(
    estimator=pipeline,
    param_grid=param_grid,
    scoring="roc_auc",
    cv=cv,
    n_jobs=-1,
    verbose=2,
)

# Fit the grid search to the training data
grid_search.fit(X_train, y_train)

# Get the best model from grid search
best_model = grid_search.best_estimator_
print("Best parameters found:", grid_search.best_params_)

# Save the parameters to a file
with open("best_params_test.txt", "w") as f:
    f.write(str(grid_search.best_params_))

# Get the final classifier
final_classifier = best_model.named_steps["classifier"]

# Plot the confusion matrix on the test set
cm_test = confusion_matrix(y_test, best_model.predict(X_test))
plt.rcParams["font.size"] = "16"
disp_test = ConfusionMatrixDisplay(
    confusion_matrix=cm_test,
    display_labels=best_model.named_steps["classifier"].classes_,
)
disp_test.plot(cmap=plt.cm.Blues)
plt.tight_layout()
plt.savefig("confusion_matrix_test.png")
plt.close()

# Calculate the SHAP values using the final classifier
explainer = shap.Explainer(final_classifier, X_train)
shap_values = explainer(X_test)
shap.summary_plot(shap_values, X_test)
plt.savefig("shap_summary_test.png")
plt.close()