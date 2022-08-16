#!/usr/local/bin/python3

# gradient boosted decision tree model 

# BPSM pandas lesson very good
# Machine learning workshop very good
# Random forest BA tutorial also good
# https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html
# https://www.freecodecamp.org/news/how-to-use-the-tree-based-algorithm-for-machine-learning/#:~:text=The%20random%20forest%20algorithm%20works%20by%20completing%20the,for%20every%20predicted%20result.%20...%20More%20items...%20
# above is good and also for preprocessing scaling which I think I need to do

# might need to run python3 -m pip install  sklearn --user
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn import ensemble
from sklearn import metrics
from sklearn.preprocessing import StandardScaler, MinMaxScaler, LabelEncoder

# load the data
data = pd.read_csv('PDB_results_clean.txt', sep="\t")

# separate data into features and outcome (column 7 onwards includes both the vina rmsd and over/under 2A)
features = list(data.columns)[9:]   	
x = data[features]
y = data['SGBN']

print("x is: ", x)
print("y is: ", y)

# standardise the independent features
scaler = StandardScaler()
x_scaled = scaler.fit_transform(x)
print("x_scaled is: ", x_scaled)

# split into test and train
x_train, x_test, y_train, y_test = train_test_split(x_scaled, y, train_size = 0.7, random_state =  42)

# make random forest model 
np.random.seed(321)
#mdl = ensemble.RandomForestClassifier()		# play around with some of the params (max_depth, min_samples_leaf) as default could be too heavy 
mdl = ensemble.GradientBoostingClassifier()
mdl = mdl.fit(x_train, y_train)
print("Model trained")

# use model to make predictions for the test set
y_pred = mdl.predict(x_test)

# Calculate Model Accuracy
print("Accuracy:", metrics.accuracy_score(y_test, y_pred))

# Print the confusion matrix
print(metrics.confusion_matrix(y_test, y_pred))
metrics.ConfusionMatrixDisplay.from_predictions(y_test, y_pred)
plt.show()

# Print the precision and recall, among other metrics
print(metrics.classification_report(y_test, y_pred, digits=3))

# check Important features
feature_importances_df = pd.DataFrame(
    {"feature": list(x.columns), "importance": mdl.feature_importances_}
).sort_values("importance", ascending=False)

# Display top 10 most important features
print(feature_importances_df.iloc[0:10])

# Can then take only the top 100 or so most important features and repeat.
# Maybe also create a bar plot of this
