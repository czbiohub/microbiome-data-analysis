#!/usr/bin/env python3
import pandas as pd
import pickle
import sys
from sklearn.model_selection import train_test_split
# import the class
from sklearn.linear_model import LogisticRegression
from sklearn import metrics
# import required modules
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
%matplotlib inline

# bam_summary_file=sys.argv[1]
bam_summary_file = "/Users/sunit.jain/Research/SyntheticCommunities/StrainAbundance/Bacteroides-coprophilus-DSM-18228/Bacteroides-coprophilus-DSM-18228.processed.sortedByCoord.csv"
# load dataset
bam_summ = pd.read_csv(bam_summary_file, header=0)

feature_cols = ["align_len","query_len","aln_cov","quality","perc_id","aln_score","mate_score","mismatches","gap_open","gap_ext","is_dup","is_primary","is_supp"]
X = bam_summ[feature_cols]
y = bam_summ.is_truth

X_train,X_test,y_train,y_test=train_test_split(X,y,test_size=0.25,random_state=0)

# instantiate the model (using the default parameters)
logreg = LogisticRegression(solver="lbfgs")

# fit the model with data
logreg.fit(X_train,y_train)

#
y_pred=logreg.predict(X_test)

filename="bam_logistic_regression.pkl"
pickle.dump(logreg, open(filename, 'wb'))

#
cnf_matrix = metrics.confusion_matrix(y_test, y_pred)
# cnf_matrix

class_names=[0,1] # name  of classes
fig, ax = plt.subplots()
tick_marks = np.arange(len(class_names))
plt.xticks(tick_marks, class_names)
plt.yticks(tick_marks, class_names)
# create heatmap
sns.heatmap(pd.DataFrame(cnf_matrix), annot=True, cmap="YlGnBu" ,fmt='g')
ax.xaxis.set_label_position("top")
plt.tight_layout()
plt.title('Confusion matrix', y=1.1)
plt.ylabel('Actual label')
plt.xlabel('Predicted label')

print("Accuracy:",metrics.accuracy_score(y_test, y_pred))
print("Precision:",metrics.precision_score(y_test, y_pred))
print("Recall:",metrics.recall_score(y_test, y_pred))

y_pred_proba = logreg.predict_proba(X_test)[::,1]
fpr, tpr, _ = metrics.roc_curve(y_test,  y_pred_proba)
auc = metrics.roc_auc_score(y_test, y_pred_proba)
plt.plot(fpr,tpr,label="data 1, auc="+str(auc))
plt.legend(loc=4)
plt.show()