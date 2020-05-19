import data_handling as dh
import dimred as dm
import train as tr

from os.path import dirname, join as pjoin, exists
from sklearn.preprocessing import StandardScaler
from sklearn import manifold
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA, QuadraticDiscriminantAnalysis as QDA
from sklearn.neighbors import KNeighborsClassifier as KNN
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split, StratifiedKFold, RandomizedSearchCV, learning_curve, ShuffleSplit
from sklearn.metrics import confusion_matrix, accuracy_score, auc, roc_auc_score, roc_curve, log_loss, f1_score
import numpy as np
from scipy import interp
import scipy.io as sio 
from scipy.spatial.distance import chebyshev, correlation, hamming, minkowski
import matplotlib.pyplot as plt
from matplotlib import interactive, is_interactive
import sys
import warnings
from sklearn.exceptions import ConvergenceWarning
import time


# # Number of trees in random forest
# n_estimators = [int(x) for x in np.linspace(start = 20, stop = 500, num = 20)]
# # Number of features to consider at every split
# max_features = ['auto', 'sqrt']
# # Maximum number of levels in tree
# max_depth = [int(x) for x in np.linspace(10, 110, num = 11)]
# max_depth.append(None)
# # Minimum number of samples required to split a node
# min_samples_split = [2, 5, 10]
# # Minimum number of samples required at each leaf node
# min_samples_leaf = [1, 2, 4]
# # Method of selecting samples for training each tree
# bootstrap = [True, False]
# # Create the random grid
# random_grid = {'n_estimators': n_estimators,
#                'max_features': max_features,
#                'max_depth': max_depth,
#                'min_samples_split': min_samples_split,
#                'min_samples_leaf': min_samples_leaf,
#                'bootstrap': bootstrap}
# print(random_grid)

# # Use the random grid to search for best hyperparameters
# # First create the base model to tune
# #rf = RandomForestClassifier(class_weight={0:1.5, 1:1.0})
# rf = RandomForestClassifier()
# # Random search of parameters, using 3 fold cross validation, 
# # search across 100 different combinations, and use all available cores
# rf_random = RandomizedSearchCV(estimator = rf, param_distributions = random_grid, n_iter = 100, cv = 3, verbose=2, random_state=2, n_jobs = -1)

# input_set = 'unique'
# test_indexes = dh.get_test_indexes(input_set)
# print("Test indexes", test_indexes)
# train_indexes = [x for x in   dh.subset_indexes(input_set) if x not in test_indexes]
# print("Train indexes", train_indexes)
# feature_set = 'spect+mlbp'
# print('\n')

# train_data, train_labels, test_data, test_labels = get_scaled_data(train_indexes, test_indexes, 'spect+mlbp', 'PCA', \
# 	20, 'ICA', 20, False)

# # Fit the random search model
# #rf_random.fit(train_data, train_labels)
# #print(rf_random.best_params_)

# rf = RandomForestClassifier(n_estimators= 20, n_jobs=-1, random_state=2)

# rf.fit(train_data, train_labels)
# predictions = rf.predict(test_data)
# scores = rf.predict_proba(test_data)
# classifier_stats = get_classifier_stats(test_labels, predictions)
# print(classifier_stats)

#with harsh_penalty, 3 scale lbp 
# {'n_estimators': 2000, 'min_samples_split': 5, 'min_samples_leaf': 2, 'max_features': 'auto', 'max_depth': 50, 'bootstrap': False}
# Pred [1 0 1 1 1 1 0 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1]
# Truth [1 1 1 1 1 1 0 0 1 1 1 1 1 1 0 0 1 1 1 1 0 0]

#without harsh penalty, 3 scale lbp 
#f = RandomForestClassifier(n_estimators= 45, min_samples_split= 10, min_samples_leaf=24, max_features= 'sqrt', \
#	max_depth= 40, bootstrap= False, n_jobs=-1, random_state=2)


def plot_learning_curve(estimator, title, X, y, ylim=None, cv=None,
                        n_jobs=None, train_sizes=np.linspace(.1, 1.0, 5)):
    plt.figure()
    plt.title(title)
    if ylim is not None:
        plt.ylim(*ylim)
    plt.xlabel("Training examples")
    plt.ylabel("Score")
    train_sizes, train_scores, test_scores = learning_curve(
        estimator, X, y, cv=cv, n_jobs=n_jobs, train_sizes=train_sizes)
    train_scores_mean = np.mean(train_scores, axis=1)
    train_scores_std = np.std(train_scores, axis=1)
    test_scores_mean = np.mean(test_scores, axis=1)
    test_scores_std = np.std(test_scores, axis=1)
    plt.grid()

    plt.fill_between(train_sizes, train_scores_mean - train_scores_std,
                     train_scores_mean + train_scores_std, alpha=0.1,
                     color="r")
    plt.fill_between(train_sizes, test_scores_mean - test_scores_std,
                     test_scores_mean + test_scores_std, alpha=0.1, color="g")
    plt.plot(train_sizes, train_scores_mean, 'o-', color="r",
             label="Training score")
    plt.plot(train_sizes, test_scores_mean, 'o-', color="g",
             label="Cross-validation score")

    plt.legend(loc="best")
    return plt

def get_testing_classifiers():
	return ["Random Forest-100-gini-log2-none-balanced",
			"Random Forest-100-gini-log2-none-harsh_penalty",
			"Random Forest-100-entropy-sqrt-none-balanced",
			"Random Forest-100-entropy-sqrt-none-harsh_penalty",
			"Random Forest-200-entropy-sqrt-none-balanced",
			"Random Forest-200-entropy-sqrt-none-harsh_penalty",
			"Random Forest-20-entropy-sqrt-none-balanced",
			"SVM-rbf-3.0-auto-True-harsh_penalty",
			"KNN-5-correlation-distance", 
			"KNN-3-chebyshev-distance"]

def do_testing(classifier_name, train_indexes, test_indexes, input_set, with_rgb, \
	feature_set, reduction1, reduction2, n_components1, n_components2, 	print_results=False):

		classifier, dimred1, dimred2, scaler1, scaler2 = tr.train_classifier(train_indexes, classifier_name,\
		 feature_set, reduction1, n_components1, reduction2, n_components2, with_rgb)
		pred_labels, scor, test_labels = tr.test_classifier(test_indexes, classifier, feature_set, dimred1,\
		 dimred2, scaler1, scaler2, with_rgb)
		classifier_stats = tr.get_classifier_stats(test_labels, pred_labels) 
		balAc = classifier_stats[bal_accuracy_id()]
		acc = classifier_stats[accuracy_id()]
		auc = classifier_stats[auc_id()]
		specificity = classifier_stats[specificity_id()] 
		sensitivity = classifier_stats[sensitivity_id()]
		if (specificity < 1 and sensitivity < 1):
			doh = specificity * sensitivity / (1 - sensitivity) / (1 - specificity) 
		else:
			doh = 0
		conf = '-'.join([classifier_name, feature_set, reduction1, str(n_components1), reduction2,  str(n_components2)])
		malignancy_prob = scor[:,1].transpose()

		if (print_results):
			print('Configuration:' , conf)
			print('AUC', auc, 'BalAcc', balAc , 'sensitivity', sensitivity, 'specificity', specificity, 'DOH', doh,  'Accuracy', acc)
			print('Scores:', malignancy_prob)
			print('Predictions:', pred_labels)
			print('True labels:', test_labels)
			print('\n')

		return conf, auc, balAc, sensitivity, specificity, doh, acc, pred_labels

def compare_testing_performance(input_set, with_rgb, feature_set, reduction1, reduction2, n_components1, n_components2, classifiers=None):

	print_results = False
	get_learning_curves = True 
	prefix  = "RGB" if with_rgb else "MSI"
	testing_filename = pjoin(out_dir, prefix + "_testing_" + dh.get_scales_in_use() + "scales.csv")

	if classifiers == None:
		classifiers =  get_testing_classifiers()

	if len(classifiers) > 1:
		dh.write_file(testing_filename, ",".join(["Configuration", "AUC", "BalAcc", "Accuracy", "Sensitivity", "Specificity", "DOH", '\n'])) 

	test_indexes = dh.get_test_indexes(input_set)
	train_indexes = [x for x in   dh.subset_indexes(input_set) if x not in test_indexes]

	for classifier_name in classifiers:
		 conf, auc, balAc, sensitivity, specificity, doh, acc, pred_labels = do_testing(classifier_name, train_indexes, test_indexes, \
		 	input_set, with_rgb, feature_set, reduction1, reduction2, n_components1, n_components2, print_results)
		
		if len(classifiers) > 1:
			dh.append_file(testing_filename, ",".join([conf, str(auc), str(balAc) , str(acc), str(sensitivity),  str(specificity), str(doh), '\n']))

		if (get_learning_curves and (acc > 0.8 or doh > 4 ) and any(pred_labels)):

			cv = ShuffleSplit(n_splits=50, test_size=0.2, random_state=0)
			clf = tr.get_classifier(classifier_name)
			train_data, train_labels, test_data, test_labels = tr.get_scaled_data(train_indexes, test_indexes, feature_set, reduction1, \
			n_components1, reduction2, n_components2, with_rgb)
			plot_learning_curve(clf, classifier_name, train_data, train_labels, ylim=(0.5, 1.01), cv=cv, n_jobs=4)
			plt.show()
			#plt.show('hold')
			plt.savefig(pjoin(out_dir, prefix + "testing_scales" + dh.get_scales_in_use(), conf + '.png'), dpi=300, bbox_inches='tight')



compare_testing_performance('unique', False, 'spect', 'ICA', 'ICA',  20, 20, {"Random Forest-100-gini-log2-none-balanced"})
compare_testing_performance('unique', False, 'spect+clbp', 'ICA', 'ICA',  20, 20, {"Random Forest-100-gini-log2-none-balanced"})
compare_testing_performance('unique', False, 'spect+mlbp', 'ICA', 'ICA',  20, 20, {"Random Forest-100-gini-log2-none-balanced"})
compare_testing_performance('unique', False, 'spect+slbp', 'ICA', 'ICA',  20, 20, {"Random Forest-100-gini-log2-none-balanced"})
compare_testing_performance('unique', True, 'spect', 'ICA', 'ICA',  20, 20, {"Random Forest-100-gini-log2-none-balanced"})
compare_testing_performance('unique', True, 'spect+slbp', 'ICA', 'ICA',  20, 20, {"Random Forest-100-gini-log2-none-balanced"})

compare_testing_performance('unique', False, 'spect', 'PCA', 'None',  20, None, {"KNN-3-chebyshev-distance"})
compare_testing_performance('unique', False, 'spect+clbp', 'PCA', 'None',  20, None, {"KNN-3-chebyshev-distance"})
compare_testing_performance('unique', False, 'spect+mlbp', 'PCA', 'None',  20, None, {"KNN-3-chebyshev-distance"})
compare_testing_performance('unique', False, 'spect+slbp', 'PCA', 'None',  20, None, {"KNN-3-chebyshev-distance"})
compare_testing_performance('unique', True, 'spect', 'PCA', 'None',  20, None, {"KNN-3-chebyshev-distance"})
compare_testing_performance('unique', True, 'spect+slbp', 'PCA', 'None',  20, None, {"KNN-3-chebyshev-distance"})

compare_testing_performance('unique', False, 'spect', 'None', 'PCA',  None, 20, {"SVM-rbf-3.0-auto-True-harsh_penalty"})
compare_testing_performance('unique', False, 'spect+clbp', 'None', 'PCA',  None, 20, {"SVM-rbf-3.0-auto-True-harsh_penalty"})
compare_testing_performance('unique', False, 'spect+mlbp', 'None', 'PCA',  None, 20, {"SVM-rbf-3.0-auto-True-harsh_penalty"})
compare_testing_performance('unique', False, 'spect+slbp', 'None', 'PCA',  None, 20, {"SVM-rbf-3.0-auto-True-harsh_penalty"})
compare_testing_performance('unique', True, 'spect', 'None', 'PCA',  None, 20, {"SVM-rbf-3.0-auto-True-harsh_penalty"})
compare_testing_performance('unique', True, 'spect+slbp', 'None', 'PCA',  None, 20, {"SVM-rbf-3.0-auto-True-harsh_penalty"})