import data_handling as dh
import dimred as dm

from os.path import dirname, join as pjoin, exists
from sklearn.preprocessing import StandardScaler
from sklearn import manifold
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA, QuadraticDiscriminantAnalysis as QDA
from sklearn.neighbors import KNeighborsClassifier as KNN
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split, StratifiedKFold, RandomizedSearchCV
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

#warnings.filterwarnings('error')

interactive(True)

out_dir = dh.get_out_dir()

#sys.stdout = open(classification_log, 'w+')

####################classification#######################
def get_classifier_stats(y_test, y_pred):
	tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()

	sensitivity =  tp / (tp + fn) if (tp + fp) > 0 else float('nan')
	specificity = tn / (tn + fp) if (tn + fp) > 0 else float('nan')
	precision = tp / (tp + fp) if (tp + fp) > 0 else float('nan')
	false_positive_rate = fp / (fp + tn) if (fp + tn) > 0 else float('nan')
	false_negative_rate = fn / (fn + tp) if (fn + tp) > 0 else float('nan')
	try:
		f1 = f1_score(y_test, y_pred)
	except:
		f1 = float('nan')
	accuracy = accuracy_score(y_test, y_pred)  #(tp + tn) / (tp + tn + fp + fn)
	bal_accuracy = (sensitivity + specificity) / 2
	auc = roc_auc_score(y_test, y_pred)
	classifier_stats = [accuracy, auc, specificity, sensitivity, false_positive_rate, false_negative_rate, bal_accuracy, f1]
	return classifier_stats

def accuracy_id():
	return 0
def auc_id():
	return 1
def specificity_id():
	return 2
def sensitivity_id():
	return 3
def false_positive_rate_id():
	return 4
def false_negative_rate_id():
	return 5
def bal_accuracy_id():
	return 6
def f1_id():
	return 7


def plot_cv_roc_curve(fpr_f, tpr_f, tprs, aucs):

	for i in range(len(aucs)):
		plt.plot(fpr, tpr, lw=1, alpha=0.3, label='ROC fold %d (AUC = %0.2f)' % (i+1, aucs[i]))

	plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r', label='Chance', alpha=.8)

	mean_tpr = np.mean(tprs, axis=0)
	mean_tpr[-1] = 1.0
	mean_auc = auc(mean_fpr, mean_tpr)
	std_auc = np.std(aucs)
	plt.plot(mean_fpr, mean_tpr, color='b', label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc), lw=2, alpha=.8)

	std_tpr = np.std(tprs, axis=0)
	tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
	tprs_lower = np.maximum(mean_tpr - std_tpr, 0)

	plt.xlim([-0.05, 1.05])
	plt.ylim([-0.05, 1.05])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	plt.title('Receiver operating characteristic curve for CrossVal')
	plt.legend(loc="lower right")
	plt.show()

def plot_accuracy_and_auc(performance_stats, performance_names, plot_title='Classification Performance Comparison'):
	fig, ax1 = plt.subplots()
	x = np.arange(len(performance_stats))
	accuracy = performance_stats[:,accuracy_id()]
	auc_score = performance_stats[:,auc_id()]

	plt.title(plot_title)
	plt.xticks(x, performance_names, rotation='vertical')
	ax1.plot(x, accuracy, 'b-')
	ax1.set_xlabel('classifier')
	# Make the y-axis label, ticks and tick labels match the line color.
	ax1.set_ylabel('Accuracy', color='b')
	ax1.tick_params('y', colors='b')

	ax2 = ax1.twinx()
	ax2.plot(x, auc_score, 'r.')
	ax2.set_ylabel('AUC', color='r')
	ax2.tick_params('y', colors='r')

	fig.tight_layout()
	plt.show()

	dh.create_directory(pjoin(out_dir, 'classification'))
	img_title = plot_title.replace(' ', '_').lower()
	plt.savefig(pjoin(out_dir, 'classification', img_title + '.png'), bbox_inches='tight')

def get_classifier(classifier_name):

	if "/" in classifier_name:
		classifier_name = classifier_name[:classifier_name.find("/")]

	if 'SVM' in classifier_name:
		a, b, c, d, e, f = classifier_name.split('-') 
		c = float(c)
		d = d if (d == 'auto' or d == 'scale') else float(d)
		if f == 'harsh_penalty':
			f = {0:1.5, 1:1.0}
		else:
			f = 'balanced'
		clf = SVC(kernel=b, C=c, gamma=d, shrinking=(e==True), probability=True, random_state=2, class_weight=f)

	elif 'KNN' in classifier_name:
		a, b, c, d =  classifier_name.split('-')
		clf = KNN(n_neighbors=int(b), weights=d, algorithm = 'auto', metric=c, n_jobs=-1)

	elif 'QDA' in classifier_name:
		clf = QDA()

	elif 'LDA' in classifier_name:
		clf = LDA()

	elif 'Decision Tree' in classifier_name:
		a, b=  classifier_name.split('-')
		clf = DecisionTreeClassifier(criterion=b)

	elif 'Random Forest' in classifier_name:
		a, b, c, d, e, f = classifier_name.split('-') 
		d = None if (d == 'none') else d
		e = None if (e == 'none') else int(e)
		if f == 'harsh_penalty':
			f = {0:1.5, 1:1.0}
		else:
			f = 'balanced'
		clf = RandomForestClassifier(n_estimators=int(b), criterion=c, max_features=d, max_depth=e, class_weight=f, n_jobs=-1, random_state=2)

	elif 'Naive Bayes' in classifier_name:
		clf = GaussianNB()

	else:
		print(classifier_name, 'Not yet implemented.')
		return

	return clf

def train_classifier(train_indexes, classifier_name, feature_set, reduction1='None', n_components1=10, reduction2='None', n_components2=10, rgb=False):
	labels = dh.get_labels()
	reference = dh.get_measured_spectra()
	if rgb:
		data = dh.get_reconstructed_spectra_rgb()
	else:
		data = dh.get_reconstructed_spectra()


	clf = get_classifier(classifier_name)
	dimred1 = dm.get_dimred(reduction1, n_components1)
	dimred2 = dm.get_dimred(reduction2, n_components2)

	scaler1 = StandardScaler()
	reference_data = scaler1.fit_transform(reference[train_indexes, :])
	train_spect, train_labels, train_lbp = dh.get_scaled_subset_with_index(train_indexes, data, labels, scaler1, feature_set, rgb)

	scaler2 = StandardScaler()
	if "spect" in feature_set and reduction1 is not 'None':
		dimred1 = dimred1.fit(reference_data, train_labels)
		train_spect = dimred1.transform(train_spect)
		#print("PCA-explained variance for specrum", dimred.explained_variance_ratio_ * 100)
	else:
		dimred1=None

	if "spect" not in feature_set:
		train_spect=None

	if "lbp" in feature_set and reduction2 is not 'None':
		dimred2 = dimred2.fit(train_lbp, train_labels)
		train_lbp = dimred2.transform(train_lbp)
			#print("PCA-explained variance for texture", dimred2.explained_variance_ratio_ * 100)
	else:
		dimred2=None

	if "lbp" not in feature_set:
		train_lbp=None

	train_data = dh.concat_features(train_spect, train_lbp)
	train_data = scaler2.fit_transform(train_data)

	#print(clf) 
	clf.fit(train_data, train_labels)

	return clf, dimred1, dimred2, scaler1, scaler2


def test_classifier(test_indexes, classifier, feature_set, dimred1, dimred2, scaler1, scaler2,rgb=False):
	if rgb: 
		data = dh.get_reconstructed_spectra_rgb()
	else:
		data = dh.get_reconstructed_spectra()
	labels = dh.get_labels()

	test_spect, test_labels, test_lbp = dh.get_scaled_subset_with_index(test_indexes, data, labels, scaler1, feature_set, rgb)
	if "spect" in feature_set and dimred1 is not None:
		test_spect = dimred1.transform(test_spect)
	if "spect" not in feature_set:
		test_spect = None 

	if "lbp" in feature_set and dimred2 is not None:
		test_lbp = dimred2.transform(test_lbp)
	if "lbp" not in feature_set:
		test_lbp = None

	test_data = dh.concat_features(test_spect, test_lbp)
	test_data = scaler2.transform(test_data)
	predictions = classifier.predict(test_data)
	scores = classifier.predict_proba(test_data)

	classifier_params = classifier.get_params()
	if 'n_neighbors' in classifier_params:
		dist, neighbor_ids  = classifier.kneighbors(test_data, classifier_params['n_neighbors'])
		dist = np.mean(dist, axis=1)
		scores[:,0] = dist.transpose()
		scores[:,1] = (1 - dist).transpose()

	return predictions, scores, test_labels

def get_scaled_data(train_indexes, test_indexes, feature_set, reduction1, n_components1, reduction2, n_components2, with_rgb=False):	
	labels = dh.get_labels()
	reference = dh.get_measured_spectra()
	if with_rgb:
		data = dh.get_reconstructed_spectra_rgb()
	else:
		data = dh.get_reconstructed_spectra()


	dimred1 = dm.get_dimred(reduction1, n_components1)
	dimred2 = dm.get_dimred(reduction2, n_components2)

	scaler1 = StandardScaler()
	reference_data = scaler1.fit_transform(reference[train_indexes, :])
	train_spect, train_labels, train_lbp = dh.get_scaled_subset_with_index(train_indexes, data, labels, scaler1, feature_set, with_rgb)

	scaler2 = StandardScaler()
	if "spect" in feature_set and reduction1 is not 'None':
		dimred1 = dimred1.fit(reference_data, train_labels)
		train_spect = dimred1.transform(train_spect)
	else:
		dimred1=None

	if "spect" not in feature_set:
		train_spect=None

	if "lbp" in feature_set and reduction2 is not 'None':
		dimred2 = dimred2.fit(train_lbp, train_labels)
		train_lbp = dimred2.transform(train_lbp)
	else:
		dimred2=None

	if "lbp" not in feature_set:
		train_lbp=None

	train_data = dh.concat_features(train_spect, train_lbp)
	train_data = scaler2.fit_transform(train_data)

	test_spect, test_labels, test_lbp = dh.get_scaled_subset_with_index(test_indexes, data, labels, scaler1, feature_set, with_rgb)
	if "spect" in feature_set and dimred1 is not None:
		test_spect = dimred1.transform(test_spect)
	if "spect" not in feature_set:
		test_spect = None 

	if "lbp" in feature_set and dimred2 is not None:
		test_lbp = dimred2.transform(test_lbp)
	if "lbp" not in feature_set:
		test_lbp = None

	test_data = dh.concat_features(test_spect, test_lbp)
	test_data = scaler2.transform(test_data)

	return train_data, train_labels, test_data, test_labels

def get_fold_data(subset_name, stratified, feature_set, reduction1, n_components1, reduction2, n_components2, with_rgb=False, ignore_test=False):
	if stratified: 
		fold_indexes, folds = dh.get_fold_indexes_stratified(subset_name, 5, ignore_test)
	else:
		fold_indexes, folds = dh.get_fold_indexes(subset_name, 5, ignore_test)

	train_data_per_fold = []
	train_labels_per_fold = []
	test_data_per_fold = []
	test_labels_per_fold = []

	for f in range(folds): 
		test_indexes = fold_indexes[f]
		train_indexes = [y for x in fold_indexes if x != test_indexes for y in x ]

		train_data, train_labels, test_data, test_labels = get_scaled_data(train_indexes, test_indexes, feature_set, reduction1, \
			n_components1, reduction2, n_components2, with_rgb)

		train_data_per_fold.append(train_data) 
		train_labels_per_fold.append(train_labels) 
		test_data_per_fold.append(test_data) 
		test_labels_per_fold.append(test_labels) 

	return folds, train_data_per_fold, train_labels_per_fold, test_data_per_fold, test_labels_per_fold

def cross_validate(classifier_name, subset_name, stratified, feature_set, reduction1, n_components1, reduction2, n_components2, with_rgb=False, ignore_test=False):

	tprs = []
	fpr_f = []
	tpr_f = []
	aucs = []
	mean_fpr = np.linspace(0, 1, 100)
	stats = []

	if stratified: 
		fold_indexes, folds = dh.get_fold_indexes_stratified(subset_name, 5, ignore_test)
	else:
		fold_indexes, folds = dh.get_fold_indexes(subset_name, 5, ignore_test)

	for f in range(folds): 
		test_indexes = fold_indexes[f]
		train_indexes = [y for x in fold_indexes if x != test_indexes for y in x ]

		classifier, dimred1, dimred2, scaler1, scaler2 = train_classifier(train_indexes, classifier_name, feature_set, reduction1, n_components1, reduction2, n_components2, with_rgb)
		pred_labels, scores, test_labels = test_classifier(test_indexes, classifier, feature_set, dimred1, dimred2, scaler1, scaler2, with_rgb)
		classifier_stats = get_classifier_stats(test_labels, pred_labels)
		stats.append(classifier_stats)

		# Compute ROC curve and area the curve
		fpr, tpr, thresholds = roc_curve(test_labels, scores[:, 1], pos_label=1)
		fpr_f.append(fpr)
		tpr_f.append(tpr)
		tprs.append(interp(mean_fpr, fpr, tpr))
		tprs[-1][0] = 0.0
		aucs.append(classifier_stats[auc_id()])

	mean_tpr = np.mean(tprs, axis=0)
	mean_tpr[-1] = 1.0
	std_auc = np.std(aucs)
	mean_stats = np.mean(stats, axis=0)
	return mean_stats

def cross_validate_with_prepared_folds(classifier_name, folds, train_data_per_fold, train_labels_per_fold, test_data_per_fold, test_labels_per_fold):

	tprs = []
	fpr_f = []
	tpr_f = []
	aucs = []
	mean_fpr = np.linspace(0, 1, 100)
	stats = []

	for f in range(folds): 
		train_data = train_data_per_fold[f]
		train_labels = train_labels_per_fold[f]
		test_data = test_data_per_fold[f]
		test_labels = test_labels_per_fold[f]
		classifier = get_classifier(classifier_name)
		classifier.fit(train_data, train_labels)
		pred_labels = classifier.predict(test_data)
		scores = classifier.predict_proba(test_data)

		classifier_params = classifier.get_params()
		if 'n_neighbors' in classifier_params:
			dist, neighbor_ids  = classifier.kneighbors(test_data, classifier_params['n_neighbors'])
			dist = np.mean(dist, axis=1)
			scores[:,0] = dist.transpose()
			scores[:,1] = (1 - dist).transpose()

		classifier_stats = get_classifier_stats(test_labels, pred_labels)
		stats.append(classifier_stats)

		# Compute ROC curve and area the curve
		fpr, tpr, thresholds = roc_curve(test_labels, scores[:, 1], pos_label=1)
		fpr_f.append(fpr)
		tpr_f.append(tpr)
		tprs.append(interp(mean_fpr, fpr, tpr))
		tprs[-1][0] = 0.0
		aucs.append(classifier_stats[auc_id()])

	mean_tpr = np.mean(tprs, axis=0)
	mean_tpr[-1] = 1.0
	std_auc = np.std(aucs)
	mean_stats = np.mean(stats, axis=0)
	return mean_stats

def apply_cross_validation(classifier_names, subset_name, stratified, feature_set, reduction1, n_components1, reduction2, n_components2, with_rgb=False, ignore_test=False):
	performance_stats = np.empty((0, 8))
	performance_names = []

	for classifier_name in classifier_names:

		try:
			cv_stats = cross_validate(classifier_name, subset_name, stratified, feature_set, reduction1, n_components1, reduction2, n_components2, with_rgb, ignore_test) 
			performance_stats = np.append(performance_stats, [cv_stats], axis=0)
			performance_names.append((",").join([classifier_name,subset_name,feature_set,reduction1, str(n_components1), reduction2, str(n_components2)]))
		except UserWarning:
			print("Variables are collinear in discriminant_analysis.")
		except RuntimeWarning: 
		 	print("Overflow or other errors.")
		except ValueError:
			print("Too few data for Dimension Reduction.")
	return performance_stats, performance_names

def apply_cross_validation_with_prepared_folds(classifier_names, subset_name, stratified, feature_set, reduction1, n_components1, reduction2, \
	n_components2, with_rgb=False, ignore_test=False):
	performance_stats = np.empty((0, 8))
	performance_names = []

	try:
		folds, train_data_per_fold, train_labels_per_fold, test_data_per_fold, test_labels_per_fold = get_fold_data(subset_name, stratified, \
			feature_set, reduction1, n_components1, reduction2, n_components2, with_rgb, ignore_test)

		for classifier_name in classifier_names:

			try:
				cv_stats = cross_validate_with_prepared_folds(classifier_name, folds, train_data_per_fold, train_labels_per_fold, test_data_per_fold, test_labels_per_fold) 
				performance_stats = np.append(performance_stats, [cv_stats], axis=0)
				performance_names.append((",").join([classifier_name,subset_name,feature_set,reduction1, str(n_components1), reduction2, str(n_components2)]))
			except RuntimeWarning: 
			 	print("Overflow or other errors.")
			except ValueError:
				print("Too few data for Dimension Reduction.")
			except UserWarning:
				print("Variables are collinear in discriminant_analysis.")
		
	except ConvergenceWarning:
		print('Dimred did not converge.')
	except UserWarning:
		print("n_components is too large.")


	return performance_stats, performance_names

def get_knn_classifiers_names(near_neighbors=None, distances=None, weights=None):
	if near_neighbors is None: 
		near_neighbors = {1, 3, 5}
	if distances is None: 
		distances = {'minkowski', 'euclidean', 'correlation', 'chebyshev', 'hamming', 'cityblock'}
	if weights is None: 
		weights = {'uniform', 'distance'}
	names = []
	for neighbors in near_neighbors: 
		for distance in distances:
			for weight in weights:
				names.append('-'.join(['KNN', str(neighbors), distance, weight]))
	return names

def get_svm_classifier_names(kernels=None, Cs=None, gammas=None, shrinkings=None, penalties=None):
	if kernels is None:
		kernels =  {'rbf', 'poly', 'linear', 'sigmoid'}
	if Cs is None: 
		Cs = {0.25, 0.5, 1.0, 2.0, 5.0}
	if gammas is None:
		gammas = {'auto', 'scale', 0.5, 1, 2}
	if shrinkings is None: 
		shrinkings = {True, False}
	if penalties is None:
		penalties = {'no_penalty', 'with_penalty', 'harsh_penalty'}
	names = []
	for kernel in kernels:
		for C in Cs:
			for gamma in gammas:
				for shrinking in shrinkings:
					for has_penalty in penalties:
						if not(kernel == 'linear' and gamma != 'auto'):
							names.append('-'.join(['SVM', kernel, str(C), str(gamma), str(shrinking), has_penalty]))
	return names

def get_rf_classifier_names(estimators=None, criteria=None, features=None, depth = None, penalties=None): 
	if estimators is None:
		estimators = {20, 50, 100}
	if criteria is None:
		criteria = {'gini', 'entropy'}
	if features is None: 
		features = {'none', 'sqrt', 'log2'}
	if depth is None: 
		depth = {'none', 20}
	if penalties is None:
		penalties = {'no_penalty', 'with_penalty', 'harsh_penalty'}
	names = []
	for n_estimators in estimators:
		for criterion in criteria:
				for max_features in features:
					for max_depth in depth:
						for has_penalty in penalties:
							names.append('-'.join(['Random Forest', str(n_estimators), criterion, max_features, str(max_depth), has_penalty]))
	return names


def get_validation_classifiers():
	return ["SVM-rbf-1-False",
			"KNN-5-correlation",
			"Random Forest-20-gini-sqrt",
			"SVM-sigmoid-0.5-True",
			"KNN-3-correlation",
			"Random Forest-50-entropy-auto",
			"SVM-sigmoid-1-True",
			"KNN-1-correlation",
			"Random Forest-50-gini-log2",
			"SVM-sigmoid-2-True",
			"KNN-3-minkowski",
			"Random Forest-50-entropy-sqrt", 
			"LDA", 
			"QDA"]

def get_validation_classifiers_noRF():
	clfs = get_svm_classifier_names({'rbf','linear', 'sigmoid'}, {0.5, 1.0, 2.0, 5.0}, {'auto', 'scale'}, {True, False}, {'no_penalty', 'with_penalty', 'harsh_penalty'})
	clfs.extend(get_knn_classifiers_names({1, 3, 5}, {'correlation', 'minkowski'}))
	clfs.extend({"LDA", "QDA"})
	return clfs

def get_validation_classifiers_withRF():
	clfs = get_svm_classifier_names({ 'rbf', 'poly'}, { 0.5, 1.0, 2.0}, {'auto'}, {True}, {'balanced','harsh_penalty'})
	clfs.extend(get_knn_classifiers_names({3, 5},  { 'euclidean', 'correlation', 'chebyshev'}, {'uniform', 'distance'}))
	clfs.extend(get_rf_classifier_names({20, 100}, {'gini', 'entropy'}, {'sqrt', 'log2'} , {'none'}, {'balanced','harsh_penalty'}))
	#clfs.append({"LDA", "QDA"})
	return clfs

def get_validation_classifiers_additional():
	clfs = get_svm_classifier_names({ 'rbf', 'linear'}, { 0.25, 0.5, 1.0, 2.0, 3.0}, {'auto'}, {True}, {'balanced','harsh_penalty'})
	clfs.extend(get_rf_classifier_names({ 200, 500}, {'gini', 'entropy'}, {'sqrt', 'log2'} , {'none'}, {'balanced','harsh_penalty'}))
	#clfs.append({"LDA", "QDA"})
	return clfs

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

def get_comparison_sets(input_sets=None, feature_sets=None, reduction1_sets=None, \
	reduction2_sets=None, n_components1_sets=None, n_components2_sets=None):

	if (input_sets == None):
		input_sets = {'unique', 'unique_unfixed', 'unique_fixed'}
	if (feature_sets == None):
		feature_sets = {'spect+clbp', 'spect+slbp', 'spect+mlbp', 'spect'}
	if (reduction1_sets == None):
		reduction1_sets = {'None', 'PCA', 'ICA'}		
	if (reduction2_sets == None):
		reduction2_sets = {'None', 'PCA', 'ICA'}
	if (n_components1_sets == None):
		n_components1_sets = {10, 20, 30, None}
	if (n_components2_sets == None):
		n_components2_sets = {10, 20, 30, None}

	print(input_sets, feature_sets, reduction1_sets, reduction2_sets, n_components1_sets, n_components2_sets)
	return input_sets, feature_sets, reduction1_sets, reduction2_sets, n_components1_sets, n_components2_sets


def compare_validation_performance(classifier_names, title, stratified=True, with_rgb=False, ignore_test=False, \
	input_sets=None, feature_sets=None, reduction1_sets=None, reduction2_sets=None, n_components1_sets=None, n_components2_sets=None):
	print(classifier_names)
	
	input_sets, feature_sets, reduction1_sets, reduction2_sets, n_components1_sets, n_components2_sets = get_comparison_sets(input_sets, feature_sets, \
		reduction1_sets,reduction2_sets, n_components1_sets, n_components2_sets)

	for input_set in input_sets: 
		for feature_set in feature_sets: 
			for reduction1 in reduction1_sets:
				for n_components1 in n_components1_sets: 
					for reduction2 in reduction2_sets:
						for n_components2 in n_components2_sets:
							if (not(reduction1 is 'None' and n_components1 is not None)) \
							 and (not(reduction2 is 'None' and n_components2 is not None)) \
							 and (not(input_set is 'spect' and reduction2 is not 'None' and n_components2 is not None)):
								start_time = time.time()									
								clf_stats, clf_names = apply_cross_validation_with_prepared_folds(classifier_names, input_set, stratified, feature_set, reduction1, n_components1, reduction2, n_components2, with_rgb, ignore_test)
								#clf_stats, clf_names = apply_cross_validation(classifier_names, input_set, stratified, feature_set, reduction1, n_components1, reduction2, n_components2, with_rgb, ignore_test)
								[dh.write_log(classification_log, ",".join([x, ",".join([str(yy) for yy in y]) ,"\n"])) for (x, y) in zip(clf_names, clf_stats)]
								print("--- %s seconds ---" % (time.time() - start_time))


								# except UserWarning:
								#  	print('Too many components.')

	print('Finished.')



def compare_testing_performance(input_set="unique", with_rgb=False,\
	feature_set=None, reduction1=None, reduction2=None, n_components1=None, n_components2=None):


	test_indexes = dh.get_test_indexes(input_set)
	print("Test indexes", test_indexes)
	train_indexes = [x for x in   dh.subset_indexes(input_set) if x not in test_indexes]
	print("Train indexes", train_indexes)
	print('\n')

	for classifier_name in get_testing_classifiers():
		classifier, dimred1, dimred2, scaler1, scaler2 = train_classifier(train_indexes, classifier_name, feature_set, reduction1, n_components1, reduction2, n_components2, with_rgb)
		pred_labels, scor, test_labels = test_classifier(test_indexes, classifier, feature_set, dimred1, dimred2, scaler1, scaler2, with_rgb)
		classifier_stats = get_classifier_stats(test_labels, pred_labels) 
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
		scr = scor[:,1].transpose()
		pred = pred_labels
		dh.write_log(classification_log, ",".join([conf, str(auc), str(balAc) , str(acc), str(sensitivity),  str(specificity), str(doh), '\n']))
		print(conf)
		print('AUC', auc, 'BalAcc', balAc , 'Accuracy', acc, 'sensitivity', sensitivity, 'specificity', specificity, 'DOH', doh)
		print('Scores', scr)
		print('Predictions', pred)
		print('True labels', test_labels, '\n')


current_run_case = 'Validation MSI'
column_names = (",").join([current_run_case, 'Input', 'Feature', 'Dimred1', 'NComp1', 'Dimred2', 'NComp2', 'Accuracy', 'AUC', \
				'Specificity', 'Sensitivity', 'FPR', 'FNR', 'BalancedAccuracy', 'F1', '\n' ])

classification_log = dh.get_log_file()
#dh.write_log(classification_log, column_names)
#get_validation_classifiers_RF_SVM
#compare_validation_performance(get_validation_classifiers_additional(), 'Validation_Classifiers', True, False, False, \
#	{ 'unique_unfixed', 'unique_fixed'},  { 'spect+clbp', 'spect+slbp', 'spect+mlbp', 'spect' }, \
#	{ 'PCA', 'ICA'}, { 'PCA', 'ICA'}, {10, 20}, {10, 20})

#classification_log = dh.get_log_file()
#dh.write_log(classification_log, column_names)
#compare_validation_performance(get_validation_classifiers_additional(), 'Validation_Classifiers', True, True, False, \
#	{'unique', 'unique_unfixed', 'unique_fixed'},  { 'spect+slbp', 'spect' }, \
#	{ 'PCA', 'ICA'}, { 'PCA', 'ICA'}, {10, 20}, {10, 20})

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

compare_testing_performance('unique', False, 'spect+mlbp', 'PCA', 'ICA',  20, 20)
print('\n\n')
compare_testing_performance('unique', False, 'spect+mlbp', 'PCA', 'PCA',  20, 20)
print('\n\n')
compare_testing_performance('unique', False, 'spect+mlbp', 'ICA', 'ICA',  20, 20)
print('\n\n')
compare_testing_performance('unique', False, 'spect+mlbp', 'ICA', 'PCA',  20, 20)
print('\n\n')

