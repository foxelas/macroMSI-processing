import data_handling as dh
import dimred as dm
import classifier_settings as cs

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
import datetime

####################Settings#######################
warnings.filterwarnings('error')

interactive(True)
out_dir = dh.get_out_dir()
isStratified = True  # False 
isRGB = False  
ignore_test = False 

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

def get_fold_data(subset_name, stratified, feature_set, reduction1, n_components1, reduction2, n_components2, with_rgb=False, ignore_test=False, split_to_folds = 5):
	if stratified: 
		fold_indexes, folds = dh.get_fold_indexes_stratified(subset_name, split_to_folds, ignore_test)
	else:
		fold_indexes, folds = dh.get_fold_indexes(subset_name, split_to_folds, ignore_test)

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

# def cross_validate(classifier_name, subset_name, stratified, feature_set, reduction1, n_components1, reduction2, n_components2, with_rgb=False, ignore_test=False):

# 	tprs = []
# 	fpr_f = []
# 	tpr_f = []
# 	aucs = []
# 	mean_fpr = np.linspace(0, 1, 100)
# 	stats = []

# 	if stratified: 
# 		fold_indexes, folds = dh.get_fold_indexes_stratified(subset_name, 5, ignore_test)
# 	else:
# 		fold_indexes, folds = dh.get_fold_indexes(subset_name, 5, ignore_test)

# 	for f in range(folds): 
# 		test_indexes = fold_indexes[f]
# 		train_indexes = [y for x in fold_indexes if x != test_indexes for y in x ]

# 		classifier, dimred1, dimred2, scaler1, scaler2 = train_classifier(train_indexes, classifier_name, feature_set, reduction1, n_components1, reduction2, n_components2, with_rgb)
# 		pred_labels, scores, test_labels = test_classifier(test_indexes, classifier, feature_set, dimred1, dimred2, scaler1, scaler2, with_rgb)
# 		classifier_stats = get_classifier_stats(test_labels, pred_labels)
# 		stats.append(classifier_stats)

# 		# Compute ROC curve and area the curve
# 		fpr, tpr, thresholds = roc_curve(test_labels, scores[:, 1], pos_label=1)
# 		fpr_f.append(fpr)
# 		tpr_f.append(tpr)
# 		tprs.append(interp(mean_fpr, fpr, tpr))
# 		tprs[-1][0] = 0.0
# 		aucs.append(classifier_stats[auc_id()])

# 	mean_tpr = np.mean(tprs, axis=0)
# 	mean_tpr[-1] = 1.0
# 	std_auc = np.std(aucs)
# 	mean_stats = np.mean(stats, axis=0)
# 	return mean_stats

def cross_validate_with_prepared_folds(classifier_name, folds, train_data_per_fold, train_labels_per_fold, test_data_per_fold, test_labels_per_fold, save_fold_stats=False):

	if save_fold_stats: 	
		prefix  = "RGB" if isRGB else "MSI"
		folds_filename = pjoin(out_dir, "classification_logs", \
			prefix + "_validation_folds_" + dh.get_scales_in_use() + "scales" + datetime.datetime.now().strftime("%Y-%m-%d %H_%M") + ".csv")

		column_names = [ "Accuracy", "AUC", "Specificity", "Sensitivity", "FPR", "FNR", "BalancedAccuracy", "F1" ]
		dh.write_csv(folds_filename, column_names)

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

		if isinstance(classifier_name, str):
			classifier = cs.get_classifier(classifier_name)
		else:
			classifier = classifier_name

		classifier.fit(train_data, train_labels)
		pred_labels = classifier.predict(test_data)
		scores = classifier.predict_proba(test_data)

		classifier_stats = get_classifier_stats(test_labels, pred_labels)
		stats.append(classifier_stats)

		# Compute ROC curve and area the curve
		fpr, tpr, thresholds = roc_curve(test_labels, scores[:, 1], pos_label=1)
		fpr_f.append(fpr)
		tpr_f.append(tpr)
		tprs.append(interp(mean_fpr, fpr, tpr))
		tprs[-1][0] = 0.0
		aucs.append(classifier_stats[auc_id()])

		if save_fold_stats:
			dh.append_csv(folds_filename, classifier_stats)
			plot_cv_roc_curve(fpr_f, tpr_f, tprs, aucs)


	mean_tpr = np.mean(tprs, axis=0)
	mean_tpr[-1] = 1.0
	std_auc = np.std(aucs)
	mean_stats = np.mean(stats, axis=0)
	return mean_stats

# def apply_cross_validation(classifier_names, subset_name, stratified, feature_set, reduction1, n_components1, reduction2, n_components2, with_rgb=False, ignore_test=False):
# 	performance_stats = np.empty((0, 8))
# 	performance_names = []

# 	for classifier_name in classifier_names:

# 		try:
# 			cv_stats = cross_validate(classifier_name, subset_name, stratified, feature_set, reduction1, n_components1, reduction2, n_components2, with_rgb, ignore_test) 
# 			performance_stats = np.append(performance_stats, [cv_stats], axis=0)
# 			if not isinstance(classifier_name, str):
# 				if "criterion" in classifier_name.get_params():
# 					classifier_name = "Random Forest"
# 				elif "kernel" in classifier_name.get_params():
# 					classifier_name = "SVM"
# 				elif "n_neighbors" in classifier_name.get_params():
# 					classifier_name = "KNN"
			
# 			performance_names.append((",").join([classifier_name,subset_name,feature_set,reduction1, str(n_components1), reduction2, str(n_components2)]))
# 		except UserWarning:
# 			print("Variables are collinear in discriminant_analysis.")
# 		except RuntimeWarning: 
# 		 	print("Overflow or other errors.")
# 		except ValueError:
# 			print("Too few data for Dimension Reduction.")
# 	return performance_stats, performance_names

def apply_cross_validation_with_prepared_folds(classifier_names, subset_name, stratified, feature_set, reduction1, n_components1, reduction2, \
	n_components2, with_rgb=False, ignore_test=False, split_to_folds = None, save_fold_stats = False):
	performance_stats = np.empty((0, 8))
	performance_names = []

	try:
		folds, train_data_per_fold, train_labels_per_fold, test_data_per_fold, test_labels_per_fold = get_fold_data(subset_name, stratified, \
			feature_set, reduction1, n_components1, reduction2, n_components2, with_rgb, ignore_test, split_to_folds)

		for classifier_name in classifier_names:

			try:
				cv_stats = cross_validate_with_prepared_folds(classifier_name, folds, train_data_per_fold, train_labels_per_fold,\
				 test_data_per_fold, test_labels_per_fold, save_fold_stats) 
				performance_stats = np.append(performance_stats, [cv_stats], axis=0)
				if not isinstance(classifier_name, str):
					if "criterion" in classifier_name.get_params():
						classifier_name = "Random Forest"
					elif "kernel" in classifier_name.get_params():
						classifier_name = "SVM"
					elif "n_neighbors" in classifier_name.get_params():
						classifier_name = "KNN"
			
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

	

def find_optimal_hyperparameters_with_grid_search(classifier_name, train_data, train_labels):
	grid_search_randomizer = cs.get_grid_search_randomizer(classifier_name)
	grid_search_randomizer.fit(train_data, train_labels)
	print("Best params for ", classifier_name)
	print(grid_search_randomizer.best_params_)

def run_find_optimal_hyperparameters_with_grid_search(input_sets=None, feature_sets=None, reduction1_sets=None, reduction2_sets=None, \
	n_components1_sets=None, n_components2_sets=None):
	#test_indexes = dh.get_test_indexes(input_set)
	#train_indexes = [x for x in   dh.subset_indexes(input_set) if x not in test_indexes]

	input_sets, feature_sets, reduction1_sets, reduction2_sets, n_components1_sets, n_components2_sets = get_comparison_sets(input_sets, feature_sets, \
		reduction1_sets,reduction2_sets, n_components1_sets, n_components2_sets)

	for input_set in input_sets: 
		for feature_set in feature_sets: 
			for reduction1 in reduction1_sets:
				for n_components1 in n_components1_sets: 
					for reduction2 in reduction2_sets:
						for n_components2 in n_components2_sets:
							if not((reduction1 is "None" and n_components1 is not None)\
								or (reduction2 is "None" and n_components2 is not None)\
								or (input_set is "spect" and reduction2 is not "None" and n_components2 is not None)):


								train_indexes, folds = dh.get_fold_indexes(input_set, 1, ignore_test)
								train_indexes = train_indexes[0]
								#print(train_indexes)
								test_indexes = [1]
								train_data, train_labels, test_data, test_labels = get_scaled_data(train_indexes, test_indexes, feature_set, reduction1, \
								n_components1, reduction2, n_components2, isRGB)

								for classifier_name in ["Random Forest", "SVM", "KNN"]: 
									print("Optimizing for ", feature_set, reduction1, n_components1, reduction2, n_components2, classifier_name )
									find_optimal_hyperparameters_with_grid_search(classifier_name, train_data, train_labels)



def compare_validation_performance(classifier_names, title, stratified=True, with_rgb=False, ignore_test=False, \
	input_sets=None, feature_sets=None, reduction1_sets=None, reduction2_sets=None, n_components1_sets=None, n_components2_sets=None, split_to_folds = None):
	
	input_sets, feature_sets, reduction1_sets, reduction2_sets, n_components1_sets, n_components2_sets = get_comparison_sets(input_sets, feature_sets, \
		reduction1_sets,reduction2_sets, n_components1_sets, n_components2_sets)
	for input_set in input_sets: 
		for feature_set in feature_sets: 
			for reduction1 in reduction1_sets:
				for n_components1 in n_components1_sets: 
					for reduction2 in reduction2_sets:
						for n_components2 in n_components2_sets:
							if not((reduction1 is "None" and n_components1 is not None)\
								or (reduction2 is "None" and n_components2 is not None)\
								or (input_set is "spect" and reduction2 is not "None" and n_components2 is not None)):

								start_time = time.time()			 
								actual_classifier_names = classifier_names
								clf_stats, clf_names = apply_cross_validation_with_prepared_folds(actual_classifier_names, input_set, stratified, feature_set, \
								 reduction1, n_components1, reduction2, n_components2, with_rgb, ignore_test, split_to_folds)
								
								[dh.append_csv(validation_filename, [x] + y.tolist()) for (x,y) in zip(clf_names, clf_stats)]
								print("--- %s seconds ---" % (time.time() - start_time))


								# except UserWarning:
								#  	print('Too many components.')

	print('Finished.')

def compare_validation_performance_of_optimized(classifier_names, title, stratified=True, with_rgb=False, ignore_test=False, \
	input_sets=None, feature_sets=None, reduction1_sets=None, reduction2_sets=None, n_components1_sets=None, n_components2_sets=None, \
	split_to_folds=None):
	
	isStringName = isinstance(classifier_names[0], str)
	input_sets, feature_sets, reduction1_sets, reduction2_sets, n_components1_sets, n_components2_sets = get_comparison_sets(input_sets, feature_sets, \
		reduction1_sets,reduction2_sets, n_components1_sets, n_components2_sets, split_to_folds)
	i = 0
	for input_set in input_sets: 
		for feature_set in feature_sets: 
			for reduction1 in reduction1_sets:
				for n_components1 in n_components1_sets: 
					for reduction2 in reduction2_sets:
						for n_components2 in n_components2_sets:

							if not((reduction1 is "None" and n_components1 is not None)\
								or (reduction2 is "None" and n_components2 is not None)\
								or (input_set is "spect" and reduction2 is not "None" and n_components2 is not None)):

								for classifier_case in ["Random Forest", "SVM", "KNN"]: 

									start_time = time.time()			
									if isStringName : 
										actual_classifier_names = classifier_names
									else :						
										actual_classifier_names = [ classifier_names[i] ]
										i += 1
									clf_stats, clf_names = apply_cross_validation_with_prepared_folds(actual_classifier_names, input_set, stratified, feature_set, reduction1, n_components1, reduction2, n_components2, with_rgb, ignore_test)
									clf_names.extend(clf_stats[0])
									dh.append_csv(validation_filename, clf_names)
									print("--- %s seconds ---" % (time.time() - start_time))


								# except UserWarning:
								#  	print('Too many components.')

	print('Finished.')

def get_fold_info(classifier_name, stratified, with_rgb, ignore_test, \
	input_set, feature_set, reduction1, n_components1, reduction2, n_components2, split_to_folds = None):

	start_time = time.time()			 
	clf_stats, clf_names = apply_cross_validation_with_prepared_folds([classifier_name], input_set, stratified, feature_set, \
	 reduction1, n_components1, reduction2, n_components2, with_rgb, ignore_test, split_to_folds)
	
	[dh.append_csv(validation_filename, [x] + y.tolist()) for (x,y) in zip(clf_names, clf_stats)]
	print("--- %s seconds ---" % (time.time() - start_time))

	print('Finished.')


# run_find_optimal_hyperparameters_with_grid_search(["unique", "unique_unfixed", "unique_fixed"],  [ "spect+clbp", "spect+slbp", "spect+mlbp", "spect" ], \
# 	[ "None", "PCA", "ICA"], [ "None", "PCA", "ICA"], [20, None], [20, None])

prefix  = "RGB" if isRGB else "MSI"
validation_filename = pjoin(out_dir, "classification_logs", \
	prefix + "_validation_" + dh.get_scales_in_use() + "scales" + datetime.datetime.now().strftime("%Y-%m-%d %H_%M") + ".csv")

column_names = ["Classifier", "Input", "Feature", "Dimred1", "NComp1", "Dimred2", "NComp2", "Accuracy", "AUC", \
				"Specificity", "Sensitivity", "FPR", "FNR", "BalancedAccuracy", "F1" ]

dh.write_csv(validation_filename, column_names)
print(validation_filename)

# compare_validation_performance_of_optimized(cs.get_optimized_classifiers(), 'Validation_Classifiers', isStratified, isRGB, ignore_test, \
# 	["unique", "unique_unfixed", "unique_fixed"],  [ "spect+clbp", "spect+slbp", "spect+mlbp", "spect" ], \
#  	[ "None", "PCA", "ICA"], [ "None", "PCA", "ICA"], [20, None], [20, None], 2)

# compare_validation_performance(get_validation_classifiers_withRF(), 'Validation_Classifiers', isStratified, isRGB, ignore_test, \
# 	["unique", "unique_unfixed", "unique_fixed"],  [ "spect+mlbp", "spect" ], \
#  	[ "None", "PCA", "ICA"], [ "None", "PCA", "ICA"], [20, None], [20, None])

#get_fold_info("Random Forest-100-gini-log2-none-balanced", isStratified, isRGB, ignore_test, "unique", "spect+mlbp", "ICA", 20, "ICA", 20, 5)

#classification_log = dh.get_log_file()
#dh.write_log(classification_log, column_names)
#compare_validation_performance(get_validation_classifiers_additional(), 'Validation_Classifiers', True, True, False, \
#	{'unique', 'unique_unfixed', 'unique_fixed'},  { 'spect+slbp', 'spect' }, \
#	{ 'PCA', 'ICA'}, { 'PCA', 'ICA'}, {10, 20}, {10, 20})