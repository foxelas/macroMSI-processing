import data_handling as dh
import dimred as dm

from os.path import dirname, join as pjoin, exists
from sklearn.preprocessing import StandardScaler
from sklearn import manifold
import sklearn.metrics
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA, QuadraticDiscriminantAnalysis as QDA
from sklearn.neighbors import KNeighborsClassifier as KNN
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.metrics import confusion_matrix, accuracy_score, auc, roc_auc_score, roc_curve
import numpy as np
from scipy import interp
import scipy.io as sio 
from scipy.spatial.distance import chebyshev, correlation, hamming, minkowski
import matplotlib.pyplot as plt
from matplotlib import interactive, is_interactive
import sys

interactive(True)

out_dir = dh.get_out_dir()
classification_log = dh.get_log_file()
#sys.stdout = open(classification_log, 'w+')

####################classification#######################

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
	accuracy = performance_stats[:,0]
	auc_score = performance_stats[:,1]

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
		a, b, c, d = classifier_name.split('-') 
		c = c if (c == 'auto' or c == 'scale') else float(c)
		clf = SVC(kernel=b, gamma=c, shrinking=(d==True), probability=True)

	elif 'KNN' in classifier_name:
		a, b, c =  classifier_name.split('-')
		clf = KNN(n_neighbors=int(b), algorithm = 'auto', metric=c)

	elif 'QDA' in classifier_name:
		clf = QDA()

	elif 'LDA' in classifier_name:
		clf = LDA()

	elif 'Decision Tree' in classifier_name:
		a, b=  classifier_name.split('-')
		clf = DecisionTreeClassifier(criterion=b)

	elif 'Random Forest' in classifier_name:
		clf = RandomForestClassifier(n_estimators=10)

	elif 'Naive Bayes' in classifier_name:
		clf = GaussianNB()

	else:
		print(classifier_name, 'Not yet implemented.')
		return

	return clf

def train_classifier(train_indexes, classifier_name, has_texture=False, n_components=10):
	data = dh.get_reconstructed_spectra()
	labels = dh.get_labels()
	reference = dh.get_measured_spectra()

	clf = get_classifier(classifier_name)
	if '/' in classifier_name:
		dimred =  dm.get_dimred(classifier_name[classifier_name.find("/")+1:], n_components)
		dimred2 =  dm.get_dimred(classifier_name[classifier_name.find("/")+1:], 20)
	else:
		dimred = None

	scaler1 = StandardScaler()
	reference_data = scaler1.fit_transform(reference[train_indexes, :])
	train_data, train_labels, train_lbp = dh.get_scaled_subset_with_index(train_indexes, data, labels, scaler1, has_texture)

	scaler2 = StandardScaler()
	if dimred is not None: 
		dimred = dimred.fit(reference_data)
		train_data = dimred.transform(train_data)
		#print("PCA-explained variance for specrum", dimred.explained_variance_ratio_ * 100)

		if has_texture:
			dimred2 = dimred2.fit(train_lbp)
			train_lbp = dimred2.transform(train_lbp)
			#print("PCA-explained variance for texture", dimred2.explained_variance_ratio_ * 100)
	else: 
		dimred2 = None

	train_data = dh.concat_features(train_data, train_lbp)
	train_data = scaler2.fit_transform(train_data)

	clf.fit(train_data, train_labels)

	return clf, dimred, scaler1, scaler2, dimred2

def test_classifier(test_indexes, classifier, dimred, scaler1, scaler2, dimred2, has_texture=False):
	data = dh.get_reconstructed_spectra()
	labels = dh.get_labels()

	test_data, test_labels, test_lbp = dh.get_scaled_subset_with_index(test_indexes, data, labels, scaler1, has_texture)
	print(test_indexes)

	if dimred is not None: 
		test_data = dimred.transform(test_data)
		if has_texture:
			test_lbp = dimred2.transform(test_lbp)

	test_data = dh.concat_features(test_data, test_lbp)

	test_data = scaler2.transform(test_data)

	predictions = classifier.predict(test_data)
	scores = classifier.predict_proba(test_data)

	return predictions, scores, test_labels

def run_classification_test(subset_name, classifier_name, has_texture=False, n_components=10):
	fold_indexes, folds = dh.get_fold_indexes(subset_name, 10)
	train_indexes = [y for x in fold_indexes for y in x ]
	classifier, dimred, scaler1, scaler2, dimred2 = train_classifier(train_indexes, classifier_name, has_texture, n_components)
	test_indexes = dh.get_test_indexes(subset_name)
	pred_labels, scores, test_labels = test_classifier(test_indexes, classifier, dimred, scaler1, scaler2, dimred2, has_texture)

	acc = accuracy_score(test_labels, pred_labels)
	auc = roc_auc_score(test_labels, pred_labels)
	print("True labels", test_labels)
	print("Pred labels", pred_labels)
	print("Pred scores", scores)
	print("Accuracy", acc)
	print("A.U.C", auc)
	return pred_labels, scores, acc, auc

def cross_validate(classifier_name, subset_name, stratified=True, n_components=10, has_texture=False):

	tprs = []
	aucs = []
	accs = []
	fpr_f = []
	tpr_f = []
	mean_fpr = np.linspace(0, 1, 100)

	if stratified: 
		fold_indexes, folds = dh.get_fold_indexes_stratified(subset_name, 10)
	else:
		fold_indexes, folds = dh.get_fold_indexes(subset_name, 10)

	for f in range(folds): 
		test_indexes = fold_indexes[f]
		train_indexes = [y for x in fold_indexes if x != test_indexes for y in x ]

		classifier, dimred, scaler1, scaler2, dimred2 = train_classifier(train_indexes, classifier_name, has_texture, n_components)
		pred_labels, scores, test_labels = test_classifier(test_indexes, classifier, dimred, scaler1, scaler2, dimred2, has_texture)

		accs.append(accuracy_score(test_labels, pred_labels))
		# Compute ROC curve and area the curve
		fpr, tpr, thresholds = roc_curve(test_labels, scores[:, 1])
		fpr_f.append(fpr)
		tpr_f.append(tpr)
		tprs.append(interp(mean_fpr, fpr, tpr))
		tprs[-1][0] = 0.0
		roc_auc = auc(fpr, tpr)
		#aucs.append(roc_auc)
		aucs.append(roc_auc_score(test_labels, pred_labels))

	mean_tpr = np.mean(tprs, axis=0)
	mean_tpr[-1] = 1.0
	mean_auc = auc(mean_fpr, mean_tpr)
	std_auc = np.std(aucs)

	mean_acc = np.mean(accs)

	return mean_acc, mean_auc 

def classify_many(classifier_names, subset_name, stratified=True, n_components=10, has_texture=False):
	performance_stats = np.empty((0, 2))
	performance_names = []

	data = dh.get_reconstructed_spectra()
	labels = dh.get_labels()
	reference = dh.get_measured_spectra()

	for classifier_name in classifier_names:
		mean_acc, mean_auc = cross_validate(classifier_name, subset_name, stratified, n_components, has_texture)
		print(classifier_name, ' with accuracy ', mean_acc,', auc ', mean_auc, '\n', flush=True)
		performance_stats = np.append(performance_stats, [(mean_acc, mean_auc)], axis=0)

	return performance_stats


def get_best_classifiers(performance_stats, performance_names, title, best_n=10):
	
	mutual_performance_stats = [ x[0] + x[1] for x in performance_stats]
	print(mutual_performance_stats)
	best_performance_indexes = np.argsort(mutual_performance_stats)[-best_n:]
	best_performance_names = [performance_names[i] for i in best_performance_indexes]
	best_performance_stats = performance_stats[ best_performance_indexes, :]

	if 'Input' in title: 
		#sort by last word in classifier name 
		best_classifier_ids_sorted = np.argsort([x[::-1] for x in best_performance_names])
		best_performance_names = [best_performance_names[x] for x in best_classifier_ids_sorted]
		best_performance_stats = best_performance_stats[best_classifier_ids_sorted,:]
	
	plot_accuracy_and_auc(best_performance_stats, best_performance_names, title)

	print('Best ' + title, flush=True)
	[print('Accuracy:', y, ', AUC:', z, '--', x, flush=True) for (x,y,z) in zip(best_performance_names, best_performance_stats[:,0], best_performance_stats[:,1])]

	return best_performance_names

def get_various_classifiers_names():
	return ["KNN-3-minkowski", "KNN-1-correlation", "SVM-linear-0.025-True", "SVM-linear-auto-True", "SVM-rbf-scale-False", "Decision Tree-gini",
         "Decision Tree-entropy", "Random Forest", "Naive Bayes", "LDA", "QDA"]

def get_knn_classifiers_names():
	return ["KNN-1-minkowski","KNN-3-minkowski","KNN-5-minkowski",
		"KNN-1-cityblock","KNN-3-cityblock","KNN-5-cityblock",
		"KNN-1-euclidean","KNN-3-euclidean","KNN-5-euclidean",
		"KNN-1-chebyshev","KNN-3-chebyshev","KNN-5-chebyshev",
		"KNN-1-hamming","KNN-3-hamming","KNN-5-hamming",
		"KNN-1-correlation","KNN-3-correlation","KNN-5-correlation"]

def get_svm_classifier_names(): 
	names = []
	for kernel in ('rbf', 'poly', 'linear', 'sigmoid'):
		for gamma in {'auto', 'scale', 0.5, 1, 2}:
				for shrinking in (True, False):
					if not(kernel == 'linear' and gamma != 'auto'):
						names.append('-'.join(['SVM', kernel, str(gamma), str(shrinking)]))
	return names

def compare_classifiers(names, subset_name, title='', n_components=10, stratified=True, has_texture=False):

	classifier_names = names + ['/'.join([c, 'PCA']) for c in names]
	performance_stats = classify_many(classifier_names, subset_name, stratified, n_components, has_texture)
	best_classifier_names = get_best_classifiers(performance_stats, classifier_names, title, 20)

	return best_classifier_names

def compare_input_sets(classifier_names, title, n_components=10, stratified=True, has_texture=False):

	performance_stats = np.empty((0, 2))
	performance_names = []
	for input_set in ('unique_unfixed', 'unique_fixed', 'unique_cut', 'unique'):
		dh.write_file(classification_log, "Input set: " + input_set)
		pf_stats1 = classify_many(classifier_names, input_set, stratified, n_components, has_texture)
		pf_names1 = [ ("/").join([x, input_set]) for x in classifier_names]
		performance_stats = np.concatenate((performance_stats, pf_stats1), axis=0)
		performance_names = performance_names + pf_names1

	best_classifier_names = get_best_classifiers(performance_stats, performance_names, title, 20)

def run_comparison():
	input_set = "unique"
	print("Classification comparison", '\n')
	compare_classifiers(get_svm_classifier_names(), input_set, 'SVM Classification Performance Comparison' + '(' + input_set + ')', 10, True)
	compare_classifiers(get_knn_classifiers_names(), input_set, 'KNN Classification Performance Comparison' + '(' + input_set + ')', 10, True)
	compare_classifiers(get_various_classifiers_names(), input_set, 'Classification Performance Comparison' + '(' + input_set + ')', 10, True)
	compare_input_sets(get_various_classifiers_names(), 'Input Set Comparison', 10, True)

run_comparison()

# img_spectra_mat = dh.loadmat( pjoin(out_dir, 'img_spectra.mat'))
# print(img_spectra_mat.keys())
# img_spectra = img_spectra_mat['estimated_2D']
# predictions, scores = get_predictions("KNN-3-minkowski/PCA", 'unique', img_spectra)
# dh.savemat(pjoin(out_dir, 'predictions.mat'), mdict = {"predictions": predictions, "scores": scores})

predictions, scores, acc, auc = run_classification_test('unique', "KNN-3-minkowski/PCA", True)

