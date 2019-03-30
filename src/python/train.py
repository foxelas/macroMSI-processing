import data_handling as dh
import dimred as dm

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


interactive(True)

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

	create_directory(pjoin(out_dir, 'classification'))
	img_title = plot_title.replace(' ', '_').lower()
	plt.savefig(pjoin(out_dir, 'classification', img_title + '.png'), bbox_inches='tight')

def get_classifier(classifier_name):
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

def cross_validate(classifier, data, labels, subset_name, reference=None, dimred=None, stratified=True):

	if reference is None: 
		reference = data

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
		train_indexes = []
		train_indexes = [y for x in fold_indexes if x != test_indexes for y in x ]

		sc = StandardScaler()
		reference_data = sc.fit_transform(reference[train_indexes, :])
		train_data, train_labels = dh.get_scaled_subset_with_index(train_indexes, data, labels, sc)
		test_data, test_labels = dh.get_scaled_subset_with_index(test_indexes, data, labels, sc)

		if dimred is not None: 
			dimred = dimred.fit(reference_data)
			train_data = dimred.transform(train_data)
			test_data = dimred.transform(test_data)

			#needed or not?
			sc2 = StandardScaler()
			train_data = sc2.fit_transform(train_data)
			test_data = sc2.fit(test_data)

		classifier.fit(train_data, train_labels)
		pred_labels = classifier.predict(test_data)
		accs.append(accuracy_score(test_labels, pred_labels))
		pred_probas = classifier.predict_proba(test_data)
		# Compute ROC curve and area the curve
		fpr, tpr, thresholds = roc_curve(test_labels, pred_probas[:, 1])
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


def classify_many(names, subset_name, n_components, stratified):
	performance_stats = np.empty((0, 2))
	performance_names = []

	data = dh.get_reconstructed_spectra()
	labels = dh.get_labels()
	reference = dh.get_measured_spectra()

	for name in names:
		clf = get_classifier(name)
		if name.contains("/"): 
			dimred =  dm.get_dimred(classifier_name[classifier_name.find("/"):], n_components)

		mean_acc, mean_auc = cross_validate(clf, data, labels, subset_name, reference, dimred, stratified)
		performance_stats = np.append(performance_stats, [(mean_acc, mean_auc)], axis=0)

	return performance_stats


def get_best_classifiers(performance_stats, performance_names, title, best_n=10):
	
	print('\n------------------------------------------------------')
	print(performance_names, ': Best performance \n')
	max_accuracy_index = (np.argsort(performance_stats[:,0]))[-best_n:]
	[print(-i, 'th best Accuracy ', x,  ' by ',y) for i,x,y in zip(range(-len(max_accuracy_index), 0), performance_stats[max_accuracy_index,0], [performance_names[i] for i in max_accuracy_index])]
	print(' ')
	max_auc_index = np.argsort(performance_stats[:,1])[-best_n:]
	[print(-i, 'th best A.U.C. ', x,  ' by ',y) for i,x,y in zip(range(-len(max_accuracy_index), 0), performance_stats[max_auc_index,1], [performance_names[i] for i in max_auc_index])]
	print('------------------------------------------------------')

	#concatenate best AUC and best accuracy performances
	best_classifier_names = [performance_names[i] for i in  max_accuracy_index] + [performance_names[i] for i in max_auc_index]
	best_classifier_stats = np.concatenate(([performance_stats[i,:] for i in  max_accuracy_index], [performance_stats[i,:] for i in max_auc_index]), axis=0)
	#return only unique classifications
	best_classifier_stats, unique_ids = np.unique(best_classifier_stats, return_index=True, axis=0) 
	best_classifier_names = [best_classifier_names[i] for i in unique_ids]

	if 'Input' in title: 
		#sort by last word in classifier name 
		best_classifier_ids_sorted = np.argsort([x[::-1] for x in best_classifier_names])
		best_classifier_names = [best_classifier_names[x] for x in best_classifier_ids_sorted]
		best_classifier_stats = best_classifier_stats[best_classifier_ids_sorted,:]
		plot_accuracy_and_auc(best_classifier_stats, best_classifier_names, title)
	else:
		plot_accuracy_and_auc(performance_stats, performance_names, title)

	print('Best ' + title)
	[print('Accuracy:', y, ', AUC:', z, '--', x) for (x,y,z) in zip(best_classifier_names, best_classifier_stats[:,0], best_classifier_stats[:,1])]

	return best_classifier_names

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

def compare_classifiers(names, subset_name, title='', n_components=10, stratified=True):

	classifier_names = names
	classifier_names = [classifier_names.append(("/").join([clf, "PCA"])) for clf in names]

	performance_stats = classify_many(classifier_names, subset_name, n_components, stratified)
	best_classifier_names = get_best_classifiers(performance_stats, classifier_names, title, 10)

	return best_classifier_names

def compare_input_sets(names, data, labels, title):
	classifiers = [get_classifier(c) for c in names]

	performance_stats = np.empty((0, 2))
	performance_names = []
	for input_set in ('unique_unfixed', 'unique_fixed', 'unique_cut', 'unique'):
		data_s, labels_s, fixation_s = get_subset(input_set, data, labels)
		print('Subset ', input_set,' contains ', len(labels_s), ' observations')
		pf_stats1, pf_names1 = classify_many([ '-'.join([x, input_set]) for x in names], classifiers, data_s, labels_s, '', 5)

		performance_stats = np.concatenate((performance_stats, pf_stats1), axis=0)
		performance_names = performance_names + pf_names1

	best_classifier_names = get_best_classifiers(performance_stats, performance_names, title, 15)

def split_and_scale(data, labels):
	sc = StandardScaler()
	train_data, test_data, train_labels, test_labels = train_test_split(data, labels, test_size=0.25, random_state=1)
	train_data = sc.fit_transform(train_data)
	test_data = sc.transform(test_data)
	return train_data, test_data, train_labels, test_labels, sc

def apply_on_dataset(input_data, input_labels, input_subset, input_name=''):
	data, labels, fixation_s = get_subset(input_subset, input_data, input_labels)
	print('Subset contains ', len(labels), ' observations')

	compare_all_dimension_reduction(data, labels, True, input_name)
	#dimension_reduction(data, labels, 'PCA', True, 2)

	#train_data, test_data, train_labels, test_labels, sc = split_and_scale(data, labels)
	best_knn_confs = compare_knn_classifiers(data, labels, 'KNN Classification Performance Comparison' + '(' + input_name + ')')
	#best_knns = [fit_classifier(c, train_data, train_labels) for c in best_knn_confs]
	best_svm_confs = compare_svm_classifiers(data, labels, 'SVM Classification Performance Comparison' + '(' + input_name + ')')
	#best_svms = [fit_classifier(c, train_data, train_labels) for c in best_svm_confs]
	best_clf_confs = compare_classifiers(data, labels, 'Classification Performance Comparison' + '(' + input_name + ')')

	best_classifier_names = compare_input_sets(best_clf_confs, input_data, input_labels, 'Input Set Comparison' + '(' + input_name + ')')

#scaler = dh.get_scaler(dh.get_measured_spectra())

####################apply on unique estimated spectra dataset#######################
#current_data, current_labels, current_fixation, current_indexes = dh.get_scaled_subset('unique', dh.get_reconstructed_spectra(), dh.get_labels(), scaler)
#best_clf_confs_estimated = apply_on_dataset(estimated_spectra, positive_labels, 'unique', 'estimated')

acc, auc = cross_validate(get_classifier("SVM-linear-auto-True"), dh.get_reconstructed_spectra(), dh.get_labels(), 'unique', dh.get_measured_spectra(), None)
print(acc, auc)
