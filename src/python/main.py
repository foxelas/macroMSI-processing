from os import mkdir, makedirs
from os.path import dirname, join as pjoin, exists
from sklearn import manifold
import sklearn.metrics
from sklearn.decomposition import PCA, FactorAnalysis, TruncatedSVD, FastICA 
from sklearn.preprocessing import StandardScaler
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA, QuadraticDiscriminantAnalysis as QDA
from sklearn.neighbors import KNeighborsClassifier as KNN
from sklearn.feature_selection import SelectFromModel
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.metrics import confusion_matrix, accuracy_score, auc, roc_auc_score, roc_curve
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy import interp
import scipy.io as sio 
from scipy.spatial.distance import chebyshev, correlation, hamming, minkowski
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import interactive, is_interactive
import time 

interactive(True)
# exec(open("main.py").read())

def loadmat(filename):
    '''
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    
    from: `StackOverflow <http://stackoverflow.com/questions/7008608/scipy-io-loadmat-nested-structures-i-e-dictionaries>`_
    '''
    data = sio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)

def _check_keys(dict):
    '''
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    '''
    for key in dict:
        if isinstance(dict[key], sio.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict        

def _todict(matobj):
    '''
    A recursive function which constructs from matobjects nested dictionaries
    '''
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, sio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict

get_indexes_equal = lambda x, xs: [i for (y, i) in zip(xs, range(len(xs))) if x == y]
get_indexes = lambda x, xs: [i for (y, i) in zip(xs, range(len(xs))) if y in  x]

def create_directory(dirName):
	if not exists(dirName):
	    makedirs(dirName)

####################load matfiles#######################

base_dir  = '/media/sf_research/input/'
data_dir = 'saitama_v5_min_region'
out_dir = pjoin('/media/sf_research/output/', data_dir)
create_directory(out_dir)
matfile = pjoin(base_dir, data_dir, 'in.mat')
in_mat = loadmat(matfile)
print('Loaded mat file with headers', in_mat.keys())

####################prepare data structs#######################
complete_spectra = in_mat['CompleteSpectra']
darkIs = in_mat['DarkIs']
msi_names = in_mat['MSINames']
msis = in_mat['MSIs']
masks = in_mat['Masks']
spectra = in_mat['Spectra']
spectra_names = in_mat['SpectraNames']
whiteIs = in_mat['WhiteIs']

msiN, lambdaN = complete_spectra.shape
msi_dim = msis[0].shape
bandN = msi_dim[0]
rgbN = msi_dim[3]

matfile = pjoin(base_dir, data_dir, 'ID.mat')
id_mat = loadmat(matfile)
print('Loaded mat file with headers', id_mat.keys())
id_struct = id_mat['ID']
is_benign = np.array([id_struct[x].IsBenign for x in range(msiN) ])
is_cut = np.array([id_struct[x].IsCut for x in range(msiN) ])
is_fixed = np.array([id_struct[x].IsFixed for x in range(msiN) ])
positive_labels = np.array([int((x + 1) % 2) for x in is_benign])
label_dict = {0: 'Benign', 1: 'Malignant'}

fixation = []
for i in range(msiN):
	if is_fixed[i] and not(is_cut[i]):
		fixation.append('fixed')	
	elif is_cut[i]:
		fixation.append('cut')
	else:
		fixation.append('unfixed')

def subset_indexes(name, data, fixation_type):
	subset_ids = []
	if name == 'unique':
		data_s, subset_ids = np.unique(data, return_index=True, axis=0)

	elif name == 'unfixed' or name == 'fixed' or name == 'cut':
		subset_ids = np.array(get_indexes_equal(name, fixation_type))

	elif 'unique' in name:
		name1, name2 = name.split('_')
		subset_ids1 = subset_indexes(name1, data, fixation_type)
		subset_ids2 = subset_indexes(name2, data, fixation_type)
		subset_ids = np.intersect1d(subset_ids1, subset_ids2)

	else: 
		print('Not implemented yet.')
		return 

	return subset_ids

def get_subset(name, data, labels, fixation_type):
	subset_ids = subset_indexes(name, data, fixation_type)
	data_s = data[subset_ids,:]
	labels_s = labels[subset_ids]
	fixation_type_s = [fixation_type[x] for x in subset_ids]

	return data_s, labels_s, fixation_type_s

####################get unique measured spectra dataset#######################
data, labels, fixation_s = get_subset('unique', spectra, positive_labels, fixation)
measuredN = len(data)
print('Subset contains ', measuredN, ' observations')


####################dimension reduction#######################
sc = StandardScaler()

def plot_da(X, X_labels, title):
	fig = plt.figure()
	tag = title.strip('Analysis')
	n_components = X.shape[1]

	if n_components == 1 : #When X is reduced to 1 dimension
		positives_end = np.sum(X_labels==0) + 1
		negatives_end = np.sum(X_labels==1) + 1
		for label, marker, color, sample_start, sample_end in zip(range(2), ('s', 'o'), ('green', 'red'), (1, positives_end + 1), (positives_end, positives_end + negatives_end)):
			plt.scatter(x=np.array(range(sample_start, sample_end)),  
						y=X[:,0][X_labels == label], 
						marker=marker,
						color=color,
						alpha=0.7,
						label=label_dict[label])
		plt.xlabel('Samples')
		plt.ylabel(tag)

	elif n_components == 2 : #when X is reduced to 2 dimensions
		for label, marker, color in zip(range(2), ('s', 'o'), ('green', 'red')):
			plt.scatter(x=X[:,0][X_labels == label], 
						y=X[:,1][X_labels == label] * -1, #To flip the plot
						marker=marker,
						color=color,
						alpha=0.7,
						label=label_dict[label])
		plt.xlabel(tag + '1')
		plt.ylabel(tag + '2')

	elif n_components == 3 : #when X is reduced to 3 dimensions
		ax = fig.add_subplot(111, projection='3d')
		for label, marker, color in zip(range(2), ('s', 'o'), ('green', 'red')):
			ax.scatter(xs=X[:,0][X_labels == label], 
						ys=X[:,1][X_labels == label] * -1, #To flip the plot
						zs=X[:,2][X_labels == label] * -1, #To flip the plot
						marker=marker,
						color=color,
						alpha=0.7,
						label=label_dict[label])
		ax.set_xlabel(tag + '1')
		ax.set_ylabel(tag + '2')
		ax.set_zlabel(tag + '3')

	else:
		print('Unsupported plot type for input with greater than 3 dimensions.')
		return 

	leg = plt.legend(loc='upper right', fancybox=True)
	leg.get_frame().set_alpha(0.7)
	plt.title(title)
	plt.grid()
	plt.tight_layout
	plt.show()
	plt.savefig(pjoin(out_dir, 'dimension_reduction', title + '.png'), bbox_inches='tight')


def plot_da_components(X, X_labels, title):
	plt.figure()
	plt.title(title)
	plt.scatter(X[:,0], X[:,1], label='1st and 2nd Component')
	plt.scatter(X[:,1], X[:,2], label='2nd and 3rd Component')
	plt.scatter(X[:,2], X[:,0], label='1st and 3rd Component')
	leg = plt.legend(loc='upper right', fancybox=True)
	leg.get_frame().set_alpha(0.7)
	plt.show()
	plt.savefig(pjoin(out_dir, 'dimension_reduction', title + '.png'), bbox_inches='tight')

def plot_feature_importance(features, importances, indices, title):
	plt.figure()
	plt.title(title)
	plt.barh(range(len(indices)), importances[indices], color='b', align='center')
	plt.yticks(range(len(indices)), [features[i] for i in indices])
	plt.xlabel('Relative Importance')
	plt.show()
	plt.savefig(pjoin(out_dir, 'dimension_reduction', title + '.png'), bbox_inches='tight')

def dimension_reduction(data, data_labels, method, show_figures=False, components=2):
	create_directory(pjoin(out_dir, 'dimension_reduction'))

	if method == 'RandomForest':		
		##################RandomForest########
		rf = RandomForestRegressor(random_state=1, max_depth=10, n_estimators=10)
		rf.fit(data,data_labels)
		features = np.linspace(380, 780, 81)
		importances = rf.feature_importances_
		indices = np.argsort(importances)[-9:]  # top 10 features
		mean_importance = np.mean(importances)
		print('Random Forest-Average wavelength importance: ', mean_importance)
		if show_figures:
			plot_feature_importance(features, importances, indices, 'Wavelength Importances by Random Forest')

		dimred = SelectFromModel(rf).fit(data, data_labels)
		data_transformed = dimred.transform(data)

	elif method == 'PCA': 
		##################PCA#################
		dimred = PCA(n_components=components).fit(data)
		data_transformed = dimred.transform(data)
		if show_figures:
			print("PCA-explained variance ", dimred.explained_variance_ratio_)
			plot_da(data_transformed, data_labels,'Principal Component Analysis')

	elif method == 'SVD': 
		##################SVD########
		dimred = TruncatedSVD(n_components=components, random_state=42).fit(data)
		data_transformed = dimred.transform(data)
		if show_figures:
			plot_da(data_transformed, data_labels, 'SVD Component Analysis') 
			plot_da_components(data_transformed, data_labels, 'SVD Components')

	elif method == 'FactorAnalysis':
		##################FactorAnalysis########
		dimred = FactorAnalysis(n_components = components).fit(data, data_labels)
		data_transformed = dimred.transform(data) #with labels 
		if show_figures:
			plot_da(data_transformed, data_labels, 'Factor Analysis') 
			plot_da_components(data_transformed, data_labels, 'Factor Analysis Components')

	elif method == 'ICA':
		##################ICA########
		dimred = FastICA(n_components=components, random_state=1, max_iter=500, tol=0.01).fit(data)
		data_transformed=dimred.transform(data)
		if show_figures:
			plot_da(data_transformed, data_labels, 'Independent Component Analysis')

	elif method == 'ISOMAP':
		##################ISOMAP########
		dimred = manifold.Isomap(n_neighbors=5, n_components=components, n_jobs=-1).fit(data)
		data_transformed = dimred.transform(data)
		if show_figures:
			plot_da_components(data_transformed, data_labels, 'ISOMAP Components')
			plot_da(data_transformed, data_labels, 'ISOMAP Component Analysis') 

	elif method == 't-SNE':
		##################t-SNE########
		dimred = manifold.TSNE(n_components=components, n_iter=300).fit(data)
		data_transformed = dimred.fit_transform(data)
		if show_figures:
			plot_da_components(data_transformed, data_labels, 't-SNE Components')
			plot_da(data_transformed, data_labels, 't-SNE Component Analysis') 

	elif method == 'LDA':
		##################LDA#################
		dimred = LDA(n_components=components, solver="svd", store_covariance=True).fit(data, data_labels)
		data_transformed = dimred.transform(data)
		if show_figures:
			plot_da(data_transformed, data_labels, 'Linear Discriminant Analysis')

	elif method == 'QDA':
		##################QDA#################
		dimred = QDA().fit(data, labels)
		data_transformed = np.array([[x] for x in dimred.decision_function(data)])
		plot_da(data_transformed, data_labels, 'Quadratic Discriminant Analysis')

	else:
		print('Method not implemented.')
		return

	if show_figures:
		print(method + '-Reduced dimensions: ', data_transformed.shape)

	return data_transformed, dimred


def apply_dimension_reduction(train_data, test_data, train_labels, test_labels, dimred_method='', components=2):
	if dimred_method!='':
		train_dimred, dimred_obj = dimension_reduction(train_data, train_labels, dimred_method, False, components)
		test_dimred = dimred_obj.transform(test_data)
	else:
		train_dimred = train_data
		test_dimred = test_data

	return train_dimred, test_dimred

def compare_all_dimension_reduction(data, data_labels, show_figures=False):
	dimension_reduction(data, labels, 'RandomForest', show_figures, 2)
	dimension_reduction(data, labels, 'PCA', show_figures, 2)
	dimension_reduction(data, labels, 'SVD', show_figures, 3)
	dimension_reduction(data, labels, 'FactorAnalysis', show_figures, 3)
	dimension_reduction(data, labels, 'ICA', show_figures, 3)
	dimension_reduction(data, labels, 'ISOMAP', show_figures, 3)
	dimension_reduction(data, labels, 't-SNE', show_figures, 3)
	dimension_reduction(data, labels, 'LDA', show_figures, 1)
	dimension_reduction(data, labels, 'QDA', show_figures, 1)

compare_all_dimension_reduction(data, labels, True)
#dimension_reduction(data, labels, 'PCA', True, 2)

####################classification#######################

def show_classifier_performance_stats(test_labels, pred_labels):
	cm = confusion_matrix(test_labels, pred_labels)
	accuracy = accuracy_score(test_labels, pred_labels)
	auc_score = roc_auc_score(test_labels, pred_labels)

	#print('Confusion matrix: ', cm)
	#print('Accuracy: ', accuracy)
	#print('AUC: ', auc_score)

	return accuracy, auc_score

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

def fit_classifier(classifier_name, train_data, train_labels):
	if 'SVM' in classifier_name:
		a, b, c, d = classifier_name.split('-') 
		d = d if not('PCA' in d) else (d.split(' '))[0] 
		c = c if (c == 'auto' or c == 'scale') else float(c)
		clf = SVC(kernel=b, gamma=c, shrinking=(d==True))
	elif 'KNN' in classifier_name:
		a, b, c =  classifier_name.split('-')
		c = c.strip(' ') if not('PCA' in c) else (c.split(' '))[0] 
		clf = KNN(n_neighbors=int(b), metric=c)
	else:
		print('Not yet implemented.')
		return

	if 'PCA' in classifier_name: 
		train_data, dimred_obj = dimension_reduction(train_data, train_labels, 'PCA', show_figures=False, components=10)
	clf.fit(train_data, train_labels)
	return clf

def cross_validate(classifier, data, labels, folds=5, method=''):
	cv = StratifiedKFold(n_splits=folds, shuffle=True, random_state=1)

	tprs = []
	aucs = []
	accs = []
	fpr_f = []
	tpr_f = []
	mean_fpr = np.linspace(0, 1, 100)

	for train, test in cv.split(data, labels):

		train_f = sc.fit_transform(data[train])
		test_f = sc.transform(data[test])

		train_f, test_f = apply_dimension_reduction(train_f, test_f, labels[train], labels[test], method, 10)
		classifier.fit(train_f, labels[train])
		pred_labels = classifier.predict(test_f)
		accs.append(accuracy_score(labels[test], pred_labels))

		pred_probas = classifier.predict_proba(test_f)
		# Compute ROC curve and area the curve
		fpr, tpr, thresholds = roc_curve(labels[test], pred_probas[:, 1])
		fpr_f.append(fpr)
		tpr_f.append(tpr)
		tprs.append(interp(mean_fpr, fpr, tpr))
		tprs[-1][0] = 0.0
		roc_auc = auc(fpr, tpr)
		aucs.append(roc_auc)

	mean_tpr = np.mean(tprs, axis=0)
	mean_tpr[-1] = 1.0
	mean_auc = auc(mean_fpr, mean_tpr)
	std_auc = np.std(aucs)

	mean_acc = np.mean(accs)

	return mean_acc, mean_auc 

def get_best_classifiers(performance_stats, performance_names, title, best_n=10):
	plot_accuracy_and_auc(performance_stats, performance_names, title)
	print('\n------------------------------------------------------')
	print(performance_names, ': Best performance \n')
	max_accuracy_index = (np.argsort(performance_stats[:,0]))[-best_n:]
	[print(-i, 'th best Accuracy ', x,  ' by ',y) for i,x,y in zip(range(-len(max_accuracy_index), 0), performance_stats[max_accuracy_index,0], [performance_names[i] for i in max_accuracy_index])]
	print(' ')
	max_auc_index = np.argsort(performance_stats[:,1])[-best_n:]
	[print(-i, 'th best A.U.C. ', x,  ' by ',y) for i,x,y in zip(range(-len(max_accuracy_index), 0), performance_stats[max_auc_index,1], [performance_names[i] for i in max_auc_index])]
	print('------------------------------------------------------')

	best_classifier_names = [performance_names[i] for i in  (np.argsort(performance_stats[:,0]))[-best_n:]] + [performance_names[i] for i in np.argsort(performance_stats[:,1])[-best_n:]]
	return best_classifier_names

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


def classify_many(names, classifiers, data, labels, method=''):
	performance_stats = np.empty((0, 2))
	performance_names = []

	#original spectum as input
	for name, clf in zip(names, classifiers):
		mean_acc, mean_auc = cross_validate(clf, data, labels, 10, method)
		performance_stats = np.append(performance_stats, [(mean_acc, mean_auc)], axis=0)
		performance_names.append(' '.join([name,method]))

	return performance_stats, performance_names

def compare_classifiers(data, labels, title=''):

	names = ["KNN-3-minkowski", "KNN-1-correlation", "SVM-linear-0.025", "SVM-linear-auto-True", "SVM-rbf-scale-False", "Decision Tree-gini",
         "Decision Tree-entropy", "Random Forest", "Naive Bayes", "LDA", "QDA"]
	classifiers = [
	    KNN(n_neighbors =5, algorithm = 'auto', metric='cityblock'),
	    KNN(n_neighbors =1, algorithm = 'auto', metric='correlation'),
	    SVC(kernel="linear", C=0.025, probability=True),
	    SVC(kernel="linear", gamma="auto", shrinking=True, probability=True),
	    SVC(kernel="rbf", gamma="scale", shrinking=False, probability=True),
	    DecisionTreeClassifier(criterion="gini"),
	    DecisionTreeClassifier(criterion="entropy"),
	    RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1),
	    GaussianNB(),
	    LDA(),
	    QDA()]

#original spectum as input
	pf_stats1, pf_names1 = classify_many(names, classifiers, data, labels)
#10-pc PCA-projected data as input
	pf_stats2, pf_names2 = classify_many(names, classifiers, data, labels, 'PCA')

	performance_stats = np.concatenate((pf_stats1, pf_stats2), axis=0)
	performance_names = pf_names1 + pf_names2

	best_classifier_names = get_best_classifiers(performance_stats, performance_names, title, 10)
	return best_classifier_names

def compare_knn_classifiers(data, labels, title=''):
	performance_stats = np.empty((0, 2))
	performance_names = []
	names = ["KNN-1-minkowski","KNN-3-minkowski","KNN-5-minkowski",
		"KNN-1-cityblock","KNN-3-cityblock","KNN-5-cityblock",
		"KNN-1-euclidean","KNN-3-euclidean","KNN-5-euclidean",
		"KNN-1-chebyshev","KNN-3-chebyshev","KNN-5-chebyshev",
		"KNN-1-hamming","KNN-3-hamming","KNN-5-hamming",
		"KNN-1-correlation","KNN-3-correlation","KNN-5-correlation"]
	classifiers = [
	    KNN(n_neighbors =1, algorithm = 'auto', metric='minkowski'),
	    KNN(n_neighbors =3, algorithm = 'auto', metric='minkowski'),
	    KNN(n_neighbors =5, algorithm = 'auto', metric='minkowski'),
	    KNN(n_neighbors =1, algorithm = 'auto', metric='cityblock'),
	    KNN(n_neighbors =3, algorithm = 'auto', metric='cityblock'),
	    KNN(n_neighbors =5, algorithm = 'auto', metric='cityblock'),
	    KNN(n_neighbors =1, algorithm = 'auto', metric='euclidean'),
	    KNN(n_neighbors =3, algorithm = 'auto', metric='euclidean'),
	    KNN(n_neighbors =5, algorithm = 'auto', metric='euclidean'),
	    KNN(n_neighbors =1, algorithm = 'auto', metric='chebyshev'),
	    KNN(n_neighbors =3, algorithm = 'auto', metric='chebyshev'),
	    KNN(n_neighbors =5, algorithm = 'auto', metric='chebyshev'),
	    KNN(n_neighbors =1, algorithm = 'auto', metric='hamming'),
	    KNN(n_neighbors =3, algorithm = 'auto', metric='hamming'),
	    KNN(n_neighbors =5, algorithm = 'auto', metric='hamming'),
	    KNN(n_neighbors =1, algorithm = 'auto', metric='correlation'),
	    KNN(n_neighbors =3, algorithm = 'auto', metric='correlation'),
	    KNN(n_neighbors =5, algorithm = 'auto', metric='correlation')
	]
#original spectum as input
	pf_stats1, pf_names1 = classify_many(names, classifiers, data, labels)
#10-pc PCA-projected data as input
	pf_stats2, pf_names2 = classify_many(names, classifiers, data, labels, 'PCA')
	
	performance_stats = np.concatenate((pf_stats1, pf_stats2), axis=0)
	performance_names = pf_names1 + pf_names2

	performance_stats, unique_ids = np.unique(performance_stats, return_index=True, axis=0)
	performance_names = [performance_names[i] for i in unique_ids]

	best_classifier_names = get_best_classifiers(performance_stats, performance_names, title, 2)
	return best_classifier_names

def compare_svm_classifiers(data, labels, title=''):
	performance_stats = np.empty((0, 2))
	performance_names = []
	names = []
	classifiers = []

	for kernel in ('rbf', 'poly', 'linear', 'sigmoid'):
		for gamma in {'auto', 'scale', 0.5, 1, 2}:
				for shrinking in (True, False):
					if not(kernel == 'linear' and gamma != 'auto'):
						names.append('-'.join(['SVM', kernel, str(gamma), str(shrinking)]))
						classifiers.append(SVC(kernel=kernel, gamma=gamma, shrinking=shrinking, probability=True))

	performance_stats, performance_names = classify_many(names, classifiers, data, labels)

	performance_stats, unique_ids = np.unique(performance_stats, return_index=True, axis=0)
	performance_names = [performance_names[i] for i in unique_ids]

	best_classifier_names = get_best_classifiers(performance_stats, performance_names, title, 2)
	return best_classifier_names

def get_classifier(classifier_name):
	classifier_name = classifier_name.replace('PCA', '').strip(' ')
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

def compare_input_sets(names, data, labels, fixation_type):
	classifiers = [get_classifier(c) for c in names]

	performance_stats = np.empty((0, 2))
	performance_names = []
	for input_set in ('unique_unfixed', 'unique_fixed', 'unique_cut', 'unique'):
		data_s, labels_s, fixation_s = get_subset(input_set, data, labels, fixation)
		pf_stats1, pf_names1 = classify_many([ '-'.join([x, input_set]) for x in names], classifiers, data_s, labels_s)

		performance_stats = np.concatenate((performance_stats, pf_stats1), axis=0)
		performance_names = performance_names + pf_names1

	best_performance_names = [performance_names[i] for i in  (np.argsort(performance_stats[:,0]))[-15:]] + [performance_names[i] for i in np.argsort(performance_stats[:,1])[-15:]]
	best_performance_stats = np.concatenate(([performance_stats[i,:] for i in  (np.argsort(performance_stats[:,0]))[-15:]], [performance_stats[i,:] for i in np.argsort(performance_stats[:,1])[-15:]]), axis=0)
	best_performance_stats, unique_ids = np.unique(best_performance_stats, return_index=True, axis=0)
	best_performance_names = [best_performance_names[i] for i in unique_ids]
	best_performance_ids = np.argsort([x[::-1] for x in best_performance_names])
	best_performance_names = [best_performance_names[x] for x in best_performance_ids]
	best_performance_stats = best_performance_stats[best_performance_ids,:]
	print(best_performance_names)
	print(best_performance_stats)

	plot_accuracy_and_auc(best_performance_stats, best_performance_names, 'Compare Input Sets')

train_data, test_data, train_labels, test_labels = train_test_split(data, labels, test_size=0.25, random_state=1)
train_data = sc.fit_transform(train_data)
test_data = sc.transform(test_data)

best_knn_confs = compare_knn_classifiers(data, labels, 'KNN Classification Performance Comparison')
best_knns = [fit_classifier(c, train_data, train_labels) for c in best_knn_confs]
print('Best KNN classifier configurations:\n', best_knn_confs)
best_svm_confs = compare_svm_classifiers(data, labels, 'SVM Classification Performance Comparison')
best_svms = [fit_classifier(c, train_data, train_labels) for c in best_svm_confs]
print('Best SVM classifier configurations:\n', best_svm_confs)

best_clf_confs = compare_classifiers(data, labels, 'Classification Performance Comparison')
print('Best overall classifier configurations:\n', best_clf_confs)

best_classifier_names = compare_input_sets(best_clf_confs, data, labels, fixation)
print('Best classsifiers based on input set:\n', best_classifier_names)