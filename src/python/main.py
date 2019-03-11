from os import mkdir, makedirs
from os.path import dirname, join as pjoin, exists
from sklearn import manifold
import sklearn.metrics
from sklearn.decomposition import PCA, FactorAnalysis, TruncatedSVD, FastICA 
from sklearn.preprocessing import StandardScaler
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis as QDA
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

get_indexes = lambda x, xs: [i for (y, i) in zip(xs, range(len(xs))) if x == y]


def create_directory(dirName):
	if not exists(dirName):
	    makedirs(dirName)

####################load matfiles#######################

base_dir  = '/media/sf_research/input/'
data_dir = 'saitama_v5_min_region'
out_dir = pjoin('/media/sf_research/output/', data_dir)
create_directory(out_dir)
matfile = pjoin(base_dir, data_dir, 'in.mat')
#mat_contents = sio.loadmat(matfile)
in_mat = loadmat(matfile)
print('Loaded mat file with headers', in_mat.keys())

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

###################Measured Data######

data, unique_ids = np.unique(spectra, return_index=True, axis=0)
measuredN = len(data)
labels = is_benign[unique_ids]
label_dict = {0: 'Malignant', 1: 'Benign'}
benign_ids = get_indexes(1, is_benign[unique_ids])
malignant_ids = get_indexes(0, is_benign[unique_ids])

train_data, test_data, train_labels, test_labels = train_test_split(data, labels, test_size=0.2, random_state=0)
sc = StandardScaler()
train_data = sc.fit_transform(train_data)
test_data = sc.transform(test_data)


def plot_lda(X, X_labels, title):
	plt.figure()
	positives_end = np.sum(X_labels==0) + 1
	negatives_end = np.sum(X_labels==1) + 1
	for label, marker, color, sample_start, sample_end in zip(range(2), ('s', 'o'), ('red', 'green'), (1, positives_end + 1), (positives_end, positives_end + negatives_end)):
		plt.scatter(x=np.array(range(sample_start, sample_end)),  
					y=X[:,0][X_labels == label], 
					marker=marker,
					color=color,
					alpha=0.7,
					label=label_dict[label])

	plt.xlabel('Samples')
	plt.ylabel('Linear Discriminant')
	leg = plt.legend(loc='upper right', fancybox=True)
	leg.get_frame().set_alpha(0.7)
	plt.title(title)
	plt.grid()
	plt.tight_layout
	plt.show()
	plt.savefig(pjoin(out_dir, 'dimension_reduction', title + '.png'), bbox_inches='tight')

def plot_da(X, X_labels, title):
	tag = title.strip('Analysis')

	plt.figure()
	for label, marker, color in zip(range(2), ('s', 'o'), ('red', 'green')):
		plt.scatter(x=X[:,0][X_labels == label], 
					y=X[:,1][X_labels == label] * -1, #To flip the plot
					marker=marker,
					color=color,
					alpha=0.7,
					label=label_dict[label])
	plt.xlabel(tag + '1')
	plt.ylabel(tag + '2')
	leg = plt.legend(loc='upper right', fancybox=True)
	leg.get_frame().set_alpha(0.7)
	plt.title(title)
	plt.grid()
	plt.tight_layout
	plt.show()
	plt.savefig(pjoin(out_dir, 'dimension_reduction', title + '.png'), bbox_inches='tight')

def plot_da3(X, X_labels, title):
	tag = title.strip('Analysis')

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	for label, marker, color in zip(range(2), ('s', 'o'), ('red', 'green')):
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

def dimension_reduction(data, data_labels, method, show_figures=False, components=2):

	create_directory(pjoin(out_dir, 'dimension_reduction'))

	if method == 'RandomForest':		
		##################RandomForest########
		rf = RandomForestRegressor(random_state=1, max_depth=10)
		rf.fit(data,data_labels)
		features = np.linspace(380, 780, 81)
		importances = rf.feature_importances_
		indices = np.argsort(importances)[-9:]  # top 10 features
		mean_importance = np.mean(importances)
		if show_figures:
			print('Average wavelength importance: ', mean_importance)
			plt.figure()
			plt.title('Wavelength Importances')
			plt.barh(range(len(indices)), importances[indices], color='b', align='center')
			plt.yticks(range(len(indices)), [features[i] for i in indices])
			plt.xlabel('Relative Importance')
			plt.show()

		method_obj = SelectFromModel(rf).fit(data, data_labels)
		data_transformed = method_obj.transform(data)
		print('Reduced dimensions: ', data_transformed.shape)
		return data_transformed, method_obj

	elif method == 'PCA': 
		##################PCA#################
		pca = PCA(n_components=components).fit(data)
		data_transformed = pca.transform(data)
		if show_figures:
			print("original dimension: ", data.shape)
			print("pca-projected dimension ", data_transformed.shape)
			print("pca explained variance ", pca.explained_variance_)
			plot_da(data_transformed, data_labels,'Principal Component Analysis')
			print('Reduced dimensions: ', data_transformed.shape)
		return data_transformed, pca

	elif method == 'SVD': 
		##################SVD########
		svd = TruncatedSVD(n_components=3, random_state=42).fit(data)
		data_transformed = svd.transform(data)
		if show_figures:
			plot_da3(data_transformed, data_labels, 'SVD Component Analysis') 
			plot_da_components(data_transformed, data_labels, 'SVD Components')
		return data_transformed, svd

	elif method == 'FactorAnalysis':
		##################FactorAnalysis########
		if show_figures:
			fa = FactorAnalysis(n_components=components).fit(data)
			data_transformed = fa.transform(data) #without labels 
			plot_da(data_transformed, data_labels, 'Factor Analysis') 
			print('Reduced dimensions: ', data_transformed.shape)

		fa = FactorAnalysis(n_components = 3).fit(data, data_labels)
		data_transformed = fa.transform(data) #with labels 
		if show_figures:
			plot_da3(data_transformed, data_labels, 'Factor Analysis') 
			print('Reduced dimensions: ', data_transformed.shape)
			plot_da_components(data_transformed, data_labels, 'Factor Analysis Components')
		return data_transformed, fa

	elif method == 'ICA':
		##################ICA########
		ica = FastICA(n_components=3, random_state=12, max_iter=500, tol=0.001).fit(data)
		data_transformed=ica.transform(data)
		if show_figures:
			plot_da(data_transformed, data_labels, 'Independent Component Analysis')
			plot_da3(data_transformed, data_labels, 'Independent Component Analysis') 
		return data_transformed, ica

	elif method == 'ISOMAP':
		##################ISOMAP########
		isomap = manifold.Isomap(n_neighbors=5, n_components=3, n_jobs=-1).fit(data)
		data_transformed = isomap.transform(data)
		if show_figures:
			plot_da_components(data_transformed, data_labels, 'ISOMAP Components')
			plot_da3(data_transformed, data_labels, 'ISOMAP Component Analysis') 
		return data_transformed, isomap

	elif method == 't-SNE':
		##################t-SNE########
		tsne = manifold.TSNE(n_components=3, n_iter=300).fit(data)
		data_transformed = tsne.transform(data)
		if show_figures:
			plot_da_components(data_transformed, data_labels, 't-SNE Components')
			plot_da3(data_transformed, data_labels, 't-SNE Component Analysis') 
		return data_transformed, tsne

	elif method == 'LDA':
		##################LDA#################
		lda = LDA(n_components=1, solver="svd", store_covariance=True).fit(data, data_labels)
		data_transformed = lda.transform(data)
		if show_figures:
			print('Reduced dimensions: ', data_transformed.shape)
			plot_lda(data_transformed, data_labels, 'Linear Discriminant Analysis')
		return data_transformed, lda

	elif method == 'QDA':
		##################QDA#################
		qda = QDA().fit(data, labels)
		data_transformed = qda.transform(data)
		#plot_da(data_transformed, 'Quadratic Discriminant Analysis')
		return data_transformed, qda

	else:
		print('Method not implemented.')
		return

def show_classifier_performance_stats(test_labels, pred_labels):
	cm = confusion_matrix(test_labels, pred_labels)
	accuracy = accuracy_score(test_labels, pred_labels)
	auc_score = roc_auc_score(test_labels, pred_labels)

	#print('Confusion matrix: ', cm)
	#print('Accuracy: ', accuracy)
	#print('AUC: ', auc_score)

	return accuracy, auc_score

def apply_dimension_reduction(train_data, test_data, train_labels, test_labels, dimred_method='', components=2):
	if dimred_method!='':
		train_dimred, dimred_obj = dimension_reduction(train_data, train_labels, dimred_method, False, components)
		test_dimred = dimred_obj.transform(test_data)
	else:
		train_dimred = train_data
		test_dimred = test_data

	return train_dimred, test_dimred

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
		train_data = dimension_reduction(train_data, train_labels, 'PCA', show_figures=False, components=10)
	clf.fit(train_data, train_labels)
	return clf

def cross_validate(classifier, data, labels, folds=5, method=''):
	cv = StratifiedKFold(n_splits=folds)

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
		mean_acc, mean_auc = cross_validate(clf, data, labels, 5, method)
		# train_dimred, test_dimred = apply_dimension_reduction(train_data, test_data, train_labels, test_labels, dimred_method=method, components=10)
		# clf.fit(train_dimred, train_labels)
		# score = clf.score(test_dimred, test_labels)
		# test_pred_labels = clf.predict(test_dimred) 
		# #print('Performance stats for ', name, ' ', method)
		#show_classifier_performance_stats(test_labels, test_pred_labels)
		performance_stats = np.append(performance_stats, [(mean_acc, mean_auc)], axis=0)
		performance_names.append(' '.join([name,method]))

	return performance_stats, performance_names

def compare_classifiers(data, labels, title=''):

	names = ["Nearest Neighbors", "KNN-1-correlation", "Linear SVM", "RBF SVM", "Decision Tree",
         "Random Forest", "Naive Bayes", "LDA", "QDA"]
	classifiers = [
	    KNN(3),
	    KNN(n_neighbors =1, algorithm = 'auto', metric='correlation'),
	    SVC(kernel="linear", C=0.025, probability=True),
	    SVC(gamma=2, C=1, probability=True),
	    DecisionTreeClassifier(max_depth=5),
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
	plot_accuracy_and_auc(performance_stats, performance_names, title)

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

	plot_accuracy_and_auc(performance_stats, performance_names, title)
	print('\n------------------------------------------------------')
	print('KNN performance \n')
	max_accuracy_index = (np.argsort(performance_stats[:,0]))[-10:]
	[print(-i, 'th best Accuracy ', x,  ' by ',y) for i,x,y in zip(range(-len(max_accuracy_index), 0), performance_stats[max_accuracy_index,0], [performance_names[i] for i in max_accuracy_index])]
	print(' ')
	max_auc_index = np.argsort(performance_stats[:,1])[-10:]
	[print(-i, 'th best A.U.C. ', x,  ' by ',y) for i,x,y in zip(range(-len(max_accuracy_index), 0), performance_stats[max_auc_index,1], [performance_names[i] for i in max_auc_index])]
	print('------------------------------------------------------')

	best_classifier_names = [performance_names[i] for i in  (np.argsort(performance_stats[:,0]))[-2:]] + [performance_names[i] for i in np.argsort(performance_stats[:,1])[-2:]]
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

	plot_accuracy_and_auc(performance_stats, performance_names, title)
	print('\n------------------------------------------------------')
	print('SVM performance \n')
	max_accuracy_index = (np.argsort(performance_stats[:,0]))[-10:]
	[print(-i, 'th best Accuracy ', x,  ' by ',y) for i,x,y in zip(range(-len(max_accuracy_index), 0), performance_stats[max_accuracy_index,0], [performance_names[i] for i in max_accuracy_index])]
	print(' ')
	max_auc_index = np.argsort(performance_stats[:,1])[-10:]
	[print(-i, 'th best A.U.C. ', x,  ' by ',y) for i,x,y in zip(range(-len(max_accuracy_index), 0), performance_stats[max_auc_index,1], [performance_names[i] for i in max_auc_index])]
	print('------------------------------------------------------')

	best_classifier_names = [performance_names[i] for i in  (np.argsort(performance_stats[:,0]))[-2:]] + [performance_names[i] for i in np.argsort(performance_stats[:,1])[-2:]]
	return best_classifier_names

best_knn_confs = compare_knn_classifiers(data, labels, 'KNN Classification Performance Comparison')
best_knns = [fit_classifier(c, train_data, train_labels) for c in best_knn_confs]
print(best_knn_confs)
best_svm_confs = compare_svm_classifiers(data, labels, 'SVM Classification Performance Comparison')
best_svms = [fit_classifier(c, train_data, train_labels) for c in best_svm_confs]
print(best_svm_confs)

compare_classifiers(data, labels, 'Classification Performance Comparison')


# # By hardic goel
# def estimate_lda_params(data):
# 	grouped = data.groupby(labels)
# 	means = ()
# 	for c in range(2): 
# 		means[c] = np.array(data[:][labels == c]).mean(axis = 0)

# 	overall_mean = np.array(data).mean(axis = 0)

# 	#between class variance 
# 	S_b = np.zeros((data.shape[1]-1, data.shape[1]-1))
# 	for c in means.keys():
# 		S_b += np.multiply(len(data[:][labels == c]), 
# 							np.outer((means[c] - overall_mean), (means[c]-overall_mean)))

# 	#within class covariance 
# 	S_w = np.zeros(S_b.shape)
# 	for c in range(2):
# 		tmp = np.subtract((data[:][labels == c]).T, np.expand_dims(means[c], axis = 1))
# 		S_w = np.add(np.dot(tmp, tmp.T), S_w)

# 	matrix = np.dot(np.linalg.pinv(S_w), S_b)
# 	eigenvals, eigenvecs = np.linalg.eig(matrix)
# 	eiglist = [(eigvals[i], eigenvecs[:,i]) for i in range(len(eigvals))]

# 	eiglist = sorted(eiglist, key = lambda x : x[0], reverse = True)

# 	w = np.array([eiglist[i][1] for i in range(2)])

# 	tot = 0
# 	for c in means.keys():
# 		tot += np.dot(w, means[c])
# 	w0 = 0.5 * tot
	
# 	c1 = means.keys()[0]
# 	c2 = means.keys()[1]
# 	mu1 = np.dot(w, means[c1])
# 	if (mu1 >= w0):
# 		c1 = 0
# 	else: 
# 		c1 = 1

# 	return w, w0, c1, means

# def calculate_lda_score(inputs, w, w0, c1, means):
# 	proj = np.dot(w, inputs.T).T
# 	c1 = means.keys()[0]
# 	c2 = means.keys()[1]
# 	if (c1 == 1): 
# 		proj = [c1 if proj[i] >= w0 else c2 for i in range(len(proj))]
# 	else:
# 		proj = [c1 if proj[i] < w0 else c2 for i in range(len(proj))]

# 	errors = (proj != labels)
# 	return sum(errors)


# w, w0, c1, means = estimate_lda_params(data)
