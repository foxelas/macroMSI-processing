from os import mkdir, makedirs
from os.path import dirname, join as pjoin, exists
from sklearn import manifold
from sklearn.decomposition import PCA, FactorAnalysis, TruncatedSVD, FastICA 
from sklearn.preprocessing import StandardScaler
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis as QDA
from sklearn.feature_selection import SelectFromModel
from sklearn.ensemble import RandomForestRegressor
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import scipy.io as sio 
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

measured_spectra, unique_ids = np.unique(spectra, return_index=True, axis=0)
sc = StandardScaler()
measured_spectra = sc.fit_transform(measured_spectra)
measuredN = len(measured_spectra)
labels = is_benign[unique_ids]
label_dict = {0: 'Malignant', 1: 'Benign'}
benign_ids = get_indexes(1, is_benign[unique_ids])
malignant_ids = get_indexes(0, is_benign[unique_ids])

def plot_lda(X, title):
	plt.figure()
	positives_end = np.sum(labels==0) + 1
	negatives_end = np.sum(labels==1) + 1
	for label, marker, color, sample_start, sample_end in zip(range(2), ('s', 'o'), ('red', 'green'), (1, positives_end + 1), (positives_end, positives_end + negatives_end)):
		print(np.array(range(sample_start, sample_end)).shape)
		print(len(X[:,0][labels == label]))
		plt.scatter(x=np.array(range(sample_start, sample_end)),  
					y=X[:,0][labels == label], 
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

def plot_da(X, title):
	tag = title.strip('Analysis')

	plt.figure()
	for label, marker, color in zip(range(2), ('s', 'o'), ('red', 'green')):
		plt.scatter(x=X[:,0][labels == label], 
					y=X[:,1][labels == label] * -1, #To flip the plot
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

def plot_da3(X, title):
	tag = title.strip('Analysis')

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	for label, marker, color in zip(range(2), ('s', 'o'), ('red', 'green')):
		ax.scatter(xs=X[:,0][labels == label], 
					ys=X[:,1][labels == label] * -1, #To flip the plot
					zs=X[:,2][labels == label] * -1, #To flip the plot
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

def plot_da_components(X, title):
	plt.figure()
	plt.title(title)
	plt.scatter(X[:,0], X[:,1], label='1st and 2nd Component')
	plt.scatter(X[:,1], X[:,2], label='2nd and 3rd Component')
	plt.scatter(X[:,2], X[:,0], label='1st and 3rd Component')
	leg = plt.legend(loc='upper right', fancybox=True)
	leg.get_frame().set_alpha(0.7)
	plt.show()
	plt.savefig(pjoin(out_dir, 'dimension_reduction', title + '.png'), bbox_inches='tight')

def dimension_reduction(method):
	if method == 'RandomForest':		
		##################RandomForest########
		rf = RandomForestRegressor(random_state=1, max_depth=10)
		rf.fit(measured_spectra,labels)
		features = np.linspace(380, 780, 81)
		importances = rf.feature_importances_
		indices = np.argsort(importances)[-9:]  # top 10 features
		mean_importance = np.mean(importances)
		print('Average wavelength importance: ', mean_importance)
		plt.figure()
		plt.title('Wavelength Importances')
		plt.barh(range(len(indices)), importances[indices], color='b', align='center')
		plt.yticks(range(len(indices)), [features[i] for i in indices])
		plt.xlabel('Relative Importance')
		plt.show()
		measured_rf = SelectFromModel(rf).fit_transform(measured_spectra, labels)
		print('Reduced dimensions: ', measured_rf.shape)
		return measured_rf

	elif method == 'PCA': 
		##################PCA#################
		pca = PCA(n_components=2)
		pca.fit(measured_spectra)
		measured_pca = pca.transform(measured_spectra)
		print("original dimension: ", measured_spectra.shape)
		print("pca-projected dimension ", measured_pca.shape)
		print("pca explained variance ", pca.explained_variance_)
		plot_da(measured_pca, 'Principal Component Analysis')
		print('Reduced dimensions: ', measured_pca.shape)
		return measured_pca

	elif method == 'SVD': 
		##################SVD########
		measured_svd = TruncatedSVD(n_components=3, random_state=42).fit_transform(measured_spectra)
		plot_da3(measured_svd, 'SVD Component Analysis') 
		plot_da_components(measured_svd, 'SVD Components')
		return measured_svd

	elif method == 'FactorAnalysis':
		##################FactorAnalysis########
		measured_fa = FactorAnalysis(n_components = 2).fit_transform(measured_spectra) #without labels 
		plot_da(measured_fa, 'Factor Analysis') 
		print('Reduced dimensions: ', measured_fa.shape)

		measured_fa = FactorAnalysis(n_components = 3).fit_transform(measured_spectra, labels) #with labels 
		plot_da3(measured_fa, 'Factor Analysis') 
		print('Reduced dimensions: ', measured_fa.shape)
		plot_da_components(measured_fa, 'Factor Analysis Components')
		return measured_fa

	elif method == 'ICA':
		##################ICA########
		ica = FastICA(n_components=3, random_state=12, max_iter=500, tol=0.001) 
		measured_ica=ica.fit_transform(measured_spectra)
		plot_da(measured_ica, 'Independent Component Analysis')
		plot_da3(measured_ica, 'Independent Component Analysis') 
		return measured_ica

	elif method == 'ISOMAP':
		##################ISOMAP########
		measured_isomap = manifold.Isomap(n_neighbors=5, n_components=3, n_jobs=-1).fit_transform(measured_spectra)
		plot_da_components(measured_isomap, 'ISOMAP Components')
		plot_da3(measured_isomap, 'ISOMAP Component Analysis') 
		return measured_isomap

	elif method == 't-SNE':
		##################t-SNE########
		measured_tsne = manifold.TSNE(n_components=3, n_iter=300).fit_transform(measured_spectra)
		plot_da_components(measured_tsne, 't-SNE Components')
		plot_da3(measured_tsne, 't-SNE Component Analysis') 
		return measured_tsne

	elif method == 'LDA':
		##################LDA#################
		lda = LDA(n_components=2, solver="svd", store_covariance=True)
		#y_pred = lda.fit(measured_spectra, labels).predict(measured_spectra)
		measured_lda = lda.fit(measured_spectra, labels).transform(measured_spectra)
		print('Reduced dimensions: ', measured_lda.shape)
		plot_lda(measured_lda, 'Linear Discriminant Analysis')
		return measured_lda

	elif method == 'QDA':
		##################QDA#################
		qda = QDA()
		measured_qda = qda.fit(measured_spectra, labels).predict(measured_spectra)
		#plot_da(measured_qda, 'Quadratic Discriminant Analysis')
		return measured_qda

	else:
		print('Method not implemented.')
		return

create_directory(pjoin(out_dir, 'dimension_reduction'))

measured_reduced = dimension_reduction('PCA')
measured_reduced = dimension_reduction('LDA')

#from sklearn.model_selection import train_test_split
#X_train, X_test, y_train, y_test = train_test_split(X,y, test_size=0.5)

# By hardic goel
def estimate_lda_params(data):
	grouped = data.groupby(labels)
	means = ()
	for c in range(2): 
		means[c] = np.array(data[:][labels == c]).mean(axis = 0)

	overall_mean = np.array(data).mean(axis = 0)

	#between class variance 
	S_b = np.zeros((data.shape[1]-1, data.shape[1]-1))
	for c in means.keys():
		S_b += np.multiply(len(data[:][labels == c]), 
							np.outer((means[c] - overall_mean), (means[c]-overall_mean)))

	#within class covariance 
	S_w = np.zeros(S_b.shape)
	for c in range(2):
		tmp = np.subtract((data[:][labels == c]).T, np.expand_dims(means[c], axis = 1))
		S_w = np.add(np.dot(tmp, tmp.T), S_w)

	matrix = np.dot(np.linalg.pinv(S_w), S_b)
	eigenvals, eigenvecs = np.linalg.eig(matrix)
	eiglist = [(eigvals[i], eigenvecs[;,i]) for i in range(len(eigvals))]_

	eiglist = sorted(eiglist, key = lambda x : x[0], reverse = True)

	w = np.array([eiglist[i][1] for i in range(2)])

	tot = 0
	for c in means.keys():
		tot += np.dot(w, means[c])
	w0 = 0.5 * tot
	
	c1 = means.keys()[0]
	c2 = means.keys()[1]
	mu1 = np.dot(w, means[c1])
	if (mu1 >= w0):
		c1 = 0
	else: 
		c1 = 1

	return w, w0, c1, means

def calculate_lda_score(inputs, w, w0, c1, means):
	proj = np.dot(w, inputs.T).T
	c1 = means.keys()[0]
	c2 = means.keys()[1]
	if (c1 == 1): 
		proj = [c1 if proj[i] >= w0 else c2 for i in range(len(proj))]
	else:
		proj = [c1 if proj[i] < w0 else c2 for i in range(len(proj))]

	errors = (proj != labels)
	return sum(errors)


w, w0, c1, means = estimate_lda_params(measured_spectra)
