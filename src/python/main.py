from os.path import dirname, join as pjoin
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
from matplotlib import interactive 
import time 

#interactive(True)

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

####################load matfiles#######################

base_dir  = '/media/sf_research/input/'
data_dir = 'saitama_v5_min_region'
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

def plot_da_components(X, title):
	plt.figure()
	plt.title(title)
	plt.scatter(X[:,0], X[:,1], label='1st and 2nd Component')
	plt.scatter(X[:,1], X[:,2], label='2nd and 3rd Component')
	plt.scatter(X[:,2], X[:,0], label='1st and 3rd Component')
	leg = plt.legend(loc='upper right', fancybox=True)
	leg.get_frame().set_alpha(0.7)
	plt.show()

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

##################PCA#################
pca = PCA(n_components=2)
pca.fit(measured_spectra)
measured_pca = pca.transform(measured_spectra)
print("original dimension: ", measured_spectra.shape)
print("pca-projected dimension ", measured_pca.shape)
print("pca explained variance ", pca.explained_variance_)
plot_da(measured_pca, 'Principal Component Analysis')
print('Reduced dimensions: ', measured_pca.shape)

##################SVD########
measured_svd = TruncatedSVD(n_components=3, random_state=42).fit_transform(measured_spectra)
plot_da3(measured_svd, 'SVD Component Analysis') 
plot_da_components(measured_svd, 'SVD Components')

##################FactorAnalysis########
measured_fa = FactorAnalysis(n_components = 2).fit_transform(measured_spectra) #without labels 
plot_da(measured_fa, 'Factor Analysis') 
print('Reduced dimensions: ', measured_fa.shape)

measured_fa = FactorAnalysis(n_components = 3).fit_transform(measured_spectra, labels) #with labels 
plot_da3(measured_fa, 'Factor Analysis') 
print('Reduced dimensions: ', measured_fa.shape)
plot_da_components(measured_fa, 'Factor Analysis Components')

##################ICA########
ica = FastICA(n_components=3, random_state=12, max_iter=500, tol=0.001) 
measured_ica=ica.fit_transform(measured_spectra)
plot_da(measured_ica, 'Independent Component Analysis')
plot_da3(measured_ica, 'Independent Component Analysis') 

##################ISOMAP########
measured_isomap = manifold.Isomap(n_neighbors=5, n_components=3, n_jobs=-1).fit_transform(measured_spectra)
plot_da_components(measured_isomap, 'ISOMAP Components')
plot_da3(measured_isomap, 'ISOMAP Component Analysis') 

##################t-SNE########
measured_tsne = manifold.TSNE(n_components=3, n_iter=300).fit_transform(measured_spectra)
plot_da_components(measured_tsne, 't-SNE Components')
plot_da3(measured_tsne, 't-SNE Component Analysis') 


##################LDA#################
lda = LDA(n_components=2, solver="svd", store_covariance=True)
#y_pred = lda.fit(measured_spectra, labels).predict(measured_spectra)
measured_lda = lda.fit(measured_spectra, labels).transform(measured_spectra)
print('Reduced dimensions: ', measured_lda.shape)
#plot_da(measured_lda, 'Linear Discriminant Analysis')
 
##################QDA#################
qda = QDA()
measured_qda = qda.fit(measured_spectra, labels).predict(measured_spectra)
#plot_da(measured_qda, 'Quadratic Discriminant Analysis')


#from sklearn.model_selection import train_test_split
#X_train, X_test, y_train, y_test = train_test_split(X,y, test_size=0.5)
