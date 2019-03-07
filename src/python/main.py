from os.path import dirname, join as pjoin
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis as QDA
from sklearn.ensemble import RandomForestRegressor
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
	if 'Principal' in title:
		tag = 'PC'
	elif 'Linear' in title:
		tag = 'LD'
	elif 'Quadratic' in title:
		tag = 'QD'
	else:
		tag = ''

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

##################PCA#################
pca = PCA(n_components=2)

#Measured Spectra PCA 
pca.fit(measured_spectra)
measured_pca = pca.transform(measured_spectra)
print("original dimension: ", measured_spectra.shape)
print("pca-projected dimension ", measured_pca.shape)
print("pca explained variance ", pca.explained_variance_)
plot_da(measured_pca, 'Principal Component Analysis')

##################RandomForest########
model = RandomForestRegressor(random_state=1, max_depth=10)
model.fit(measured_spectra,labels)
features = np.linspace(380, 780, 81)
importances = model.feature_importances_
indices = np.argsort(importances)[-9:]  # top 10 features
mean_importance = np.mean(importances)
print(mean_importance)
print(indices)
plt.figure()
plt.title('Feature Importances')
plt.barh(range(len(indices)), importances[indices], color='b', align='center')
plt.yticks(range(len(indices)), [features[i] for i in indices])
plt.xlabel('Relative Importance')
plt.show()

##################LDA#################
lda = LDA(n_components=2, solver="svd", store_covariance=True)
#y_pred = lda.fit(measured_spectra, labels).predict(measured_spectra)
measured_lda = lda.fit(measured_spectra, labels).transform(measured_spectra)
#plot_da(measured_lda, 'Linear Discriminant Analysis')
 
##################QDA#################
qda = QDA()
measured_qda = qda.fit(measured_spectra, labels).predict(measured_spectra)
#plot_da(measured_qda, 'Quadratic Discriminant Analysis')

#head, *tail = Spectra



#from sklearn.model_selection import train_test_split
#X_train, X_test, y_train, y_test = train_test_split(X,y, test_size=0.5)
