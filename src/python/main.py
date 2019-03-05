from os.path import dirname, join as pjoin
from sklearn.decomposition import PCA
import numpy as np
import scipy.io as sio 
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import interactive

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


##################PCA#################
pca = PCA(n_components=2)

#Measured Spectra PCA 
measured_spectra, unique_ids = np.unique(spectra, return_index=True, axis=0)
benign_ids = get_indexes(1, is_benign[unique_ids])
malignant_ids = get_indexes(0, is_benign[unique_ids])
pca.fit(measured_spectra)
measured_pca = pca.transform(measured_spectra)
print("original dimension: ", measured_spectra.shape)
print("pca-projected dimension ", measured_pca.shape)
print("pca explained variance ", pca.explained_variance_)

plt.figure()
#plt.scatter(measured_pca[: , 0], measured_pca[: , 1], c=is_benign[unique_ids], cmap='seismic', alpha=0.8)
plt.scatter(measured_pca[benign_ids , 0], measured_pca[benign_ids , 1], c='green', alpha=0.7, label='Benign')
plt.scatter(measured_pca[malignant_ids , 0], measured_pca[malignant_ids , 1], c='red', alpha=0.7, label='Malignant')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.title('Principal Component Analysis')
plt.legend('Benign', 'Malignant')
plt.legend()



#head, *tail = Spectra



#from sklearn.model_selection import train_test_split
#X_train, X_test, y_train, y_test = train_test_split(X,y, test_size=0.5)
