from os.path import dirname, join as pjoin, exists
from sklearn.model_selection import train_test_split, StratifiedKFold
import scipy.io as sio 
import numpy as np


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

####################load matfiles#######################

base_dir  = '/media/sf_research/input/'
data_dir = 'saitama_v5_min_region'
in_dir = pjoin(base_dir, data_dir)

in_mat = loadmat( pjoin(in_dir, 'in.mat'))
print('Finished loading input matfile.')#print('Loaded mat file with headers', in_mat.keys())

matfile = pjoin(in_dir, 'ID.mat')
id_mat = loadmat(matfile)
print('Finished loading ID matfile.')
id_struct = id_mat['ID']
msiN = id_struct.shape[0]
names = np.array([id_struct[x].SpectrumFile for x in range(msiN) ])
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

def subset_indexes(name, data):
	subset_ids = []
	if name == 'unique':
		x, subset_ids = np.unique(names, return_index=True, axis=0)

	elif name == 'unfixed' or name == 'fixed' or name == 'cut':
		subset_ids = np.array(get_indexes_equal(name, fixation))

	elif 'unique' in name:
		name1, name2 = name.split('_')
		subset_ids1 = subset_indexes(name1, data)
		subset_ids2 = subset_indexes(name2, data)
		subset_ids = np.intersect1d(subset_ids1, subset_ids2)

	else: 
		print('Not implemented yet.')
		return 

	return subset_ids

def get_subset(name, data, labels):
	subset_ids = subset_indexes(name, data)
	if (len(data.shape) > 1):
		data_s = data[subset_ids,:]
	else:
		data_s = data[subset_ids]
	labels_s = labels[subset_ids]
	fixation_s = [fixation[x] for x in subset_ids]

	print('Created subset ' + name + '')
	return data_s, labels_s, fixation_s

def head(values, n = 10):
	print('\n--------------------HEAD--------------------')
	[print(values[i]) for i in range(n)]
	print('....')

measured_names, measured_ids = np.unique(names, return_index=True, axis=0)

head(measured_names)
data, labels, fixation = get_subset('unique', names, positive_labels)
print(data.shape)
