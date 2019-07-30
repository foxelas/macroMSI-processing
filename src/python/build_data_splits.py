from os.path import dirname, join as pjoin, exists
from sklearn.model_selection import train_test_split, StratifiedKFold
import scipy.io as sio 
import numpy as np
import csv
import math

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
names = np.array(['_'.join([id_struct[x].SpectrumFile, id_struct[x].T]) for x in range(msiN) ])
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
	head(data_s)
	return data_s, labels_s, fixation_s, subset_ids

def get_subset_contents(name, data, labels):
	subset_ids = subset_indexes(name, data)
	if (len(data.shape) > 1):
		data_s = data[subset_ids,:]
	else:
		data_s = data[subset_ids]

	contents = [[x,y] for (x,y) in zip(subset_ids, data_s)]
	return contents

def head(values, n = 10):
	print('\n--------------------HEAD--------------------')
	[print(values[i]) for i in range(n)]
	print('....')

def write_file(filename, contents):

	with open(pjoin(in_dir, filename), 'w') as csvFile:
		writer = csv.writer(csvFile)
		writer.writerows(contents)
	csvFile.close()


measured_names, measured_ids = np.unique(names, return_index=True, axis=0)

write_file( 'names_mixed.csv', get_subset_contents('unique', names, positive_labels))
write_file( 'names_unfixed.csv', get_subset_contents('unfixed', names, positive_labels))
write_file( 'names_fixed.csv', get_subset_contents('fixed', names, positive_labels))
write_file( 'names_cut.csv', get_subset_contents('cut', names, positive_labels))
write_file( 'names_unique_unfixed.csv', get_subset_contents('unique_unfixed', names, positive_labels))

samples = np.unique(np.array([id_struct[x].Sample for x in range(msiN) ])).tolist()
folds = 10
subset = get_subset_contents('unique_unfixed', names, positive_labels)

start_index = 0
for i in range(folds):
	end_index = math.floor(len(samples) / folds) + start_index if i != folds - 1 else len(samples)
	fold_contents = [ subset[j][:] for j in range(len(subset)) for x in samples[start_index:end_index] if x in subset[j][1]]
	start_index = end_index
	write_file( 'names_fold' + str(i) + '.csv', fold_contents)
