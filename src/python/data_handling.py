from os import mkdir, makedirs
from os.path import dirname, join as pjoin, exists
import scipy.io as sio 
import numpy as np
import csv
import math
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split, StratifiedKFold


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
in_dir = pjoin(base_dir, data_dir)
out_dir = pjoin('/media/sf_research/output/', data_dir)
create_directory(out_dir)

in_mat = loadmat( pjoin(in_dir, 'in.mat'))
print('Finished loading input matfile.')#print('Loaded mat file with headers', in_mat.keys())
out_mat = loadmat(pjoin(out_dir, 'ReflectanceEstimationPreset', 'out.mat'))
print('Finished loading output matfile.')

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

estimated_spectra = out_mat['EstimatedSpectra']

matfile = pjoin(in_dir, 'ID.mat')
id_mat = loadmat(matfile)
print('Finished loading ID matfile.')
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

samples = np.unique(np.array([id_struct[x].Sample for x in range(msiN) ])).tolist()

###################gets#####################
def get_out_dir():
	return out_dir

def get_in_dir():
	return in_dir

def get_label_dict():
	return label_dict

def get_fixation():
	return fixation

def get_measured_spectra():
	return spectra

def get_reconstructed_spectra():
	return estimated_spectra

def get_labels():
	return positive_labels

def subset_indexes(name, data):
	subset_ids = []
	if name == 'unique':
		names = np.array(['_'.join([id_struct[x].SpectrumFile, id_struct[x].T]) for x in range(msiN) ])
		x, subset_ids = np.unique(np.array(spectra_names, dtype=str), return_index=True, axis=0)

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

def get_subset_with_index(indexes, data, labels):
	if (len(data.shape) > 1):
		data_s = data[indexes,:]
	else:
		data_s = data[indexes]
	labels_s = labels[indexes]
	fixation_s = [fixation[x] for x in indexes]

	return data_s, labels_s, fixation_s

def get_subset(name, data, labels):
	subset_ids = subset_indexes(name, data)
	data_s, labels_s, fixation_s = get_subset_with_index(subset_ids, data, labels)

	#sprint('Created subset ' + name + '')
	#head(data_s)
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


############################Data set splits############################################
test_sample_names = ['9933']

def write_folds(subset_name='unique', folds=10):
	subset = get_subset_contents(subset_name, spectra_names, positive_labels)

	start_index = 0
	for i in range(folds):
		end_index = math.floor(len(samples) / folds) + start_index if i != folds - 1 else len(samples)
		fold_contents = [ subset[j][:] for j in range(len(subset)) for x in samples[start_index:end_index] if x in subset[j][1]]
		start_index = end_index
		write_file( 'names_fold' + str(i) + '.csv', fold_contents)

def get_fold_indexes(subset_name='unique', folds=10):
	subset = get_subset_contents(subset_name, spectra_names, positive_labels)

	start_index = 0
	fold_indexes = []
	for i in range(folds):
		end_index = math.floor(len(samples) / folds) + start_index if i != folds - 1 else len(samples)
		fold_indexes.append([ subset[j][0] for j in range(len(subset)) for x in samples[start_index:end_index] if x not in test_sample_names if x in subset[j][1]])
		start_index = end_index

	fold_indexes[-2] = fold_indexes[-2] + fold_indexes[1]
	fold_indexes[1] = [] 
	fold_indexes[3] = fold_indexes[3] + fold_indexes[2]
	fold_indexes[2] = []
	fold_indexes = [x for x in fold_indexes if x]
	folds = len(fold_indexes)

	return fold_indexes, folds

def get_fold_indexes_stratified(subset_name='unique',folds=10):
	data, labels, fixation_s, subset_ids = get_subset(subset_name, spectra, positive_labels)
	cv = StratifiedKFold(n_splits=folds, shuffle=True, random_state=1)
	fold_indexes = []
	for train, test in cv.split(data, labels):
		fold_indexes.append(subset_ids[test].tolist())

	return fold_indexes, folds

def get_test_indexes(subset_name='unique'):
	subset = get_subset_contents(subset_name, spectra_names, positive_labels)
	test_indexes = [ subset[j][0] for j in range(len(subset)) for x in test_sample_names if x in subset[j][1] ]
	return test_indexes	

def get_scaler(data):
	return StandardScaler().fit(data)

def get_scaled_subset(subset_name, data, labels, scaler):
	data_s, labels_s, fixation_s, subset_ids = get_subset(subset_name, data, labels)
	data_s = scaler.transform(data_s)
	return data_s, labels_s, fixation_s, subset_ids

def get_scaled_subset_with_index(indexes, data, labels, scaler):
	data_s, labels_s, fixation_s = get_subset_with_index(indexes, data, labels)
	data_s = scaler.transform(data_s)
	return data_s, labels_s

#print(get_test_indexes('unique_unfixed'))
#print(get_fold_indexes('unique_unfixed', 10))
#write_file( 'names_mixed.csv', get_subset_contents('unique', spectra_names, positive_labels))
#print(get_subset_with_index(get_test_indexes('unique_unfixed'), spectra, positive_labels))

# ids,e = get_fold_indexes('unique', 10)
# [print(x) for x in ids]
# [print(positive_labels[x]) for x in ids]

# print("wewe")
# ids,e = get_fold_indexes_stratified('unique', 10)
# [print(x) for x in ids]
# [print(positive_labels[x]) for x in ids]
# print("wqqq")
# ids,e = get_fold_indexes_stratified('unique', 10)
# [print(x) for x in ids]
# [print(positive_labels[x]) for x in ids]