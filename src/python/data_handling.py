from os import mkdir, makedirs
from os.path import dirname, join as pjoin, exists
import scipy.io as sio 
import numpy as np
import csv
import math
import datetime
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

def savemat(filename, mdict):
	sio.savemat(filename, mdict)

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

############## directories #####################
def create_directory(dirName):
	if not exists(dirName):
	    makedirs(dirName)

def get_base_dir():
	base_dir  = '/media/sf_research/input/'
	return base_dir

def get_data_dir():
	data_dir = 'saitama_v8_min_region_bright'
	return data_dir

def get_out_dir():
	data_dir = get_data_dir()
	out_dir = pjoin('/media/sf_research/output/', data_dir)
	create_directory(out_dir)
	return out_dir

def get_in_dir():
	base_dir = get_base_dir()
	data_dir = get_data_dir()
	in_dir = pjoin(base_dir, data_dir)
	return in_dir

def get_log_dir():
	log_dir = '/media/sf_research/macroMSI-processing/logs/'
	return log_dir

in_dir = get_in_dir()
in_mat = loadmat( pjoin(in_dir, 'in.mat'))
print('Finished loading input matfile.')#print('Loaded mat file with headers', in_mat.keys())

out_dir = get_out_dir()
out_mat = loadmat(pjoin(out_dir, '8-Features', 'out.mat'))
print('Finished loading output matfile.')

in_dir = get_in_dir()
matfile = pjoin(in_dir, 'ID.mat')
id_mat = loadmat(matfile)
print('Finished loading ID matfile.')
id_struct = id_mat['ID']

###################gets#####################

def get_label_dict():
	label_dict = {0: 'Benign', 1: 'Malignant'}
	return label_dict

def load_in_mat_internal():
	return in_mat

def load_out_mat_internal():
	return out_mat

in_mat_i = load_in_mat_internal()
def load_in_mat():
	return in_mat_i

out_mat_i = load_out_mat_internal()
def load_out_mat():
	return out_mat_i

def get_measured_spectra_internal():
	in_mat = load_in_mat()
	spectra = in_mat['Spectra']
	return spectra

def get_reconstructed_spectra_internal():
	out_mat = load_out_mat()
	estimated_spectra = out_mat['EstimatedSpectra']
	return estimated_spectra

def get_reconstructed_spectra_rgb_internal():
	out_mat = load_out_mat()
	estimated_spectra = out_mat['EstimatedRGBSpectra']
	return estimated_spectra

def get_conc_lbp_internal():
	out_mat = load_out_mat()
	lbp_features = out_mat['ConcatLbpFeatures']
	return lbp_features

def get_sum_lbp_internal():
	out_mat = load_out_mat()
	lbp_features = out_mat['SumLbpFeatures']
	return lbp_features

def get_multispectral_lbp_internal():
	out_mat = load_out_mat()
	lbp_features = out_mat['MMLbpFeatures']
	return lbp_features

def get_rgb_lbp_internal():
	out_mat = load_out_mat()
	lbp_features = out_mat['RgbLbpFeatures']
	return lbp_features

def get_complete_spectra_internal():
	in_mat = load_in_mat()
	complete_spectra = in_mat['CompleteSpectra']
	return complete_spectra

# def get_darkIs_internal(): 
# 	in_mat = load_in_mat()
# 	darkIs = in_mat['DarkIs']
# 	return darkIs

# def get_msi_names_internal(): 
# 	in_mat = load_in_mat()
# 	msi_names = in_mat['MSINames']
# 	return msi_names

# def get_msis_internal(): 
# 	in_mat = load_in_mat()
# 	msis = in_mat['MSIs']
# 	return msis

# def get_masks_internal(): 
# 	in_mat = load_in_mat()
# 	masks = in_mat['Masks']
# 	return masks

def get_spectra_names_internal():
	in_mat = load_in_mat()
	spectra_names = in_mat['SpectraNames']
	return spectra_names

# def get_whiteIs_internal():
# 	in_mat = load_in_mat()
# 	whiteIs = in_mat['WhiteIs']
# 	return whiteIs

#######################Use of internals#########################
estimated_spectra_i = get_reconstructed_spectra_internal()
def get_reconstructed_spectra():
	return estimated_spectra_i

measured_spectra_i = get_measured_spectra_internal()
def get_measured_spectra():
	return measured_spectra_i

estimated_spectra_rgb_i = get_reconstructed_spectra_rgb_internal()
def get_reconstructed_spectra_rgb():
	return estimated_spectra_rgb_i

conc_lbp_features_i = get_conc_lbp_internal()
def get_conc_lbp():
	return conc_lbp_features_i

sum_lbp_features_i = get_sum_lbp_internal()
def get_sum_lbp():
	return sum_lbp_features_i

m_lbp_features_i = get_multispectral_lbp_internal()
def get_multispectral_lbp():
	return m_lbp_features_i

rgb_lbp_features_i = get_rgb_lbp_internal()
def get_rgb_lbp():
	return rgb_lbp_features_i

complete_spectra_i = get_complete_spectra_internal()
def get_complete_spectra():
	return complete_spectra_i

# darkIs_i = get_darkIs_internal()
# def get_darkIs(): 
# 	return darkIs_i

# msi_names_i = get_msi_names_internal()
# def get_msi_names(): 
# 	return msi_names_i

# msis_i = get_msis_internal()
# def get_msis(): 
# 	return msis_i

# masks_i = get_masks_internal()
# def get_masks(): 
# 	return masks_i

spectra_names_i = get_spectra_names_internal()
def get_spectra_names():
	return spectra_names_i

# whiteIs_i = get_whiteIs_internal()
# def get_whiteIs():
# 	return whiteIs_i

def concat_features(feat1, feat2=None, feat3=None):
	if feat1 is None and feat2 is not None:
		return concat_features(feat2, feat3)
	if feat3 is not None: 
		return np.array([np.concatenate((x,y,z),0) for (x,y,z) in zip(feat1, feat2, feat3)])
	if feat2 is not None:
		return np.array([np.concatenate((x,y),0) for (x,y) in zip(feat1, feat2)])
	else:
		return np.array(feat1)

def get_concat_lbp(feature_set, rgb=False):
	if rgb:
		lbp_features = get_rgb_lbp()

	elif "slbp" in feature_set:
		lbp_features = get_sum_lbp()

	elif "clbp" in feature_set:
		lbp_features = get_conc_lbp()

	elif "mlbp" in feature_set:
		lbp_features = get_multispectral_lbp()

	concat_lbp = concat_features(lbp_features[0], lbp_features[1], lbp_features[2])
	#concat_lbp = concat_features(lbp_features[0], lbp_features[1], lbp_features[2])
	return concat_lbp


def get_data_dimensions():
	complete_spectra = get_complete_spectra()
	msis = get_msis()
	msiN, lambdaN = complete_spectra.shape
	msi_dim = msis[0].shape
	bandN = msi_dim[0]
	rgbN = msi_dim[3]
	return msiN, lambdaN, bandN, rgbN

def get_labels():
	id_struct = get_id_struct()
	is_benign = np.array([id_struct[x].IsBenign for x in range(len(id_struct)) ])
	positive_labels = np.array([int((x + 1) % 2) for x in is_benign])
	return positive_labels

def get_log_file():
	return pjoin(get_log_dir(), datetime.datetime.now().strftime("%Y-%m-%d %H_%M") + '.log')

def get_id_struct():
	return id_struct

def get_fixation():
	id_struct = get_id_struct()
	is_cut = np.array([id_struct[x].IsCut for x in range(len(id_struct)) ])
	is_fixed = np.array([id_struct[x].IsFixed for x in range(len(id_struct)) ])

	fixation = []
	for i in range(len(id_struct)):
		if is_fixed[i] and not(is_cut[i]):
			fixation.append('fixed')	
		elif is_cut[i]:
			fixation.append('cut')
		else:
			fixation.append('unfixed')
	return fixation 

def get_samples():
	id_struct = get_id_struct()
	samples = np.unique(np.array([id_struct[x].Sample for x in range(len(id_struct)) ])).tolist()
	return samples

def get_samples_all():
	id_struct = get_id_struct()
	samples = [id_struct[x].Sample for x in range(len(id_struct)) ]
	return samples

def subset_indexes(name):
	subset_ids = []
	if name == 'unique':
		id_struct = get_id_struct()
		names = np.array(['_'.join([id_struct[x].SpectrumFile, id_struct[x].T]) for x in range(len(id_struct)) ])
		spectra_names = get_spectra_names()
		x, subset_ids = np.unique(np.array(spectra_names, dtype=str), return_index=True, axis=0)

	elif name == 'unfixed':
		fixation = get_fixation()
		subset_ids = np.array(get_indexes_equal(name, fixation))

	elif name == 'fixed':
		fixation = get_fixation()
		subset_ids_f = np.array(get_indexes_equal('fixed', fixation))
		subset_ids_c = np.array(get_indexes_equal('cut', fixation))
		subset_ids = np.concatenate((subset_ids_f, subset_ids_c), axis=0) #axis=0

	elif 'unique' in name:
		name1, name2 = name.split('_')
		subset_ids1 = subset_indexes(name1)
		subset_ids2 = subset_indexes(name2)
		subset_ids = np.intersect1d(subset_ids1, subset_ids2)

	else: 
		print('Not implemented yet.')
		return 

	return subset_ids

def get_subset_with_index(indexes, data, labels, feature_set="spect", rgb=False):

	if (len(data.shape) > 1):
		data_s = data[indexes,:]
	else:
		data_s = data[indexes]
	labels_s = labels[indexes]
	fixation = get_fixation()
	fixation_s = [fixation[x] for x in indexes]
	if "lbp" in feature_set:
		lbp_features = get_concat_lbp(feature_set, rgb)
		lbp_features_s = np.array([lbp_features[x] for x in indexes])
	else:
		lbp_features_s = None

	return data_s, labels_s, fixation_s, lbp_features_s

def get_subset(name, data, labels):
	subset_ids = subset_indexes(name, data)
	data_s, labels_s, fixation_s = get_subset_with_index(subset_ids, data, labels)

	#sprint('Created subset ' + name + '')
	#head(data_s)
	return data_s, labels_s, fixation_s, subset_ids

def get_subset_contents(name, data, labels):
	subset_ids = subset_indexes(name)
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

def write_file(filename, contents, write_options="w+"):
	with open(filename, write_options) as csv_file:
		writer = csv.writer(csv_file)
		writer.writerows(contents)
	csv_file.close()

def write_log(filename, contents, write_options="a+"):
	with open(filename, write_options) as log_file:
		log_file.write(contents)
	log_file.close()

############################Data set splits############################################
test_sample_names = ['9933']

def write_folds(subset_name='unique', folds=10):
	spectra_names = get_spectra_names()
	positive_labels = get_labels()
	subset = get_subset_contents(subset_name, spectra_names, positive_labels)

	samples = get_samples()
	start_index = 0
	for i in range(folds):
		end_index = math.floor(len(samples) / folds) + start_index if i != folds - 1 else len(samples)
		fold_contents = [ subset[j][:] for j in range(len(subset)) for x in samples[start_index:end_index] if x in subset[j][1]]
		start_index = end_index
		write_file( 'names_fold' + str(i) + '.csv', fold_contents)

def get_fold_indexes(subset_name='unique', folds=10, ignore_test=False):
	spectra_names = get_spectra_names()
	positive_labels = get_labels()
	subset = get_subset_contents(subset_name, spectra_names, positive_labels)

	samples = get_samples()
	start_index = 0
	fold_indexes = []
	on_hold_indexes = []
	for i in range(folds):
		end_index = math.floor(len(samples) / folds) + start_index if i != folds - 1 else len(samples)
		if ignore_test:
			current_fold_indexes = [ subset[j][0] for j in range(len(subset)) for x in samples[start_index:end_index] if x in subset[j][1]]
		else:
			current_fold_indexes = [ subset[j][0] for j in range(len(subset)) for x in samples[start_index:end_index] if x not in test_sample_names if x in subset[j][1]]
		current_labels = [positive_labels[x] for x in current_fold_indexes ]
		if all(v == 0 for v in current_labels) or all(v == 1 for v in current_labels):
			on_hold_indexes.append(current_fold_indexes)
		else:
			fold_indexes.append(current_fold_indexes)

		start_index = end_index

	fold_indexes = [x for x in fold_indexes if x]
	on_hold_indexes = [x for x in on_hold_indexes if x]
	shortest_indexes = sorted(range(len(fold_indexes)), key=lambda k: len(fold_indexes[k]))
	shortest_on_hold_indexes = sorted(range(len(on_hold_indexes)), key=lambda k: len(on_hold_indexes[k]), reverse=True)
	for i in range(len(on_hold_indexes)):
		fold_indexes[shortest_indexes[i]] = fold_indexes[shortest_indexes[i]] + on_hold_indexes[shortest_on_hold_indexes[i]]
	folds = len(fold_indexes)

	return fold_indexes, folds

def get_fold_indexes_stratified(subset_name='unique',folds=7, ignore_test=False):
	spectra = get_measured_spectra()
	positive_labels = get_labels()
	subset_ids_s = subset_indexes(subset_name) 

	samples_all = get_samples_all()
	if ignore_test:
		subset_ids = [x for x in subset_ids_s]
	else:
		subset_ids = [x for x in subset_ids_s if samples_all[x] not in test_sample_names]
	data, labels, fixation, lbp_features = get_subset_with_index(subset_ids, spectra, positive_labels)


	cv = StratifiedKFold(n_splits=folds, shuffle=True, random_state=1)
	fold_indexes = []
	for train, test in cv.split(data, labels):
		fold_indexes.append([subset_ids[x] for x in test])

	return fold_indexes, folds



def get_test_indexes(subset_name='unique'):
	spectra_names = get_spectra_names()
	positive_labels = get_labels()
	subset = get_subset_contents(subset_name, spectra_names, positive_labels)
	test_indexes = [ subset[j][0] for j in range(len(subset)) for x in test_sample_names if x in subset[j][1] ]
	return test_indexes	

def get_scaler(data):
	return StandardScaler().fit(data)

def get_scaled_subset(subset_name, data, labels, scaler):
	data_s, labels_s, fixation_s, subset_ids = get_subset(subset_name, data, labels)
	data_s = scaler.transform(data_s)
	return data_s, labels_s, fixation_s, subset_ids

def get_scaled_subset_with_index(indexes, data, labels, scaler, feature_set, rgb=False):
	data_s, labels_s, fixation_s, lbp_features_s = get_subset_with_index(indexes, data, labels, feature_set, rgb)
	data_s = scaler.transform(data_s)
	return data_s, labels_s, lbp_features_s


#print(get_test_indexes('unique_unfixed'))
#print(get_fold_indexes('unique_unfixed', 10))
#write_file( 'names_mixed.csv', get_subset_contents('unique', spectra_names, positive_labels))
#print(get_subset_with_index(get_test_indexes('unique_unfixed'), spectra, positive_labels))

# ids,e = get_fold_indexes('unique', 10)
# [print(x) for x in ids]
# [print(positive_labels[x]) for x in ids]

# print("wqqq")
# ids,e = get_fold_indexes('unique', 10)
# [print(x) for x in ids]
# positive_labels = get_labels()
# [print(positive_labels[x]) for x in ids]