import data_handling as dh

from os import mkdir, makedirs
from os.path import dirname, join as pjoin, exists
from sklearn import manifold
from sklearn.decomposition import PCA, FactorAnalysis, TruncatedSVD, FastICA 
from sklearn.preprocessing import StandardScaler
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA, QuadraticDiscriminantAnalysis as QDA
from sklearn.feature_selection import SelectFromModel
from sklearn.ensemble import RandomForestRegressor
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib import interactive, is_interactive, rcParams
import matplotlib.pyplot as plt
rcParams.update({'font.size': 20})


interactive(True)

out_dir = dh.get_out_dir()
label_dict = dh.get_label_dict()

####################dimension reduction#######################
# Dimension reduction is trained based on the measured set 

def plot_da(X, X_labels, tag, title):
	fig = plt.figure()
	#tag = title.strip('Analysis')
	n_components = X.shape[1]

	if n_components == 1 : #When X is reduced to 1 dimension
		for label, marker, color in zip(range(2), ('s', 'o'), ('green', 'red')):
			plt.scatter(x=np.where(X_labels == label)[0],  
						y=X[:,0][X_labels == label], 
						marker=marker,
						color=color,
						alpha=0.7,
						label=label_dict[label])
		plt.xlabel('Samples')
		plt.ylabel(tag)

	elif n_components == 2  or n_components > 3: #when X is reduced to 2 dimensions
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
		print('Unsupported plot type for these dimensions.')
		return 

	leg = plt.legend(loc='upper right', fancybox=True)
	leg.get_frame().set_alpha(0.7)
	#plt.title(title)
	plt.grid()
	plt.tight_layout
	plt.show()
	plt.savefig(pjoin(out_dir, 'dimension_reduction', title + tag + '.png'), dpi=200, bbox_inches='tight')


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
	print(out_dir)
	print(range(len(indices)), [features[i] for i in indices])
	print(importances[indices])
	plt.savefig(pjoin(out_dir, 'dimension_reduction', title + '.png'), bbox_inches='tight')

def get_dimred(method, components=2):
	if method == 'RF':		
		##################RandomForest########
		dimred_obj = RandomForestRegressor(random_state=2, max_depth=10, n_estimators=10)

	elif method == 'PCA': 
		##################PCA#################
		if (components == -1):
			dimred_obj = PCA(random_state=2)
		else:
			dimred_obj = PCA(n_components=components, random_state=2)

	elif method == 'SVD': 
		##################SVD########
		dimred_obj = TruncatedSVD(n_components=components, random_state=2)

	elif method == 'FA':
		##################FactorAnalysis########
		dimred_obj = FactorAnalysis(n_components = components)

	elif method == 'ICA':
		##################ICA########
		if (components == -1):
			dimred_obj = FastICA(random_state=2, max_iter=500, tol=0.1)
		else:
			dimred_obj = FastICA(n_components=components, random_state=2, max_iter=500, tol=0.1)

	elif method == 'ISOMAP':
		##################ISOMAP########
		dimred_obj = manifold.Isomap(n_neighbors=5, n_components=components, n_jobs=-1)

	elif method == 't-SNE':
		##################t-SNE########
		dimred_obj = manifold.TSNE(n_components=components, n_iter=300)

	elif method == 'LDA':
		##################LDA#################
		dimred_obj = LDA(n_components=components, solver="svd", store_covariance=True)

	elif method == 'QDA':
		##################QDA#################
		dimred_obj = QDA()

	elif method == 'None':
		dimred_obj = None

	else:
		print(method, ' not implemented.')
		return	

	return dimred_obj

def _get_ica_map(ica):
    """Get ICA topomap for components"""
    fast_dot = np.dot
    maps = fast_dot(ica.components_, ica.mixing_)
    return maps

def dimension_reduction(data, data_labels, method=None, show_figures=False, components=2, name=''):

	if method is None: 
		return data

	dh.create_directory(pjoin(out_dir, 'dimension_reduction'))

	dimred_obj = get_dimred(method, components)
	if method == 'RF':		
		##################RandomForest########
		rf = dimred_obj.fit(data,data_labels)
		features = np.linspace(380, 780, 81)
		importances = rf.feature_importances_
		indices = np.argsort(importances)[-9:]  # top 10 features
		mean_importance = np.mean(importances)
		print('Random Forest-Average wavelength importance: ', mean_importance)
		if show_figures:
			plot_feature_importance(features, importances, indices, 'Wavelength Importances by Random Forest' + '(' + name + ')')

		fitted_dimred = SelectFromModel(rf).fit(data, data_labels)
		data_transformed = fitted_dimred.transform(data)

	elif method == 'PCA': 
		##################PCA#################
		fitted_dimred = dimred_obj.fit(data)
		data_transformed = fitted_dimred.transform(data)
		if show_figures:
			print("PCA-explained variance ", fitted_dimred.explained_variance_ratio_)
			plot_da(data_transformed, data_labels,'PC', name)

	elif method == 'SVD': 
		##################SVD########
		fitted_dimred = dimred_obj.fit(data)
		data_transformed = fitted_dimred.transform(data)
		if show_figures:
			plot_da(data_transformed, data_labels, 'SVD C', name) 
			plot_da_components(data_transformed, data_labels, 'SVD Component Analysis' + '(' + name + ')')

	elif method == 'FA':
		##################FactorAnalysis########
		fitted_dimred = dimred_obj.fit(data, data_labels)
		data_transformed = fitted_dimred.transform(data) #with labels 
		if show_figures:
			plot_da(data_transformed, data_labels, 'FA' , name ) 
			plot_da_components(data_transformed, data_labels, 'Factor Analysis Components' + '(' + name + ')')

	elif method == 'ICA':
		##################ICA########
		fitted_dimred = dimred_obj.fit(data)
		data_transformed=fitted_dimred.transform(data)
		if show_figures:
			#print(_get_ica_map(fitted_dimred))
			#print(fitted_dimred.components_ , fitted_dimred.components_.shape)
			#print(fitted_dimred.mixing_)
			#print(fitted_dimred.get_params())
			plot_da(data_transformed, data_labels, 'IC', name)

	elif method == 'ISOMAP':
		##################ISOMAP########
		fitted_dimred = dimred_obj.fit(data)
		data_transformed = fitted_dimred.transform(data)
		if show_figures:
			plot_da_components(data_transformed, data_labels, 'ISOMAP Components' + '(' + name + ')')
			plot_da(data_transformed, data_labels, 'ISOMAP C', name) 

	elif method == 't-SNE':
		##################t-SNE########
		fitted_dimred = dimred_obj.fit(data)
		data_transformed = fitted_dimred.fit_transform(data)
		if show_figures:
			plot_da_components(data_transformed, data_labels, 't-SNE Components' + '(' + name + ')')
			plot_da(data_transformed, data_labels, 't-SNE C', name) 

	elif method == 'LDA':
		##################LDA#################
		fitted_dimred = dimred_obj.fit(data, data_labels)
		data_transformed = fitted_dimred.transform(data)
		if show_figures:
			plot_da(data_transformed, data_labels, 'LD', name)

	elif method == 'QDA':
		##################QDA#################
		fitted_dimred = dimred_obj.fit(data, data_labels)
		data_transformed = np.array([[x] for x in fitted_dimred.decision_function(data)])
		plot_da(data_transformed, data_labels, 'QD', name)

	elif method == 'None':
		fitted_dimred = None
		data_transformed = data
	else:
		print(method, ' not implemented.')
		return

	if show_figures:
		print(method + '-Reduced dimensions: ', data_transformed.shape)

	return data_transformed, fitted_dimred


def reduce(data, fitted_dimred=None):
	if fitted_dimred is None:
		return data
	else:
		return fitted_dimred.transform(data)

def fit_and_reduce(train, test, train_labels, method=None, n_components=2):
	train_transformed, fitted_dimred = dimension_reduction(train, train_labels, method, False, n_components)
	test_transformed = reduce(test, fitted_dimred)
	return train_transformed, test_transformed

def compare_dimension_reduction(data, labels, show_figures=False, name=''):
	dimension_reduction(data, labels, 'RF', show_figures, 2, name)
	dimension_reduction(data, labels, 'PCA', show_figures, 10, name)
	dimension_reduction(data, labels, 'SVD', show_figures, 3, name)
	dimension_reduction(data, labels, 'FA', show_figures, 3, name)
	dimension_reduction(data, labels, 'ICA', show_figures, 10, name)
	dimension_reduction(data, labels, 'ISOMAP', show_figures, 3, name)
	dimension_reduction(data, labels, 't-SNE', show_figures, 3, name)
	dimension_reduction(data, labels, 'LDA', show_figures, 1, name)
	dimension_reduction(data, labels, 'QDA', show_figures, 1, name)

def compare_pca_ica(data, labels, show_figures=False, name=''):
	dimension_reduction(data, labels, 'PCA', show_figures, 20, name)
	dimension_reduction(data, labels, 'ICA', show_figures, 20, name)

# scaler = dh.get_scaler(dh.get_measured_spectra())
# current_data, current_labels, current_fixation, current_indexes = dh.get_scaled_subset('unique', dh.get_measured_spectra(), dh.get_labels(), scaler)
# print('Subset contains ', len(current_labels), ' observations')

compare_pca_ica(dh.get_measured_spectra(),  dh.get_labels(), True, 'measured')
#compare_pca_ica(dh.get_reconstructed_spectra(),  dh.get_labels(), True, 'reconstructed')
#compare_pca_ica(dh.get_concat_lbp('mmlbp'),  dh.get_labels(), True, 'mmlbp')

#dimension_reduction(dh.get_measured_spectra(), dh.get_measured_spectra(), 'RF', True, 2, 'measured')


