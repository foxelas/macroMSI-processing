import data_handling as dh

from os import mkdir, makedirs
from os.path import dirname, join as pjoin, exists
from sklearn import manifold
import sklearn.metrics
from sklearn.decomposition import PCA, FactorAnalysis, TruncatedSVD, FastICA 
from sklearn.preprocessing import StandardScaler
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA, QuadraticDiscriminantAnalysis as QDA
from sklearn.neighbors import KNeighborsClassifier as KNN
from sklearn.feature_selection import SelectFromModel
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import interactive, is_interactive
import time 
import csv
import math

interactive(True)

out_dir = dh.get_out_dir()
label_dict = dh.get_label_dict()

####################dimension reduction#######################
# Dimension reduction is trained based on the measured set 

def plot_da(X, X_labels, title):
	fig = plt.figure()
	tag = title.strip('Analysis')
	n_components = X.shape[1]

	if n_components == 1 : #When X is reduced to 1 dimension
		positives_end = np.sum(X_labels==0) + 1
		negatives_end = np.sum(X_labels==1) + 1
		for label, marker, color, sample_start, sample_end in zip(range(2), ('s', 'o'), ('green', 'red'), (1, positives_end + 1), (positives_end, positives_end + negatives_end)):
			plt.scatter(x=np.array(range(sample_start, sample_end)),  
						y=X[:,0][X_labels == label], 
						marker=marker,
						color=color,
						alpha=0.7,
						label=label_dict[label])
		plt.xlabel('Samples')
		plt.ylabel(tag)

	elif n_components == 2 : #when X is reduced to 2 dimensions
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
		print('Unsupported plot type for input with greater than 3 dimensions.')
		return 

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

def plot_feature_importance(features, importances, indices, title):
	plt.figure()
	plt.title(title)
	plt.barh(range(len(indices)), importances[indices], color='b', align='center')
	plt.yticks(range(len(indices)), [features[i] for i in indices])
	plt.xlabel('Relative Importance')
	plt.show()
	plt.savefig(pjoin(out_dir, 'dimension_reduction', title + '.png'), bbox_inches='tight')

def dimension_reduction(data, data_labels, method, show_figures=False, components=2, name=''):
	dh.create_directory(pjoin(out_dir, 'dimension_reduction'))

	if method == 'RF':		
		##################RandomForest########
		rf = RandomForestRegressor(random_state=1, max_depth=10, n_estimators=10)
		rf.fit(data,data_labels)
		features = np.linspace(380, 780, 81)
		importances = rf.feature_importances_
		indices = np.argsort(importances)[-9:]  # top 10 features
		mean_importance = np.mean(importances)
		print('Random Forest-Average wavelength importance: ', mean_importance)
		if show_figures:
			plot_feature_importance(features, importances, indices, 'Wavelength Importances by Random Forest' + '(' + name + ')')

		dimred = SelectFromModel(rf).fit(data, data_labels)
		data_transformed = dimred.transform(data)

	elif method == 'PCA': 
		##################PCA#################
		dimred = PCA(n_components=components).fit(data)
		data_transformed = dimred.transform(data)
		if show_figures:
			print("PCA-explained variance ", dimred.explained_variance_ratio_)
			plot_da(data_transformed, data_labels,'Principal Component Analysis' + '(' + name + ')')

	elif method == 'SVD': 
		##################SVD########
		dimred = TruncatedSVD(n_components=components, random_state=42).fit(data)
		data_transformed = dimred.transform(data)
		if show_figures:
			plot_da(data_transformed, data_labels, 'SVD Component Analysis' + '(' + name + ')') 
			plot_da_components(data_transformed, data_labels, 'SVD Components' + '(' + name + ')')

	elif method == 'FA':
		##################FactorAnalysis########
		dimred = FactorAnalysis(n_components = components).fit(data, data_labels)
		data_transformed = dimred.transform(data) #with labels 
		if show_figures:
			plot_da(data_transformed, data_labels, 'Factor Analysis' + '(' + name + ')') 
			plot_da_components(data_transformed, data_labels, 'Factor Analysis Components' + '(' + name + ')')

	elif method == 'ICA':
		##################ICA########
		dimred = FastICA(n_components=components, random_state=1, max_iter=500, tol=0.01).fit(data)
		data_transformed=dimred.transform(data)
		if show_figures:
			plot_da(data_transformed, data_labels, 'Independent Component Analysis' + '(' + name + ')')

	elif method == 'ISOMAP':
		##################ISOMAP########
		dimred = manifold.Isomap(n_neighbors=5, n_components=components, n_jobs=-1).fit(data)
		data_transformed = dimred.transform(data)
		if show_figures:
			plot_da_components(data_transformed, data_labels, 'ISOMAP Components' + '(' + name + ')')
			plot_da(data_transformed, data_labels, 'ISOMAP Component Analysis' + '(' + name + ')') 

	elif method == 't-SNE':
		##################t-SNE########
		dimred = manifold.TSNE(n_components=components, n_iter=300).fit(data)
		data_transformed = dimred.fit_transform(data)
		if show_figures:
			plot_da_components(data_transformed, data_labels, 't-SNE Components' + '(' + name + ')')
			plot_da(data_transformed, data_labels, 't-SNE Component Analysis' + '(' + name + ')') 

	elif method == 'LDA':
		##################LDA#################
		dimred = LDA(n_components=components, solver="svd", store_covariance=True).fit(data, data_labels)
		data_transformed = dimred.transform(data)
		if show_figures:
			plot_da(data_transformed, data_labels, 'Linear Discriminant Analysis' + '(' + name + ')')

	elif method == 'QDA':
		##################QDA#################
		dimred = QDA().fit(data, data_labels)
		data_transformed = np.array([[x] for x in dimred.decision_function(data)])
		plot_da(data_transformed, data_labels, 'Quadratic Discriminant Analysis' + '(' + name + ')')

	else:
		print(method, ' not implemented.')
		return

	if show_figures:
		print(method + '-Reduced dimensions: ', data_transformed.shape)

	return data_transformed, dimred


def reduce(data, dimred_obj=None):
	if dimred_obj is None:
		return data
	else:
		return dimred_obj.transform(data)

def compare_dimension_reduction(data, labels, show_figures=False, name=''):
	dimension_reduction(data, labels, 'RF', show_figures, 2, name)
	dimension_reduction(data, labels, 'PCA', show_figures, 2, name)
	dimension_reduction(data, labels, 'SVD', show_figures, 3, name)
	dimension_reduction(data, labels, 'FA', show_figures, 3, name)
	dimension_reduction(data, labels, 'ICA', show_figures, 3, name)
	dimension_reduction(data, labels, 'ISOMAP', show_figures, 3, name)
	dimension_reduction(data, labels, 't-SNE', show_figures, 3, name)
	dimension_reduction(data, labels, 'LDA', show_figures, 1, name)
	dimension_reduction(data, labels, 'QDA', show_figures, 1, name)


current_data, current_labels, current_fixation, current_indexes = dh.get_subset('unique', dh.get_measured_spectra(), dh.get_labels())
print('Subset contains ', len(current_labels), ' observations')

compare_dimension_reduction(current_data, current_labels, True, 'measured')
