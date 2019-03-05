from os.path import dirname, join as pjoin
import scipy.io as sio 

base_dir  = '/media/sf_research/input/'
data_dir = 'saitama_v2_min_square'
mat_fname = pjoin(base_dir, data_dir, 'in.mat')
mat_contents = sio.loadmat(mat_fname)

#from sklearn.model_selection import train_test_split
#X_train, X_test, y_train, y_test = train_test_split(X,y, test_size=0.5)