"""
Working with the Kaggle face key data (kfkd)
"""

import os
import numpy as np
from pandas.io.parsers import read_csv
from sklearn.utils import shuffle
from lasagne import layers
from lasagne.updates import nesterov_momentum
from nolearn.lasagne import NeuralNet

FTRAIN = "~/Documents/phd/phdenv/kaggle_face_tutorial/training.csv"
FTEST = "~/Documents/phd/phdenv/kaggle_face_tutorial/test.csv"

def load(test=False, cols=None):
    """
    Loads data from FTEST if 'test' is True, otherwise from FTRAIN.
    Pass a list of 'cols' if you're only interested in a subset of the target columns.
    """
    fname = FTEST if test else FTRAIN
    df = read_csv(os.path.expanduser(fname))  # Load pandas dataframe

    # The Image column has pixel values separated by a space; convert
    # the values to numpy arrays:
    df['Image'] = df['Image'].apply(lambda im: np.fromstring(im, sep=' '))

    if cols:  # Get a subset of columns
        df = df[list(cols) + ['Image']]

    print(df.count())  # Prints the number of values for each column
    df = df.dropna()  # Drop all rows that have missing values

    X = np.vstack(df['Image'].values) / 255.  # Scale pixel values to [0, 1]
    X = X.astype(np.float32)

    if not test:  # Only FTRAIN has target columns
        y = df[df.columns[:-1]].values
        y = (y - 48) / 48  # Scale target coordinates to [-1, 1]
        X, y = shuffle(X, y, random_state=42)  # Shuffle the train data
        y = y.astype(np.float32)
    else:
        y = None

    return X, y


def load2d(test=False, cols=None):
    X, y = load(test=test)
    X = X.reshape(-1, 1, 96, 96)
    return X, y


# Basic single hidden layer neural net, no normalization, augmentation, optimization, or dropout
net1 = NeuralNet(
    layers=[  # Three layers, one hidden layer
              ('input', layers.InputLayer),
              ('hidden', layers.DenseLayer),
              ('output', layers.DenseLayer),
    ],
    # Layer parameters
    input_shape=(None, 9216),  # 96x96 input pixels per batch
    hidden_num_units=100,  # number of units in hidden layer
    output_nonlinearity=None,  # output layer uses identity function
    output_num_units=30,  # 30 target values

    # Optimization model
    update=nesterov_momentum,
    update_learning_rate=0.01,
    update_momentum=0.9,

    regression=True,  # flag to indicate this is a regression problem
    max_epochs=400,  # train this many epochs
    verbose=1,
)


net2 = NeuralNet(
    layers=[
        ('input', layers.InputLayer),
        ('conv1', layers.Conv2DLayer),
        ('pool1', layers.MaxPool2DLayer),
        ('conv2', layers.Conv2DLayer),
        ('pool2', layers.MaxPool2DLayer),
        ('conv3', layers.Conv2DLayer),
        ('pool3', layers.MaxPool2DLayer),
        ('hidden4', layers.DenseLayer),
        ('hidden5', layers.DenseLayer),
        ('output', layers.DenseLayer),
    ],
    input_shape=(None, 1, 96, 96),
    conv1_num_filters=32, conv1_filter_size=(3, 3), pool1_pool_size=(2, 2),
    conv2_num_filters=64, conv2_filter_size=(2, 2), pool2_pool_size=(2, 2),
    conv3_num_filters=128, conv3_filter_size=(2, 2), pool3_pool_size=(2, 2),
    hidden4_num_units=500, hidden5_num_units=500,
    output_num_units=30, output_nonlinearity=None,

    update_learning_rate=0.01,
    update_momentum=0.9,

    regression=True,
    max_epochs=1000,
    verbose=1,
)

X, y = load2d()
# print("X.shape == {}; X.min == {:.3f}; X.max == {:.3f}".format(X.shape, X.min(), X.max()))
# print("y.shape == {}; y.min == {:.3f}; y.max == {:.3f}".format(y.shape, y.min(), y.max()))

net2.fit(X, y)

import cPickle as pickle
with open('net2.pickle', 'wb') as f:
    pickle.dump(net2, f, -1)





