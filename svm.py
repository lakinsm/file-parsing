'''
@author: Steven Lakin
Date Created: October 31st, 2015
Following along with the sklearn examples for using different kernels with SVM classifiers
'''

import numpy as np
import matplotlib.pyplot as plt
from sklearn import svm, datasets

# Load the Iris data from the sklearn datasets
iris = datasets.load_iris()
X, y = iris.data, iris.target

# We will look only at a binary classification, so we need to pull out the data and features
# for two of the three classes
X = X[y != 0, :2]
y = y[y != 0]

# How many samples are there?
n_sample = len(X)

# Randomize sample order
np.random.seed(54)
order = np.random.permutation(n_sample)
X = X[order]
y = y[order].astype(np.float)

# Put the data into train or test categories
X_train = X[:int(.9 * n_sample)]
y_train = y[:int(.9 * n_sample)]
X_test = X[int(.9 * n_sample):]
y_test = y[int(.9 * n_sample):]

# Fit the model using three different kernels
for fig_num, kernel in enumerate(('linear', 'rbf', 'poly')):
    clf = svm.SVC(kernel=kernel, gamma=10)
    clf.fit(X_train, y_train)

    plt.figure(fig_num)
    plt.clf()
    plt.scatter(X[:, 0], X[:, 1], c=y, zorder=10, cmap=plt.cm.Paired)

    # Circle out the test data
    plt.scatter(X_test[:, 0], X_test[:, 1], s=80, facecolors='none', zorder=10)

    plt.axis('tight')
    x_min = X[:, 0].min()
    x_max = X[:, 0].max()
    y_min = X[:, 1].min()
    y_max = X[:, 1].max()

    XX, YY = np.mgrid[x_min:x_max:200j, y_min:y_max:200j]
    Z = clf.decision_function(np.c_[XX.ravel(), YY.ravel()])

    # Put the result into a color plot
    Z = Z.reshape(XX.shape)
    plt.pcolormesh(XX, YY, Z > 0, cmap=plt.cm.Paired)
    plt.contour(XX, YY, Z, colors=['k', 'k', 'k'], linestyles=['--', '-', '--'],
                levels=[-.5, 0, .5])

    plt.title(kernel)
plt.show()





