import numpy as np
from matplotlib import pyplot as plt
from sklearn.neural_network import MLPRegressor

np.random.seed(0)

npt = 20

x = np.loadtxt("pair_pt_drellyan.txt", usecols=0, skiprows=0, dtype='float')
y = np.loadtxt("pair_pt_drellyan.txt", usecols=1, skiprows=0, dtype='float')
err_y = np.loadtxt("pair_pt_drellyan.txt", usecols=2, skiprows=0, dtype='float')

plt.figure()
plt.errorbar(x.T, y, yerr = err_y, fmt="r.", markersize=10)

val = []
score = []

count = 0

while count < 200:

    X = []
    Y = []

    y1 = np.zeros(npt)

    for b in range(npt):
        while y1[b] <= 0:
            y1[b] = np.random.normal(y[b], err_y[b])

    X = np.append(X, x)
    X = np.atleast_2d(X).T

    Y = np.append(Y, np.log10(y1))

    regr = MLPRegressor(max_iter=10000, solver='lbfgs', activation='logistic').fit(X, Y) #'relu' logistic softmax tanh
    model_score = regr.score(np.atleast_2d(x).T, np.log10(y))

    if model_score < 0.70: #This matters a lot
        continue
    
    count = count + 1

    xpred = np.atleast_2d(np.linspace(0.0, 5.0, 200)).T
    ypred = regr.predict(xpred)
    
    val.append(ypred)
    score.append(model_score)
    
    if count%50 == 0:
        print(count)

value = np.power(10,val)
average = np.average(value, weights=score, axis=0)
sig = np.std(value, axis=0)

x_coor = xpred
y_coor = np.atleast_2d(average).T
e_coor = np.atleast_2d(sig).T

myarray = np.concatenate((x_coor,y_coor,e_coor),axis=1)
np.savetxt("mlp_pair_pt_drellyan.txt", myarray, fmt='%.5e', delimiter='\t', newline='\n', header='', footer='', comments='# ', encoding=None)

plt.yscale('log')
#plt.ylim(1e-13, 1e-8)
plt.plot(xpred, average)

plt.fill(
    np.concatenate([x_coor, x_coor[::-1]]),
    np.concatenate([y_coor - 1.00 * e_coor, (y_coor + 1.00 * e_coor)[::-1]]),
    alpha=0.5,
    fc="r",
    ec="None",
    label="1 sigma",
)

plt.show()
