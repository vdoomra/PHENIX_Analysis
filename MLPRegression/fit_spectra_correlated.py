import numpy as np
from matplotlib import pyplot as plt
from sklearn.neural_network import MLPRegressor

np.random.seed(0)

npt = 13

x = np.loadtxt("outFile_DCA_Fits_BG.txt", usecols=0, skiprows=0, dtype='float')
y = np.loadtxt("outFile_DCA_Fits_BG.txt", usecols=1, skiprows=0, dtype='float')
err_y = np.loadtxt("outFile_DCA_Fits_BG.txt", usecols=2, skiprows=0, dtype='float')

plt.figure()
plt.errorbar(x.T, y, yerr = err_y, fmt="r.", markersize=10)

val = []
score = []

count = 0

while count < 500:

    shared_fluctuation = np.random.normal(0, 1)
    shifts = shared_fluctuation*err_y
            
    y1 = y+shifts

    if np.any(y1 <= 0):
        continue

    X = np.atleast_2d(x).T
    Y = np.log10(y1)

    regr = MLPRegressor(max_iter=5000, solver='lbfgs', activation='tanh').fit(X, Y) #'relu' logistic softmax tanh
    model_score = regr.score(np.atleast_2d(x).T, np.log10(y))

    if model_score < 0.60: #This matters a lot
        continue
    
    count = count + 1

    xpred = np.atleast_2d(np.linspace(0.9, 4.5, 200)).T
    ypred = 10**regr.predict(xpred)
    
    val.append(ypred)
    score.append(model_score)
    
    if count%50 == 0:
        print(count)

average = np.average(val, weights=score, axis=0)
sig = np.std(val, axis=0)

x_coor = xpred
y_coor = np.atleast_2d(average).T
e_coor = np.atleast_2d(sig).T

myarray = np.concatenate((x_coor,y_coor,e_coor),axis=1)
np.savetxt("mlp_DCA_Fits_BG.txt", myarray, fmt='%.5e', delimiter='\t', newline='\n', header='', footer='', comments='# ', encoding=None)

plt.yscale('log')
plt.ylim(1e-11, 1e-7)
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
