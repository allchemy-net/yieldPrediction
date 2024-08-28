import statistics
import numpy
from sklearn.linear_model import LinearRegression
from sklearn.feature_selection import r_regression


def getFinalProductConc(masses, args):
    return max(masses[0]), max(masses[1]), max(masses[2])


def getPreds(res, names, args):
    pred0 = [max(res[name][0]) for name in names]
    pred1 = [max(res[name][1]) for name in names]
    pred2 = [max(res[name][2]) for name in names]
    return pred0, pred1, pred2


def getPearson(res, expYields, args):
    names = tuple(sorted(expYields.keys()))
    names = tuple(sorted(res.keys()))
    full = dict()
    if isinstance(res[names[0]][0], (tuple, list)):
        yld = [expYields[name] for name in names]
        pred0, pred1, pred2 = getPreds(res, names, args)
        r2 = []
        full = {name: (yld[poz], pred0[poz], pred1[poz], pred2[poz]) for poz, name in enumerate(names)}
    else:
        # --propagationMode front
        raise NotImplementedError
    Y = numpy.array(yld)
    # Y = Y.reshape(-1, 1)
    for x in (pred0, pred1, pred2):
        X = numpy.array(x)
        X = X.reshape(-1, 1)
        # print("SHAPE", X.shape, Y.shape)
        reg = r_regression(X, Y, force_finite=True)[0]
        print("PEARSON", type(reg), reg)
        r2.append(reg)
    return r2, full


def getR2_old(res, expYields, args):
    names = tuple(sorted(expYields.keys()))
    names = tuple(sorted(res.keys()))
    full = dict()
    if isinstance(res[names[0]][0], (tuple, list)):
        yld = [expYields[name] for name in names]
        pred0, pred1, pred2 = getPreds(res, names, args)
        r2 = []
        full = {name: (yld[poz], pred0[poz], pred1[poz], pred2[poz]) for poz, name in enumerate(names)}
    else:
        # --propagationMode front
        raise NotImplementedError
    for x in (pred0, pred1, pred2):
        # slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, yld)
        reg = LinearRegression().fit(x, yld)
        r2.append(reg.score(x, yld))
        # r2.append(r_value**2)
    return r2, full


def getR2(res, expYields, args):
    return _getR2_linreg(res, expYields, args, intercept=True)


def getR2_0(res, expYields, args):
    return _getR2_linreg(res, expYields, args, intercept=False)


def _getR2_linreg(res, expYields, args, intercept):
    names = tuple(sorted(expYields.keys()))
    names = tuple(sorted(res.keys()))
    full = dict()
    if isinstance(res[names[0]][0], (tuple, list)):
        yld = [expYields[name] for name in names]
        pred0, pred1, pred2 = getPreds(res, names, args)
        r2 = []
        full = {name: (yld[poz], pred0[poz], pred1[poz], pred2[poz]) for poz, name in enumerate(names)}
    else:
        # --propagationMode front
        raise NotImplementedError
    for x in (pred0, pred1, pred2):
        # slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, yld)
        x = numpy.array(x).reshape((-1, 1))
        reg = LinearRegression(fit_intercept=intercept).fit(x, yld)
        r2.append(reg.score(x, yld))
        # r2.append(r_value**2)
    return r2, full


def getMAE(res, expYields, args):
    names = tuple(sorted(expYields.keys()))
    names = tuple(sorted(res.keys()))
    full = dict()
    if isinstance(res[names[0]][0], (tuple, list)):
        yld = [expYields[name] for name in names]
        pred0, pred1, pred2 = getPreds(res, names, args)
        MAEs = []
        full = {name: (yld[poz], pred0[poz], pred1[poz], pred2[poz]) for poz, name in enumerate(names)}
    else:
        # --propagationMode front
        raise NotImplementedError
    for x in (pred0, pred1, pred2):
        # slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, yld)
        maes = statistics.mean([abs(xi - yi) for xi, yi in zip(x, yld)])
        # r2.append(r_value**2)
        MAEs.append(maes)
    return MAEs, full


def getMSE(res, expYields, args):
    names = tuple(sorted(expYields.keys()))
    names = tuple(sorted(res.keys()))
    full = dict()
    if isinstance(res[names[0]][0], (tuple, list)):
        yld = [expYields[name] for name in names]
        pred0, pred1, pred2 = getPreds(res, names, args)
        MSEs = []
        full = {name: (yld[poz], pred0[poz], pred1[poz], pred2[poz]) for poz, name in enumerate(names)}
    else:
        # --propagationMode front
        raise NotImplementedError
    for x in (pred0, pred1, pred2):
        # slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, yld)
        mses = statistics.mean([(xi - yi) ** 2 for xi, yi in zip(x, yld)])
        # r2.append(r_value**2)
        MSEs.append(mses)
    return MSEs, full


def getRMSE(res, expYields, args):
    names = tuple(sorted(expYields.keys()))
    names = tuple(sorted(res.keys()))
    full = dict()
    if isinstance(res[names[0]][0], (tuple, list)):
        yld = [expYields[name] for name in names]
        pred0, pred1, pred2 = getPreds(res, names, args)
        MAEs = []
        full = {name: (yld[poz], pred0[poz], pred1[poz], pred2[poz]) for poz, name in enumerate(names)}
    else:
        # --propagationMode front
        raise NotImplementedError
    for x in (pred0, pred1, pred2):
        # slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, yld)
        maes = statistics.mean([(xi - yi) ** 2 for xi, yi in zip(x, yld)])
        # r2.append(r_value**2)
        MAEs.append(maes)
    return MAEs, full
