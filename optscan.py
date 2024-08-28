import argparse, json, traceback
from functools import partial
from multiprocessing import Queue, Process
import helpers, consts, stats, statsMetrics


# ###############################################################


def loadOptSpaceFromFile(scanfile):
    fh = open(scanfile)
    optScape = json.load(fh)
    fh.close()
    return optScape


def getArgsParamFromSpaceSampling(args, optSpace):
    newargs = argparse.Namespace(**vars(args))
    spaceOrdered = {par: getScanOrderValues(optSpace[par]) for par in optSpace}
    # newargs = setValues


def getScanOrderValues(optSpacePar):
    if 'min' in optSpacePar and 'max' in optSpacePar and 'default' in optSpacePar and 'step' in optSpacePar:
        if isinstance(optSpacePar['default'], (tuple, list)):
            allValues = getAllScanValuesForManyD(optSpacePar)
        else:
            allValues = getAllScanValuesFor1D(optSpacePar)
    else:
        raise NotImplementedError
    return allValues


def getAllScanValuesForManyD(optSpacePar):
    lists = []
    for poz in range(len(optSpacePar['default'])):
        dct = {key: optSpacePar[key][poz] for key in optSpacePar}
        lists.append(getAllScanValuesFor1D(dct))
    return list(zip(*lists))


def getAllScanValuesFor1D(optSpacePar):
    allValues = [optSpacePar['default'], ]
    allAbove = []
    itr = 1
    while True:
        newVal = optSpacePar['default'] + (itr * optSpacePar['step'])
        if newVal <= optSpacePar['max']:
            allAbove.append(newVal)
            itr += 1
        else:
            break
    allBelow = []
    itr = 1
    while True:
        newVal = optSpacePar['default'] - (itr * optSpacePar['step'])
        if newVal >= optSpacePar['max']:
            allBelow.append(newVal)
            itr += 1
        else:
            break
    lenAbove = len(allAbove)
    lenBelow = len(allBelow)
    mlen = max(lenAbove, lenBelow)
    aboveIdx = getMidIdx(lenAbove)
    belowIdx = getMidIdx(lenBelow)
    for i in range(mlen):
        if i < lenAbove:
            idx = aboveIdx[i]
            allValues.append(allAbove[idx])
        if i < lenBelow:
            idx = belowIdx[i]
            allValues.append(allBelow[idx])
    return allValues


def getMidIdx(numelem):
    allList = [0, numelem - 1, ]
    allSet = set(allList)
    for isplit in range(1, int(numelem**0.5)):
        elems = [(numelem * i) / (2 * isplit) for i in range(1, 1 + isplit)]
        for el in elems:
            if el in allSet:
                continue
            allList.append(el)
            allSet.add(el)
    for i in range(numelem):
        if i in allSet:
            continue
        allList.append(i)
    return allList


def performScan(dataDict, initStoich, adbData, args, optSpace):

    def _calcProcess(qin, qout):
        while True:
            msg = qin.get()
            if msg == 'STOP':
                break
            argsPar = msg
            dataToPrint = stats.getStats(dataDict, initStoich, adbData, argsPar)
            resultSummary = stats.doSummaryAndPprint(dataToPrint, consts.expYields, argsPar)
            # return {'r2': r2, 'results': results}
            print("RES", resultSummary)
            qout.put(resultSummary)
    qin, qout = Queue(), Queue()
    procs = [Process(target=_calcProcess, args=(qin, qout)) for _ in range(args.numprocs)]
    _ = [p.start() for p in procs]
    for argsPar in getArgsParamFromSpaceSampling(args, optSpace):
        qin.put(argsPar)
    space = None
    return space

# ###############################################################
# hyper opt


def optFunctionOrion(dataDict, initStoich, adbData, args, **kwargs):
    args = helpers.replaceArgFromDict(args, kwargs)
    try:
        res = stats.getStatsSingle(dataDict, initStoich, adbData, args)
        if args.paropt == 'R2':
            r2, full = statsMetrics.getR2(res, consts.expYields, args)
        else:
            raise NotImplementedError
        print("R2:", r2, res, args)
        r2_0, r2_1, r2_2 = r2
        score = 1 - r2_2
    except consts.NegativeMassException:
        print("NEG MASS")
        score = 1
    except:
        print("FAILED", kwargs, traceback.format_exc())
        score = 1
    return [{"name": "objective", "type": "objective", "value": score}]


def buildOrionSpaceFromDict(optSpace):
    retDict = dict()
    for param in optSpace:
        if 'min' in optSpace[param]:
            minim = optSpace[param]['min']
            maxim = optSpace[param]['max']
            setInt = False
            if isinstance(minim, int) and isinstance(maxim, int):
                setInt = True
            if isinstance(minim, (tuple, list)):
                for poz in range(len(minim)):
                    if setInt:
                        retDict[f"{param}__PART{poz}_{len(minim)}_MULTIPART"] = f"uniform({minim[poz]}, {maxim[poz]}, discrete=True)"
                    else:
                        retDict[f"{param}__PART{poz}_{len(minim)}_MULTIPART"] = f"uniform({minim[poz]}, {maxim[poz]})"
            else:
                if setInt:
                    retDict[param] = f"uniform({minim}, {maxim}, discrete=True)"
                else:
                    retDict[param] = f"uniform({minim}, {maxim})"
        elif 'values':
            retDict[param] = f"choices{optSpace[param]['values']}"
        else:
            raise NotImplementedError
    return retDict


def orionOpt(fn, optSpace, args):
    from orion.client import build_experiment
    storage = {"type": "legacy", "database": {"type": "pickleddb", "host": "./db.pkl"}}
    # space = {"x": "uniform(0, 30)"}
    space = buildOrionSpaceFromDict(optSpace)
    print("SPACE", space)
    algo = {"tpe": {"n_initial_points": 5}}
    experiment = build_experiment("random-rosenbrock", space=space, storage=storage, algorithm=algo, debug=True)
    experiment.workon(fn, max_trials=25_000)  # , n_workers=args.numprocs, )


# ## openbox #####################################


def _getMaxValue(r2, args):
    lst = _getValuesToFit(r2, args)
    return max(lst)


def _getMinValue(r2, args):
    lst = _getValuesToFit(r2, args)
    return min(lst)


def _getValuesToFit(r2, args):
    retList = []
    if '0' in args.fitTo:
        retList.append(r2[0])
    if '1' in args.fitTo:
        retList.append(r2[1])
    if '2' in args.fitTo:
        retList.append(r2[2])
    if not retList:
        raise NotImplementedError
    return retList


def optFunctionOpenBox(dataDict, initStoich, adbData, args, kwargs):
    args = helpers.replaceArgFromDict(args, kwargs)
    try:
        res = stats.getStats(dataDict, initStoich, adbData, args)
        if not res:
            score = 1
        else:
            if args.paropt == 'R2':
                r2, full = statsMetrics.getR2(res, consts.expYields, args)
            elif args.paropt == 'R2_0':
                r2, full = statsMetrics.getR2_0(res, consts.expYields, args)
            elif args.paropt == 'MAE':
                r2, full = statsMetrics.getMAE(res, consts.expYields, args)
            elif args.paropt == 'MSE':
                r2, full = statsMetrics.getMSE(res, consts.expYields, args)
            elif args.paropt == 'RMSE':
                r2, full = statsMetrics.getRMSE(res, consts.expYields, args)
            elif args.paropt == 'pearson':
                r2, full = statsMetrics.getPearson(res, consts.expYields, args)
            else:
                raise NotImplementedError
            print("R2::", *r2, kwargs, full, args)
            # r2_0, r2_1, r2_2 = r2
            if args.paropt in {'R2', 'R2_0', 'pearson'}:
                score = 1 - _getMaxValue(r2, args)
            elif args.paropt in {'MAE', 'MSE', 'RMSE'}:
                score = _getMinValue(r2, args)
            else:
                raise NotImplementedError
    except consts.NegativeMassException:
        print("NEG MASS")
        if args.paropt in {'R2', 'R2_0', 'pearson'}:
            score = 1
        elif args.paropt in {'MAE', 'RMSE'}:
            score = 100
        elif args.paropt in {'MSE'}:
            score = 10_000  # 100 * 100 max error for MSE
        else:
            raise NotImplementedError

    except:
        print("FAILED", kwargs, traceback.format_exc())
        if args.paropt in {'R2', 'pearson', 'R2_0'}:
            score = 1
        elif args.paropt in {'MAE', 'RMSE'}:
            score = 100
        elif args.paropt in {'MSE', }:
            score = 10_000
        else:
            raise NotImplementedError
    return {"objectives": score}


def buildOpenBoxSpaceFromDict(optSpace, limits):
    from openbox import space as sp
    # space = sp.Space()
    space = sp.ConditionedSpace()
    # space.add_variables([x1, x2])

    # this function define impossible combination of paramethers i.e
    def sample_condition(config):
        # use limits
        if 'scaleByADBspeed__PART0_3_MULTIPART' in config:
            if config['scaleByADBspeed__PART0_3_MULTIPART'] < config['scaleByADBspeed__PART1_3_MULTIPART']:
                return False
            if config['scaleByADBspeed__PART1_3_MULTIPART'] < config['scaleByADBspeed__PART2_3_MULTIPART']:
                return False
        if 'sidemechYield' in config:
            if limits.get('noSideMech', False):
                if config['sidemechYield'] != config['defaultYield']:
                    return False
            else:
                if config['sidemechYield'] > 0.5 * config['defaultYield'] or config['fromTargetYield'] > config['defaultYield']:
                    return False
        if 'fixAcidBaseRatio' in config:
            if config['fixAcidBaseRatio'] == 'Y' and config['addReverse'] == 'wildToMainOnly':
                return False
            if config['fixAcidBaseRatio'] == 'N' and config['setAcidBaseRatio'] != 1.0:
                # allow one value to prevent scaning when it cannot change anything
                return False
        if 'maxAllowDEtr' in config:
            if config['maxAllowDEtr'] >= config['maxAllowDE']:
                return False
        if 'forwardReversableScaler' in config and 'reverseSpeed' in config:
            if config['forwardReversableScaler'] / config['reverseSpeed'] > 5:
                return False
        return True

    retList = []
    for param in optSpace:
        if 'min' in optSpace[param]:
            minim = optSpace[param]['min']
            maxim = optSpace[param]['max']
            default = optSpace[param]['default']
            setInt = False
            if isinstance(minim, int) and isinstance(maxim, int):
                setInt = True
            if isinstance(minim, (tuple, list)):
                for poz in range(len(minim)):
                    if setInt:
                        variable = sp.Int(f"{param}__PART{poz}_{len(minim)}_MULTIPART", minim[poz], maxim[poz], default_value=default[poz])
                    else:
                        variable = sp.Real(f"{param}__PART{poz}_{len(minim)}_MULTIPART", minim[poz], maxim[poz], default_value=default[poz])
                    retList.append(variable)
            else:
                if setInt:
                    variable = sp.Int(param, minim, maxim, default_value=default)
                else:
                    variable = sp.Real(param, minim, maxim, default_value=default)
                retList.append(variable)
        elif 'values' in optSpace[param]:
            retList.append(sp.Categorical(param, optSpace[param]['values'], default_value=optSpace[param]['default']))
        else:
            raise NotImplementedError
    space.add_variables(retList)
    space.set_sample_condition(sample_condition)
    return space


def openBoxOpt(dataDict, initStoich, adbData, optSpace, args):
    # optFn = partial(optFunctionOpenBox, dataDict, initStoich, adbData, args)
    # doParOpt(dataDict, args, params, expYields, initStoich)
    print("OPT")
    optFn = partial(optFunctionOpenBox, dataDict, initStoich, adbData, args)
    # optscan.openBoxOpt(optFn, optSpace, args)
    # optscan.openBoxOpt(dataDict, initStoich, adbData, arg
    return _openBoxOpt(optFn, optSpace, args)


def _openBoxOpt(fn, optSpace, args):
    from openbox import Optimizer, logger
    # Run
    limits = getLimitsFromArgs(args)
    space = buildOpenBoxSpaceFromDict(optSpace, limits)
    loggerSetting = {'name': 'OpenBox', 'stream': False,
                     'level': 'DEBUG',  # ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
                     'logdir': 'logs/',
                     'fmt': logger.DEFAULT_FORMAT,
                     'color': True}
    # prf gp
    # print("FN", fn, space)
    opt = Optimizer(fn, space, max_runs=args.optSteps,  surrogate_type=args.openboxSurrogate, task_id='q_gp4', logger_kwargs=loggerSetting)
    # Have a try on the new HTML visualization feature!
    # visualization='advanced',   # or 'basic'. For 'advanced', run 'pip install "openbox[extra]"' first
    # auto_open_html=True,        # open the visualization page in your browser automatically
    history = opt.run()
    print("HISTORY", history)
    return history


def openBoxOptParal(fn, optSpace, args):
    from openbox import ParallelOptimizer
    limits = getLimitsFromArgs(args)
    space = buildOpenBoxSpaceFromDict(optSpace, limits)
    opt = ParallelOptimizer(fn, space, parallel_strategy='async', batch_size=4, batch_strategy='default', num_objectives=1,
                            num_constraints=0, max_runs=50, surrogate_type='gp', task_id='parallel_async')
    history = opt.run()
    print(history)


def getLimitsFromArgs(args):
    retDict = dict()
    if args.disableOpt == 'noSpeedClass':
        retDict['noSpeedClass'] = True
    elif args.disableOpt == 'noSideMech':
        retDict['noSideMech'] = True
    return retDict


def adjustOptSpace(optSpace, args):
    if args.useEvansPolanyi == 'N':
        print("dont use Evans-Polanyi contribution !!!")
        DEopt = {'scaleByDE', 'maxAllowDEtr', 'maxAllowDE'}
        for opt in DEopt:
            if opt in optSpace:
                del optSpace[opt]
                print("REMOVE FROM OPT", opt)
    if args.useMayr == 'N':
        mayrOpt = {'scaleByMayr', }
        for opt in mayrOpt:
            if opt in optSpace:
                del optSpace[opt]
                print("REMOVE FROM OPT", opt)
    return optSpace
