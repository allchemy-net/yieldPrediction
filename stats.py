import traceback, itertools
from multiprocessing import Queue, Process
import numpy
from sklearn.linear_model import LinearRegression
# from sklearn.feature_selection import r_regression
import calc, consts, statsMetrics


def getStatService(dataDict, initStoich, adbData, args):
    sbsStoich = None
    print("KEY", initStoich.keys())
    if initStoich and initStoich['cmd1']:
        sbsStoich = initStoich['cmd1']
    if args.debug:
        print("====", name, key)
    try:
        gen = [g for g in dataDict['cmd1'].keys() if g != 0]
        if len(gen) != 1:
            raise NotImplementedError
        gen = gen[0]
        masses, fullMass, rxChanged, target = calc.propagateNetwork(dataDict['cmd1'][gen], sbsStoich, adbData, dataDict['cmd1'][0]['common'], args)
        retDict = {'target': target, 'masses': masses, 'fullMass': fullMass, 'rxChanged': rxChanged}
        keyToCopy = ('targetTautomers', 'targetWildForms', 'targetTautomersWild')
        for keyToCp in keyToCopy:
            if keyToCp in dataDict['cmd1'][gen]['data']:
                retDict[keyToCp] = dataDict['cmd1'][gen]['data'][keyToCp]
    except consts.NegativeMassException:
        print("NEGATIVE MASS in calc")
        print(f"CATCHED TB: {traceback.format_exc()}")
        retDict = {'inp': ('ERROR', 'ERROR', 'ERROR')}
    except:
        print("OTHER ERROR in calc")
        print(f"CATCHED TB: {traceback.format_exc()}")
        retDict = {'inp': ('ERROR', 'ERROR', 'ERROR')}
    return retDict


def getStats(dataDict, initStoich, adbData, args):
    def _calcSingle(qin, qout):
        while True:
            msg = qin.get()
            if msg == 'STOP':
                return
            molname, name, key = msg
            sbsStoich = None
            # print("KEY", initStoich.keys())
            if initStoich and initStoich[molname]:
                sbsStoich = initStoich[molname]
            if args.debug:
                print("====", name, key)
            try:
                masses, fullMass, rxChanged, target = calc.propagateNetwork(dataDict[name][key], sbsStoich, adbData, dataDict[name][0]['common'], args)
                retDict = {'inp': msg, 'target': target, 'masses': masses, 'fullMass': fullMass, 'rxChanged': rxChanged}
                keyToCopy = ('targetTautomers', 'targetWildForms', 'targetTautomersWild')
                for keyToCp in keyToCopy:
                    if keyToCp in dataDict[name][key]['data']:
                        retDict[keyToCp] = dataDict[name][key]['data'][keyToCp]
            except consts.NegativeMassException:
                print("NEGATIVE MASS in ", msg)
                print(f"CATCHED TB: {traceback.format_exc()}")
                retDict = {'inp': ('ERROR', 'ERROR', 'ERROR')}
            except:
                print("OTHER ERROR", msg)
                print(f"CATCHED TB: {traceback.format_exc()}")
                retDict = {'inp': ('ERROR', 'ERROR', 'ERROR')}
            qout.put(retDict)
    ##########

    qin, qout = Queue(), Queue()
    maxProc = args.numprocs  # check how many point to calc and adjust max proc
    print(maxProc, 'process(es) will be started')
    procs = [Process(target=_calcSingle, args=(qin, qout)) for _ in range(maxProc)]
    _ = [p.start() for p in procs]
    dataToPrint = dict()
    msgSubmitted = 0
    for name in dataDict:
        molname = name.split('_')[0]
        for key in dataDict[name]:
            msg = (molname, name, key)
            qin.put(msg)
            msgSubmitted += 1
    for _ in procs:
        qin.put('STOP')
    qin.close()
    # isError = False
    fh = None
    fhrx = None
    for idx in range(msgSubmitted):
        retDict = qout.get()
        molname, name, key = retDict['inp']
        if args.saveDataForPlot:
            fh = open(f'{args.saveDataForPlot}_{molname}_{key}_cmds.tsv', 'w')
            fhrx = open(f'{args.saveDataForPlot}_{molname}_{key}_rxes.tsv', 'w')
        if molname == 'ERROR':
            print("ERRROR??", idx)
            # isError = True
            return None
        target = retDict['target']
        if args.debug:
            print("====", name, key)
        # masses, fullMass, rxChanged, target = propagateNetwork(dataDict[name][key], sbsStoich, adbData, dataDict[name][0]['common'], args)
        if fh:
            print(name, key, formFullCmdMass(retDict['fullMass']), sep='\t', file=fh)
            if args.debug > 100:
                dcts = dict()
                for dct in retDict['fullMass']:
                    for k in dct:
                        if k not in dcts:
                            dcts[k] = 1
                        else:
                            dcts[k] += 1
                for k in retDict['fullMass'][-1]:
                    dct[k] *= -1
                print(f"fullMass:76: {name} {key} :: {len(retDict['fullMass'])} :: {dcts}")
            print(name, key, formFullRxMass(retDict['rxChanged'], dataDict[name][key]['data']), sep='\t', file=fhrx)
            if args.debug > 100:
                print(f"fullrx:79: {name} {key} :: {[(x, len(retDict['rxChanged'][x])) for x in retDict['rxChanged']]}")
        if args.debug:
            massDct = {smi: retDict['masses'][smi] if smi != target else max(retDict['masses'][smi]) for smi in retDict['masses']}
            print("MASSES", name, key, target, massDct)
        if args.pprint or args.paropt != 'N':
            if molname not in dataToPrint:
                dataToPrint[molname] = dict()
            # ============tutaj dodac dodawania tautomerow targetu
            dataToPrint[molname][key] = _getTargetMasses(retDict, target, args)
        else:
            targetMasses = _getTargetMasses(retDict, target, args)
            print(name, key, ">>>>>", targetMasses, len(retDict['masses']))
        if fh:
            fh.close()
    qout.close()
    return dataToPrint


def _getTargetMasses(retDict, target, args):
    summed = []
    if 'targetTautomers' in retDict:
        # print("TT", target, retDict['targetTautomers'], "MODE", args.yieldDefMode)
        if args.yieldDefMode in ('maxWithTauto', 'lastWithTauto', 'lastWithTautoAndDeProt', 'maxWithTautoAndDeProt'):
            numSteps = len(retDict['masses'][target])
            lenDct = len(retDict['fullMass'])
            if not lenDct:
                print("XXXX", retDict['targetTautomers'])
                return retDict['masses'][target]
            lastVal = dict()
            # print("NS", numSteps)
            for step in range(numSteps):
                if step < lenDct:
                    dct = retDict['fullMass'][step]
                else:
                    # print("LEE", lenDct - 1, len(retDict['fullMass']))
                    dct = retDict['fullMass'][lenDct - 1]
                try:
                    vals = [dct[smi] for smi in retDict['targetTautomers']]
                    lastVal = vals
                except KeyError:
                    # if no reaction (no value) at this step try to get concentration from prev step i.e. lastVal
                    # if no prev step set, so no reation done and concentration is 0
                    vals = [dct.get(smi, lastVal.get(smi, 0)) for smi in retDict['targetTautomers']]
                    lstVal = vals
                # print("NS", step, numSteps, sum(vals))
                summed.append(sum(vals))
    # print(len(summed), max(summed))
    if 'targetWildForms' in retDict:
        appendMode = True
        if summed:
            appendMode = False
        if args.yieldDefMode in ('lastWithTautoAndDeProt', 'maxWithTautoAndDeProt'):
            numSteps = len(retDict['masses'][target])
            lenDct = len(retDict['fullMass'])
            if not lenDct:
                print("XXXX", retDict['targetTautomers'])
                return retDict['masses'][target]
            lastVal = dict()
            for step in range(numSteps):
                if step < lenDct:
                    dct = retDict['fullMass'][step]
                else:
                    # print("LEE", lenDct - 1, len(retDict['fullMass']))
                    dct = retDict['fullMass'][lenDct - 1]
                try:
                    vals = [dct[smi] for smi in retDict['targetWildForms']]
                    lastVal = vals
                except KeyError:
                    # if no reaction (no value) at this step try to get concentration from prev step i.e. lastVal
                    # if no prev step set, so no reation done and concentration is 0
                    vals = [dct.get(smi, lastVal.get(smi, 0)) for smi in retDict['targetWildForms']]
                    lstVal = vals
                if appendMode:
                    # no tautomer summation so add also main form concentration
                    summed.append(sum(vals) + retDict['masses'][target][step])
                else:
                    summed[step] += sum(vals)
    if 'targetTautomersWild' in retDict:
        appendMode = True
        if summed:
            appendMode = False
        if args.yieldDefMode in ('lastWithTautoAndDeProt', 'maxWithTautoAndDeProt'):
            numSteps = len(retDict['masses'][target])
            lenDct = len(retDict['fullMass'])
            if not lenDct:
                print("XXXX", retDict['targetTautomersWild'])
                return retDict['masses'][target]
            lastVal = dict()
            for step in range(numSteps):
                if step < lenDct:
                    dct = retDict['fullMass'][step]
                else:
                    # print("LEE", lenDct - 1, len(retDict['fullMass']))
                    dct = retDict['fullMass'][lenDct - 1]
                try:
                    vals = [dct[smi] for smi in retDict['targetTautomersWild']]
                    lastVal = vals
                except KeyError:
                    # if no reaction (no value) at this step try to get concentration from prev step i.e. lastVal
                    # if no prev step set, so no reation done and concentration is 0
                    vals = [dct.get(smi, lastVal.get(smi, 0)) for smi in retDict['targetTautomersWild']]
                    lstVal = vals
                if appendMode:
                    # no tautomer summation so add also main form concentration
                    summed.append(sum(vals) + retDict['masses'][target][step])
                else:
                    summed[step] += sum(vals)
    if summed:
        return summed
    return retDict['masses'][target]


def printSummary(data, args):
    print("RESULTS")
    print('input name', 'L2-level yield', 'L3-level yield', 'L4-level yield', sep='\t')
    for name in data:
        conc1, conc2, conc3 = statsMetrics.getFinalProductConc(data[name], args)
        print(name, round(conc1, 1), round(conc2, 1), round(conc3, 1), sep='    \t ')

def doSummaryAndPprint(data, ylds, args):
    if args.pprint:
        print("name", "yield", 'predicted gen0', 'predicted gen1', 'predicted gen2')
    y = []
    x0 = []
    x1 = []
    x2 = []
    results = dict()
    # shortNameToFullName = {fullName.split('_')[0]: fullName for fullName in dataDict}
    for name in sorted(ylds, key=lambda k: ylds[k], reverse=True):
        if args.pprint and name not in data:
            print(name, ylds[name], "NO RES")
        else:
            # fullName = shortNameToFullName[name]
            # print("DATAD", name, fullName in dataDict, dataDict[fullName].keys())
            y.append(ylds[name])
            if args.propagationMode in {'all', 'copasi'}:
                # print("NAME", name, {key: (type(data[name][key]), len(data[name][key])) for key in data[name]})
                conc1, conc2, conc3 = statsMetrics.getFinalProductConc(data[name], args)
                if args.pprint:
                    # print(name, ylds[name], max(data[name][0]), max(data[name][1]), max(data[name][2]))
                    print(name, ylds[name], conc1, conc2, conc3)
                # results[name] = {'exp': ylds[name], 'calced': (max(data[name][0]), max(data[name][1]), max(data[name][2]))}
                results[name] = {'exp': ylds[name], 'calced': (conc1, conc2, conc3)}
                x0.append(conc1)
                x1.append(conc2)
                x2.append(conc3)
            elif args.propagationMode == 'front':
                x0.append(data[name][0])
                x1.append(data[name][1])
                x2.append(data[name][2])
                if args.pprint:
                    print(name, ylds[name], data[name][0], data[name][1], data[name][2])
                results[name] = {'exp': ylds[name], 'calced': (data[name][0], data[name][1], data[name][2])}
            else:
                raise NotImplementedError
    r2 = []
    x0 = numpy.array(x0).reshape((-1, 1))
    x1 = numpy.array(x1).reshape((-1, 1))
    x2 = numpy.array(x2).reshape((-1, 1))
    max01 = [max(x00, x10) for x00, x10 in zip(x0, x1)]
    max012 = [max(x00, x10, x20) for x00, x10, x20 in zip(x0, x1, x2)]
    avg01 = [(x00 + x10) / 2 for x00, x10 in zip(x0, x1)]
    diff01 = [max(0, x10 - x00) for x00, x10 in zip(x0, x1)]
    avgdiff01 = [a + d for a, d in zip(avg01, diff01)]
    # diff
    for x in (x0, x1, x2, max01, max012, avg01, diff01, avgdiff01):
        try:
            # slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
            if args.paropt == 'R2_0':
                reg = LinearRegression(fit_intercept=False).fit(x, y)
            else:
                reg = LinearRegression().fit(x, y)
            r2.append(reg.score(x, y))
            # r2.append(r_value**2)
        except:
            r2.append("FAILED")
    if args.pprint:
        print("R2values:", r2)
    return {'r2': r2, 'results': results}


def getStatsSingle(dataDict, initStoich, adbData, args):
    dataToPrint = dict()
    fh = None
    if args.saveDataForPlot:
        fh = open(args.saveDataForPlot, 'w')

    for name in dataDict:
        molname = name.split('_')[0]
        for key in dataDict[name]:
            sbsStoich = None
            if initStoich and initStoich[molname]:
                sbsStoich = initStoich[molname]
            if args.debug:
                print("====", name, key)
            masses, fullMass, rxChanged, target = calc.propagateNetwork(dataDict[name][key], sbsStoich, adbData, dataDict[name][0]['common'], args)
            if args.debug:
                print("====", name, key)
            if fh:
                print(name, key, formFullCmdMass(fullMass), sep='\t', file=fh)
            if args.debug:
                massDct = {smi: masses[smi] if smi != target else max(masses[smi]) for smi in masses}
                print("MASSES", name, key, target, massDct)
            if args.pprint or args.paropt != 'N':
                if molname not in dataToPrint:
                    dataToPrint[molname] = dict()
                dataToPrint[molname][key] = masses[target]
            else:
                print(name, key, ">>>>>", masses[target], len(masses))
    if fh:
        fh.close()
    return dataToPrint


def formFullCmdMass(listOfDict):
    dictOfList = {smi: [] for smi in listOfDict[-1].keys()}
    for dct in listOfDict:
        for smi in dictOfList:
            dictOfList[smi].append(dct.get(smi, 0))
    return dictOfList


def formFullRxMass(dictOfRxes, dataDict):
    # dictOfList = {smi: [] for smi in DictOfRxes.keys()}
    def _tupleToStr(rxtuple):
        # RXT (('CC(CC(C)([O-])CC(C)(C)O)=[N+]1CCCCC1',), 'CC(CC(C)(O)CC(C)(C)O)=[N+]1CCCCC1', (1,))
        sbses, prod, rxidList = rxtuple
        rxsmi = '.'.join(sbses) + '>>' + prod
        rxstr = f'{rxsmi}::{".".join([str(rxid) for rxid in rxidList])}'
        return rxstr

    rxNumToFull = {dataDict['rxtupleToId'][rxtuple]: _tupleToStr(rxtuple) for rxtuple in dataDict['rxtupleToId']}
    revTmp = {f'{rxnum}_REV': rxNumToFull[rxnum]+'_REV' for rxnum in rxNumToFull}
    rxNumToFull.update(revTmp)
    retDict = {rxNumToFull[rxnum]: dictOfRxes[rxnum] for rxnum in dictOfRxes}
    emptyKey = []
    for key in retDict:
        if not retDict[key]:
            emptyKey.append(key)
    for key in emptyKey:
        del retDict[key]
    return retDict


def setParams(args, combination, orderedParams):
    for name, value in zip(orderedParams, combination):
        args.__dict__[name] = value
    return args


def rangeToList(minim, maxim, step):
    values = []
    val = minim
    while val < maxim:
        values.append(val)
        val += step
    return values


def doParOpt(dataDict, args, params, expYields, initStoich):
    orderedParams = tuple(sorted(params.keys()))
    bestsVal = [0, 0, 0]
    bestsSet = [None, None, None]
    print("PARAMS ORDER", orderedParams)
    qin, qout = Queue(), Queue()

    def _calcOneCombination(qin, qout, args):
        while True:
            combination = qin.get()
            if combination == 'STOP':
                break
            try:
                args = setParams(args, combination, orderedParams)
                res = getStats(dataDict, initStoich, args)
                if args.paropt == 'R2':
                    r2, full = statsMetrics.getR2(res, expYields, args)
                elif args.paropt == 'R2_0':
                    r2, full = statsMetrics.getR2_0(res, expYields, args)
                elif args.paropt == 'MAE':
                    r2, full = statsMetrics.getMAE(res, expYields, args)
                elif args.paropt == 'RMSE':
                    r2, full = statsMetrics.getRMSE(res, expYields, args)
                elif args.paropt == 'pearson':
                    r2, full = statsMetrics.getPearson(res, expYields, args)
                answer = (combination, r2)
            except:
                answer = (combination, 'FAILED')
            qout.put(answer)
        return 1

    procs = [Process(target=_calcOneCombination, args=(qin, qout, args)) for _ in range(args.numprocs)]
    _ = [p.start() for p in procs]

    toGet = 0
    print("Processes stated", len(procs))
    for combination in itertools.product(*[params[k] for k in orderedParams]):
        qin.put(combination)
        toGet += 1
    print("task submitted", toGet)
    for idx in range(toGet):
        msg = qout.get()
        if msg[1] == 'FAILED':
            continue
        if args.debug:
            print("COMBINATION and res", msg)
        newBest = False
        combination, r2 = msg
        for poz in range(3):
            if r2[poz] > bestsVal[poz]:
                bestsVal[poz] = r2[poz]
                bestsSet[poz] = combination
                newBest = True
        if newBest:
            print("BEST", bestsVal, bestsSet)


def doParOptSingle(dataDict, args, params, expYields, initStoich):
    orderedParams = tuple(sorted(params.keys()))
    bestsVal = [0, 0, 0]
    bestsSet = [None, None, None]
    print("PARAMS ORDER", orderedParams)
    for combination in itertools.product(*[params[k] for k in orderedParams]):
        args = setParams(args, combination, orderedParams)
        try:
            res = getStats(dataDict, initStoich, args)
            if args.paropt == 'R2':
                r2, full = statsMetrics.getR2(res, expYields, args)
            elif args.paropt == 'R2_0':
                r2, full = statsMetrics.getR2_0(res, expYields, args)
            elif args.paropt == 'MAE':
                r2, full = statsMetrics.getMAE(res, expYields, args)
            elif args.paropt == 'RMSE':
                r2, full = statsMetrics.getRMSE(res, expYields, args)
            elif args.paropt == 'pearson':
                r2, full = statsMetrics.getPearson(res, expYields, args)
            else:
                raise NotImplementedError
            if args.debug:
                print("COMBINATION", combination, "RESULT", r2)
            newBest = False
            for poz in range(3):
                if r2[poz] > bestsVal[poz]:
                    bestsVal[poz] = r2[poz]
                    bestsSet[poz] = combination
                    newBest = True
            if newBest:
                print("BEST", bestsVal, bestsSet)
        except:
            if args.debug:
                print("PROBLEM WITH PARAMS", combination)
