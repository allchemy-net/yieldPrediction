import json, hashlib
# from functools import partial
from argparse import Namespace
import argsParser, calc, consts, helpers, optscan, loader, stats, calcInitC, loaderHelper, loaderChem, export

# this should helps on multicores machine see https://github.com/PKU-DAIR/open-box/issues/73#issuecomment-1870877934
import os
NUM_THREADS = "1"
os.environ["OMP_NUM_THREADS"] = NUM_THREADS         # export OMP_NUM_THREADS=1
os.environ["OPENBLAS_NUM_THREADS"] = NUM_THREADS    # export OPENBLAS_NUM_THREADS=1
os.environ["MKL_NUM_THREADS"] = NUM_THREADS         # export MKL_NUM_THREADS=1
os.environ["VECLIB_MAXIMUM_THREADS"] = NUM_THREADS  # export VECLIB_MAXIMUM_THREADS=1
os.environ["NUMEXPR_NUM_THREADS"] = NUM_THREADS     # export NUMEXPR_NUM_THREADS=1


def applyMasterFile(args):
    fh = open(args.masterconfig)
    data = fh.read()
    fh.close()
    if "'" in data:
        data = data.replace("'", '"').replace('\n', ' ')
    print("config options from asterconfig file:", data)
    masterargs = json.loads(data)
    nargs = helpers.replaceArgFromDict(args, masterargs)
    return nargs


def applyNamespace(fname, orygArgs):
    # fnList, saveDataPlot, debug, numprocs, exportGraph):
    # .input, args.saveDataForPlot, args.debug, args.numprocs, args.exportGraph)
    fh = open(fname)
    data = fh.read()
    fh.close()
    args = eval(data)
    args.paropt = 'N'
    if orygArgs.loadpickles:
        args.loadpickles = orygArgs.loadpickles
        args.input = None
    else:
        args.input = orygArgs.input
    if orygArgs.loadrxdb:
        args.loadrxdb = orygArgs.loadrxdb
        args.adb = None
    args.saveDataForPlot = orygArgs.saveDataForPlot
    args.numprocs = orygArgs.numprocs
    args.exportGraph = orygArgs.exportGraph
    if orygArgs.debug:
        args.debug = orygArgs.debug
    if orygArgs.overwriteStoich:
        args.overwriteStoich = orygArgs.overwriteStoich
    return args

# # archive output


def doArchiveResult(dataToPrint, args):
    fname = getFilenameForArchive(args)
    data = {'args': args.__dict__}
    data.update(dataToPrint)
    fh = open(fname, 'w')
    json.dump(data, fh, indent=2)
    fh.close()


def getFilenameForArchive(args):
    msg = hashlib.sha512()
    for argname in sorted(args.__dict__):
        msg.update(argname.encode())
        msg.update(str(args.__dict__[argname]).encode())
    return msg.hexdigest()


# #################################
import pickle

if __name__ == "__main__":
    args, optSpace = argsParser.parseArgs()
    if args.masterconfig:
        args = applyMasterFile(args)

    if args.namespace:
        # read all parameters from namespace except input filenames an savaDataForPlot
        args = applyNamespace(args.namespace, args)
    if args.loadrxdb:
        adbData = export.loadrxdb(args.loadrxdb)
    else:
        adbData = loader.parseAdb(args.adb, args)
        adbData = export.exportRXDBtoJSON(adbData)
    print("calculation running with following configs/options: ", args)
    if args.loadpickles:
        dataDict = dict()
        initStoich = dict()
        for fn in args.loadpickles:
            basename = fn.split('_')[0].split('/')[-1]
            dataDict[basename], initStoich[basename] = export.loadallPickle(fn)
    else:
        heatFormation = None
        adbData = None
        if args.adb:
            adbData = loader.parseAdb(args.adb, args)
        dataDict = loader.loadPickle(args.input, adbData, args)
        allSmiles = calc._getAllSmiles(dataDict)
        if args.useMayr != 'N':
            if args.useMayr == 'exp':
                if args.mayrMatchMode == 'exact':
                    raise NotImplementedError
                mayrData = mayr.parseMayrExp(*args.MayrFilesExp)  # 'Nucleophiles.tsv', 'Electrophiles.tsv')
            elif args.useMayr == 'luosz':
                if args.mayrMatchMode != 'exact':
                    raise NotImplementedError
                mayrData = mayr.parseMayrPred(args.MayrFilesPred)
                assert mayrData['dict']
            else:
                raise NotImplementedError
            # calcRxRate(allSmiles, mayrData)
        else:
            mayrData = {'nucleo': list(), 'electro': list(), 'dict': dict()}
        mutexRxDict = loader.loadMutexRxFile(args.mutexRxFile)
        dataDict = loaderChem.addReactionInfos(dataDict, mayrData, adbData, mutexRxDict, args)
        if args.heatFormation:
            heatFormation = loader.getHeatOfFormation(args)
            dataDict = calc.addReactionDE(dataDict, heatFormation)
        overwritenInitC = None
        if args.overwriteStoich:
            overwritenInitC = loader.loadOverwritenInitC(dataDict, args.overwriteStoich, args)
        initStoich = calcInitC.getInitialStoichiometry(dataDict, adbData, overwritenInitC, args)
        if args.fixTautomerRatio != 'no':
            # add info about tautomers
            dataDict = loader.addTautomerInfoToDataset(dataDict, args)
        if args.fixAcidBaseRatio != 'N':
            if args.addReverse == 'wildToMainOnly':
                raise NotImplementedError
        dataDict = loader.addWildFormOfTargetInfoToDataset(dataDict, args)
        dataDict = loaderHelper.checkIfWaterNeedRm(dataDict, args)
        dataDict = loaderHelper.markFromTargetRx(dataDict, args)
    # do calculation
    #pickle.dump(dataDict, open('example.datadict.pickle', 'wb'))
    if args.input:
        #export.exportToJSON(dataDict)
        export.exportToPickle(dataDict, initStoich)

    elif args.paropt != 'N':
        print("ARGS", args, "OPTspace:", optSpace)
        optSpace = optscan.adjustOptSpace(optSpace, args)
        # doParOpt(dataDict, args, params, expYields, initStoich)
        # optFn = partial(optscan.optFunctionOpenBox, dataDict, initStoich, adbData, args)
        # optscan.openBoxOpt(optFn, optSpace, args)
        # dataDict, initStoich, adbData, optSpace, args
        optscan.openBoxOpt(dataDict, initStoich, adbData, optSpace, args)
    elif args.scan:
        optSpace = optscan.loadOptSpaceFromFile(args.scan)
        optscan.performScan(dataDict, initStoich, adbData, args, optSpace)
    else:
        dataToPrint = stats.getStats(dataDict, initStoich, adbData, args)
        # resultSummary = stats.doSummaryAndPprint(dataToPrint, consts.expYields, args)
        stats.printSummary(dataToPrint, args)
        if args.archiveResult:
            doArchiveResult(resultSummary, args)
