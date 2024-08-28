from rdkit import Chem
import consts


# @jit(nopython=True)
def getStoichiometry(args):
    retDict = dict()
    fh = open(args.sbsStoichFile)
    for line in fh:
        line = line.split('\t')
        molid, sbses, _, _, stoich, hasExces, isCat = line
        isCat = isCat.strip()
        try:
            molid = f'mol{int(molid)}'
        except ValueError:
            molid = f'{molid}'
        sbses = [Chem.CanonSmiles(s) for s in sbses.split('.')]
        stoich = [float(sti) for sti in stoich.split(':')]
        stoichDict = getNormalizedStoichiometryDict(sbses, stoich, hasExces, isCat, args)
        retDict[molid] = stoichDict
        if args.debug:
            print("INFO", molid, stoichDict)
    fh.close()
    return retDict


# @jit(nopython=True)
def getNormalizedStoichiometryDict(sbses, stoich, hasExces, isCat, args):
    minSbs = min(stoich)
    if isCat:
        minSbs = tuple(sorted(stoich))[1]
    scale = 100 / minSbs
    newstoi = {sbs: stoich[poz] * scale for poz, sbs in enumerate(sbses)}
    if isCat:
        # add catalysts
        argmin = min(sorted([(s, poz) for poz, s in enumerate(stoich)]))[1]
        sbsname = sbses[argmin]
        newstoi[sbsname] *= args.multiplyCat
    return newstoi


def getInitialStoichiometry(dataDict, adbData, overwritenC, args):
    if args.stoichiometry == 'fromFile':
        if args.sbsStoichFile:
            initStoich = getStoichiometry(args)
            return initStoich
        else:
            raise consts.WrongOptions
    initStoich = dict()
    for name in dataDict:
        data = dataDict[name][0]
        trueName = name.split('_')[0]
        init = [x for x in data['graph'].nodes if data['graph'].in_degree(x) == 0]
        # top = [x for x in data['graph'].nodes if data['graph'].out_degree(x) == 0]
        zeroGen = [smi for smi in data['data']['allSmiles'] if data['data']['allSmiles'][smi]['generation'] == 0]
        initSet = set(init).union(set(zeroGen))
        init = tuple(initSet)
        for smi in init:
            print("SMI", smi, _getOtherSubstrates(data['graph'], smi, initSet))
        # print("ZERO GEN", zeroGen, "INIT", init, "DD", dataDict[name][0]['data']['idToFullRx'][1].keys())
        if args.stoichiometry == 'equal':
            initConc = {smi: 100 for smi in init}
        elif args.stoichiometry == 'optimal':
            initConc = _getOptimalSbsConcentration(init, data, args)
            if len(set(initConc.values())) != 1:
                print(f"WARNING INIT CONC {name} {initConc}")
        else:
            raise NotImplementedError
        initStoich[trueName] = initConc
    if args.debug:
        print("INIT STOICHIOMETRY", initStoich)
    if args.overwriteStoich:
        initStoich = _doOverwriteInitC(initStoich, overwritenC)
        if args.debug:
            print("INIT STOICHIOMETRY (overwriten)", initStoich)
    return initStoich


def _getOtherSubstrates(graph, smi, allInitSmiles):
    retDict = dict()
    for rx in graph.succ[smi]:
        # limit only to main path reactions
        othersbs =[s for s in graph.pred[rx] if s != smi and s not in allInitSmiles]
        retDict[rx] = othersbs
    return retDict


def _doOverwriteInitC(initStoich, overwriten):
    print("INIT", initStoich)
    print("OV", overwriten)
    for graphName in overwriten:
        if graphName not in initStoich:
            continue
        for cmd in overwriten[graphName]:
            if cmd not in initStoich[graphName]:
                raise consts.WrongData
            initStoich[graphName][cmd] = overwriten[graphName][cmd]
    return initStoich


def _getOptimalSbsConcentration(initCmds, data, args):
    sbsesCount = dict()
    for rx in data['common']['mainPathInfo']['rxes']:
        sbses = rx[0].split('>>')[0].split('.')
        for sbs in sbses:
            if sbs not in sbsesCount:
                sbsesCount[sbs] = 0
            sbsesCount[sbs] += 1
    if args.debug:
        print("SBS usage", sbsesCount)
    initSbs = {smi: 100 * sbsesCount.get(smi, 1) for smi in initCmds}
    return initSbs
