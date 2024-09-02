import argparse, json, pickle, time, copy, re, itertools
import networkx
from rdkit import Chem
from rdkit.Chem import AllChem
import sys
sys.path.append('..')
import loaderCat


def parseArgs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--adb', required=True, type=str, help='file with reaction database in json format')
    parser.add_argument('-i', '--input', required=True, type=str, help='input json file to convert')
    parser.add_argument('-o', '--output', required=True, type=str, help='name of output pickle file')
    # optional
    parser.add_argument('--basename', default='path', type=str)
    parser.add_argument('--includeConditions', default='minimal', choices=['minimal', 'minimalMild'])
    parser.add_argument('--catAddTime', default='peak', choices=['peak', 'jitLate'])
    parser.add_argument('--verbose', default=0, type=int, help='verbosity level')
    parser.add_argument('--scaleSideInOldCondition', default=0.9990131459471931, type=float)
    parser.add_argument('--rmRxMode', default='N', choices=['N', 'mainAwild', 'mainAwildALA'])
    parser.add_argument('--removeLateSides', type=str, default='rm1', choices=['N', 'rm1', 'rm2'])
    parser.add_argument('--redox', type=str, default='rmRedOxWater', choices=['N', 'NAVS', 'NAsidemech', 'rm', 'rmRedOxWater'])
    parser.add_argument('--mayrTryAnionic', type=str, default='')
    parser.add_argument('--mayrFiles', nargs='+', default=['Nucleophiles.tsv', 'Electrophiles.tsv'])
    parser.add_argument('--stoichScale', type=str, default='c100', choices=['N', 'simple2', 'simple', 'uniq', 'c100'])
    parser.add_argument('--stoichiometry', type=str, default='optimal', choices=['equal', 'optimal'])
    parser.add_argument('--fixTautomerRatio', type=str, default='constEnolate', choices=['no', 'constEnolate', 'const', 'DE'])
    parser.add_argument('--tautomerProductIs', default='lessStable', choices=['lessStable', 'rightHand', 'moreStable'])
    parser.add_argument('--fixAcidBaseRatio', type=str, default='condIgnoreWeak', choices=['N', 'Y', 'cond', 'condIgnoreWeak'])
    parser.add_argument('--addReverse', type=str, default='YnoTAB', choices=['N', 'Y', 'wildToMainOnly', 'DE', 'YnoTAB'])
    parser.add_argument('--scaleWaterRx', type=str, default='whenWaterSbs', choices=['no', 'whenWaterSbs', 'whenSbsNoRev'])
    parser.add_argument('--fromTargetDef', type=str, default='sideWithTautoDeProt', choices=['strict', 'sideWithTautoDeProt'])
    args = parser.parse_args()
    return args


enolateC = Chem.MolFromSmarts('[CX3-][CX3](=[OX1])[#6,O]')
enolateO = Chem.MolFromSmarts('[CX3]=[CX3]([O-])[#6,O]')
enoleToKetone = AllChem.ReactionFromSmarts('[CX3:1]=[CX3:2]([OH:3])[#6,O:4]>>[C:1][C:2](=[OX1:3])[*:4]')
tautomerisationStrList = [('[CX3]=[CX3]([OH])[#6,O]', '[CX3:1]=[CX3:2]([OH:3])[#6,O:4]>>[C:1][C:2](=[OX1:3])[*:4]'), ]
tautomerisationDict = {Chem.MolFromSmarts(lst[0]): AllChem.ReactionFromSmarts(lst[1]) for lst in tautomerisationStrList}

condABtoInt = {'SA': -3, 'A': -2 , 'WA': -1, 'N': 0, 'WB': 1, 'B': 2, 'SB': 3}
condIntToAB = {condABtoInt[ab]: ab for ab in condABtoInt}
condTempToInt = {'VL': -2, 'L': -1, 'L+RT': -0.5, 'RT': 0, 'H': 1, 'VH': 2}
condIntToTemp = {condTempToInt[temp]: temp for temp in condTempToInt}


def main(args):
    dicts = buildDicts(args)
    adbData = json.load(open(args.adb))
    dataDict = loadPickle(dicts, adbData, args)
    dataDict, initC = formPickle(dataDict, dicts, adbData, args)
    dataToPickle = {'data': dataDict, 'initC': initC}
    fh = open(args.output, 'wb')
    pickle.dump(dataToPickle, fh)
    fh.close()


def formPickle(dataDict, dicts, adbData, args):
    # allSmiles = _getAllSmiles(dataDict)
    mutexRxDict = dict()
    dataDict = addReactionInfos(dataDict, adbData, mutexRxDict, args)
    overwritenInitC = dict()
    if dicts[0]['conditions'].get('initC', False):
        for cmd in dicts[0]['conditions']['initC']:
            overwritenInitC[cmd] = dicts[0]['conditions']['initC'][cmd]
    initStoich = getInitialStoichiometry(dataDict, adbData, overwritenInitC, args)
    if args.fixTautomerRatio != 'no':
        # add info about tautomers
        dataDict = addTautomerInfoToDataset(dataDict, adbData, args)
    if args.fixAcidBaseRatio != 'N':
        if args.addReverse == 'wildToMainOnly':
            raise NotImplementedError
    dataDict = addWildFormOfTargetInfoToDataset(dataDict, args)
    dataDict = checkIfWaterNeedRm(dataDict, args)
    dataDict = markFromTargetRx(dataDict, args)
    return dataDict, initStoich


def addTautomerInfoToDataset(dataDict, rxdbinfo, args):
    for name in dataDict:
        for gen in dataDict[name]:
            dataDict[name][gen]['data'] = _addTautomerInfo(dataDict[name][gen]['data'], dataDict[name][0]['common'], rxdbinfo, args)
    return dataDict


def _addTautomerInfo(data, commonData, rxdbinfo, args):
    tinfo = {'prod': dict(), 'sbs': dict(), 'full': set(), 'rxnum': dict()}
    isEmpty = True
    if args.fixTautomerRatio == 'DE':
        raise NotImplementedError
    # xProd = args.setTautomerRatio / (1 + args.setTautomerRatio)
    for idrx in data['idToFullRx']:
        rxid = data['idToFullRx'][idrx]['rxid']
        if isinstance(rxid, (str, int)):
            rxid = (rxid, )
        if len(rxid) != 1:
            raise NotImplementedError
        rxid = rxid[0]
        # print("FFF", data['idToFullRx'][idrx], rxdbinfo[rxid][')
        # if '.' in data['idToFullRx'][idrx]['smiles'].split('>>')[0]:#
        if rxid == 1 or not rxdbinfo['rxes'][rxid]['tautomerisation']:
            continue
        sbs, prod = data['idToFullRx'][idrx]['smiles'].split('>>')
        sbs = sbs.split('.')
        if len(sbs) != 1:
            raise NotImplementedError
        sbs = sbs[0]
        prod = prod.split('.')
        if len(prod) != 1:
            raise NotImplementedError
        prod = prod[0]
        # if args.tautomerProductIs == 'lessStable' and data['idToFullRx'][idrx]['DE'] < 0:
        #    xProd = 1 - xProd
        # tid = (sbs, prod, rxid, idrx, xProd)
        lessStable = 1
        if args.tautomerProductIs == 'lessStable' and data['idToFullRx'][idrx]['DE'] < 0:
            lessStable = 0
        elif args.tautomerProductIs == 'moreStable' and data['idToFullRx'][idrx]['DE'] > 0:
            lessStable = 0
        sbsMol = Chem.MolFromSmiles(sbs)
        hasSbsEnolC = sbsMol.HasSubstructMatch(enolateC)
        hasSbsEnolO = sbsMol.HasSubstructMatch(enolateO)
        isTautoEnolate = False
        if hasSbsEnolC and hasSbsEnolO:
            raise NotImplementedError
        elif hasSbsEnolC or hasSbsEnolO:
            prodMol = Chem.MolFromSmiles(prod)
            hasProdEnolC = prodMol.HasSubstructMatch(enolateC)
            hasProdEnolO = prodMol.HasSubstructMatch(enolateO)
            if hasProdEnolC and hasProdEnolO:
                raise NotImplementedError
            elif hasProdEnolC or hasProdEnolO:
                isTautoEnolate = True
        tid = (sbs, prod, rxid, idrx, lessStable, isTautoEnolate)
        tinfo['full'].add(tid)
        if prod not in tinfo['prod']:
            tinfo['prod'][prod] = set()
        tinfo['prod'][prod].add(tid)
        if sbs not in tinfo['sbs']:
            tinfo['sbs'][sbs] = set()
        tinfo['sbs'][sbs].add(tid)
        isEmpty = False
        tinfo['rxnum'][idrx] = tid
    if not isEmpty:
        data['tautomers'] = tinfo
        targetTautomers = _getTargetTautomers(data['targetSmiles'], tinfo)
        if targetTautomers:
            data['targetTautomers'] = targetTautomers
            tautomerWild = _getTautomersWildForms([t for t in targetTautomers if t != data['targetSmiles']], commonData)
            if args.debug:
                print("TARGET", data['targetSmiles'], "TAUTO", targetTautomers, "TARGT WILD FORM", tautomerWild)
            if tautomerWild:
                data['targetTautomersWild'] = tautomerWild
            # data['targetTautomersWild'] =
            # find all prot form
            # print("\n\nTIYYY", data['targetSmiles'], tinfo, "\n\nTTTTT", targetTautomers)
    return data


def _getTautomersWildForms(targetTautomers, commonData):
    tautoWild = []
    for tauto in targetTautomers:
        if tauto in commonData['acidBase']['mainToWild']:
            tautoWild.extend(commonData['acidBase']['mainToWild'][tauto])
    return tautoWild


def checkIfWaterNeedRm(data, args):
    # this is for reevaluate if water-needed sie should be really punished - e.g. reverse side which move sbs to main path may need water
    # if so there is no need to scale down such reaction
    for mol in data:
        if data[mol][0]['common']['mainPathInfo']['isWaterSbs'] or data[mol][0]['common']['mainPathInfo']['isWaterCondNeeded']:
            # water is needed so no action needed
            continue
        for gen in data[mol]:
            if gen == 0:
                continue
            if args.scaleWaterRx == 'whenSbsNoRev':
                for rxnum in data[mol][gen]['data']['idToFullRx']:
                    if isinstance(rxnum, str) and 'REV' not in rxnum:
                        continue
                    rxinfo = data[mol][gen]['data']['idToFullRx'][rxnum]
                    sbses, prod = rxinfo['fullRxSmiles'].split('>>')
                    sbses = sbses.split('.')
                    if 'O' in sbses:
                        data[mol][0]['common']['mainPathInfo']['isWaterSbs'] = True
                        data[mol][0]['common']['mainPathInfo']['isWaterCondNeeded'] = True
                        print("HAS WATER", mol, gen, rxnum)
                        break
    return data


# wild of target
def addWildFormOfTargetInfoToDataset(dataDict, args):
    for name in dataDict:
        for gen in dataDict[name]:
            dataDict[name][gen]['data'] = _addWildFormOfTargetInfo(dataDict[name][gen]['data'], dataDict[name][0]['common'], args)
    return dataDict


def _addWildFormOfTargetInfo(data, common, args):
    targetSmi = data['targetSmiles']
    mainToWild = common['acidBase']['mainToWild']
    if targetSmi not in mainToWild:
        return data
    data['targetWildForms'] = mainToWild[targetSmi].copy()
    return data


def _getTargetTautomers(targetSmiles, tautoInfo, maxIter=10):
    allTauto = {targetSmiles, }
    addedInLastIter = {targetSmiles, }
    for _ in range(maxIter):
        thisIterAdded = set()
        for smi in tautoInfo['prod']:
            if smi not in addedInLastIter:
                continue
            for tautorx in tautoInfo['prod'][smi]:
                sbs = tautorx[0]
                if len(sbs) < 3:
                    raise
                if sbs in allTauto:
                    continue
                thisIterAdded.add(sbs)
        for smi in tautoInfo['sbs']:
            if smi not in addedInLastIter:
                continue
            for tautorx in tautoInfo['sbs'][smi]:
                prd = tautorx[1]
                if len(prd) < 3:
                    print("T")
                    raise
                if prd in allTauto:
                    continue
                thisIterAdded.add(prd)
        if not thisIterAdded:
            if len(allTauto) == 1:
                return None
            return allTauto
        addedInLastIter = thisIterAdded
        allTauto.update(thisIterAdded)
    raise NotImplementedError


def getInitialStoichiometry(dataDict, adbData, overwritenC, args):
    initStoich = dict()
    for name in dataDict:
        data = dataDict[name][0]
        trueName = name.split('_')[0]
        init = [x for x in data['graph'].nodes if data['graph'].in_degree(x) == 0]
        # top = [x for x in data['graph'].nodes if data['graph'].out_degree(x) == 0]
        zeroGen = [smi for smi in data['data']['allSmiles'] if data['data']['allSmiles'][smi] == 0]
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
    if args.verbose > 0:
        print("INIT STOICHIOMETRY", initStoich)
    for cmd in overwritenC:
        initStoich[cmd] = overwritenC[cmd]
        if args.verbose > 0:
            print("INIT STOICHIOMETRY (overwriten)", initStoich)
    return initStoich


def _getOptimalSbsConcentration(initCmds, data, args):
    sbsesCount = dict()
    for rx in data['common']['mainPathInfo']['rxes']:
        sbses = rx[0].split('>>')[0].split('.')
        for sbs in sbses:
            if sbs not in sbsesCount:
                sbsesCount[sbs] = 0
            sbsesCount[sbs] += 1
    if args.verbose > 0:
        print("SBS usage", sbsesCount)
    initSbs = {smi: 100 * sbsesCount.get(smi, 1) for smi in initCmds}
    return initSbs


def _getOtherSubstrates(graph, smi, allInitSmiles):
    retDict = dict()
    for rx in graph.succ[smi]:
        # limit only to main path reactions
        othersbs = [s for s in graph.pred[rx] if s != smi and s not in allInitSmiles]
        retDict[rx] = othersbs
    return retDict


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


def markFromTargetRx(dataDict, args):
    if args.fromTargetDef == 'sideWithTautoDeProt':
        return _markFromTargetRxSideWithTautoDeProt(dataDict)
    elif args.fromTargetDef == 'strict':
        return dataDict
    else:
        raise NotImplementedError


def _markFromTargetRxSideWithTautoDeProt(dataDict):
    for mol in dataDict:
        mainPathSet = set(dataDict[mol][0]['common']['mainPathInfo']['rxes'])
        for gen in dataDict[mol]:
            mainInfo = dataDict[mol][gen]['data']
            allTargets = mainInfo.get('targetTautomers', set()).union(mainInfo.get('targetWildForms', set()))
            allTargets.add(mainInfo['targetSmiles'])
            for smil in mainInfo.get('targetTautomers', []):
                for smi in dataDict[mol][0]['common']['acidBase']['mainToWild'].get(smil, []):
                    allTargets.add(smi)
            for rx in mainInfo['idToFullRx']:
                thisRx = mainInfo['idToFullRx'][rx]
                thisRx['isSideFromTarget'] = False
                rxids = thisRx['rxid']
                if isinstance(rxids, (str, int)):
                    rxids = [rxids, ]
                rxhash = (thisRx['fullRxSmiles'], tuple(rxids))
                if rxhash in mainPathSet:
                    continue
                sbs, _ = thisRx['smiles'].split('>>')
                sbs = set(sbs.split('.'))
                if allTargets.intersection(sbs):
                    thisRx['isSideFromTarget'] = True
    return dataDict


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


def addReactionInfos(dataDict, adbData, mutexDict, args):
    cache = dict()
    # hasProblems = False
    # values = []
    for name in dataDict:
        # commonData = dataDict[name][0]['common']
        mainPathRxid = set([smirxid[1][0] for smirxid in dataDict[name][0]['common']['mainPathInfo']['rxes']])
        mutexToRxOnMainPath = set()
        for mutex1fs, mutex2fs in mutexDict.items():
            isec1 = mutex1fs.intersection(mainPathRxid)
            isec2 = mutex2fs.intersection(mainPathRxid)
            if isec1 and isec2:
                raise NotImplementedError
            if isec1:
                mutexToRxOnMainPath.update(mutex2fs)
            if isec2:
                mutexToRxOnMainPath.update(mutex1fs)
        if mutexToRxOnMainPath.intersection(mainPathRxid):
            print("ISEC", isec1, isec2)
            print("MUT", mutexToRxOnMainPath, "MAIN", mainPathRxid)
            raise NotImplementedError
        if args.verbose:
            print("MUTEX:3409:", name, mutexToRxOnMainPath, "MP", mainPathRxid)
        for gen in dataDict[name]:
            for rxid in dataDict[name][gen]['data']['idToFullRx']:
                print("RX", rxid, dataDict[name][gen]['data']['idToFullRx'][rxid])
                rxsmiles = dataDict[name][gen]['data']['idToFullRx'][rxid]['fullRxSmiles']
                if rxsmiles not in cache:
                    try:
                        mappedSmiles = dataDict[name][gen]['data']['idToFullRx'][rxid]['additionalInfo']['mappedRx']
                    except:
                        mappedSmiles = dataDict[name][gen]['data']['idToFullRx'][rxid]['mappedRx']
                    # fix problem with reaction with reverse order: products>>substrates
                    mappedSmiles = _correctMappedSmiles(mappedSmiles, rxsmiles)
                    rxidOryg = dataDict[name][gen]['data']['idToFullRx'][rxid]['rxid']
                    if isinstance(rxidOryg, (str, int)):
                        rxidOryg = (rxidOryg, )

                    if rxidOryg[0] == 1:
                        # condition changes
                        continue
                    rxInfo = adbData['rxes'][rxidOryg[0]]
                    tryCalcMayr = True
                    if is1c(rxInfo, rxsmiles, mappedSmiles):
                        if args.verbose > 100:
                            print(f"Mayr: no rate for 1c: {rxsmiles} {rxidOryg}")
                        tryCalcMayr = False
                    if not rxInfo['mayr']:
                        if args.verbose > 100:
                            print(f"Mayr: no rate for non AD/SUB: {rxsmiles} {rxidOryg}")
                        tryCalcMayr = False
                    try:
                        mayrRate = None
                        if tryCalcMayr and 'mayr' not in dataDict[name][gen]['data']['idToFullRx'][rxid]:
                            print('!!no mayr rate!!')
                        ringBoostFactor = _calcRingBoost(rxsmiles, mappedSmiles, rxInfo, args)
                        stoichScale = _calcStoichScale(rxsmiles, mappedSmiles, args)

                        isMutex = set(rxidOryg).intersection(mutexToRxOnMainPath)
                        cache[rxsmiles] = {'mayr': mayrRate, 'boostRing': ringBoostFactor, 'stoichScale': stoichScale, 'isMutex': isMutex}
                    except:
                        print(f"\n RXID {rxid} {rxsmiles} \n {dataDict[name][gen]['data']['idToFullRx'][rxid]}")
                        raise
                    # values.append(cache[rxsmiles])
                    if args.verbose > 100:
                        print(f"Mayr RXVAL:: {cache[rxsmiles]} {rxsmiles} {rxidOryg}")
                if cache[rxsmiles]['mayr'] is not None:
                    dataDict[name][gen]['data']['idToFullRx'][rxid]['rate'] = cache[rxsmiles]['mayr']
                if cache[rxsmiles]['boostRing'] is not None:
                    dataDict[name][gen]['data']['idToFullRx'][rxid]['boostRing'] = cache[rxsmiles]['boostRing']
                if cache[rxsmiles]['stoichScale'] is not None:
                    dataDict[name][gen]['data']['idToFullRx'][rxid]['stoichScale'] = cache[rxsmiles]['stoichScale']
                if cache[rxsmiles]['isMutex']:
                    if args.verbose > 1:
                        print(f"MUTEX:3455:: {name} {gen} {rxid} {dataDict[name][gen]['data']['idToFullRx'][rxid]['rxid']}")
                    dataDict[name][gen]['data']['idToFullRx'][rxid]['isMutex'] = cache[rxsmiles]['isMutex']
                # print(f"RATEmayr {name} {gen}  {rxsmiles} ::  {cache[rxsmiles]}")
    return dataDict


def is1c(rxInfo, rxsmiles, mappedRx):
    if rxInfo['nComponents'] != 1:
        return False
    return True


def _getAllSmiles(dataDict):
    allSmiles = set()
    for calcName in dataDict:
        for gen in dataDict[calcName]:
            smiles = set(dataDict[calcName][gen]['data']['allSmiles'].keys())
            allSmiles.update(smiles)
    return allSmiles


def _calcRingBoost(rxsmiles, mappedSmiles, oneRxInfo, args):
    if '1' not in rxsmiles:
        return None
    if oneRxInfo['rearrangement'].upper().strip().startswith('Y'):
        return None
    sbs, prod = mappedSmiles.split('>>')
    prods = prod.split('.')
    mainProdRings, mainProdIdxToNum = _getRingMapNum(prods[0])
    sbsRingIdxInfo = [_getRingMapNum(singlesbs) for singlesbs in sbs.split('.')]
    ringInSbses = []
    sbsMapNumSetList = []
    for sbsRingList, idxToMap in sbsRingIdxInfo:
        ringInSbses.extend(sbsRingList)
        sbsMapNumSetList.append(set(idxToMap.values()))
    formedRing = set(mainProdRings) - set(ringInSbses)
    if not formedRing:
        return None

    for ring in formedRing:
        if len(ring) not in {3, 5, 6, 7}:
            continue
        for sbs in sbsMapNumSetList:
            if ring.issubset(sbs):
                if args.debug > 10:
                    print(f"RXboost 1 {rxsmiles} formed:: {formedRing} mapp: {sbsMapNumSetList} MAPP: {mappedSmiles} ringIdx: {sbsRingIdxInfo}")
                # then it is multiply by args.boostIntra
                return 1
    if args.debug > 10:
        print(f"RXboost 0 {rxsmiles} formed:: {formedRing} mapp: {sbsMapNumSetList} MAPP: {mappedSmiles} ringIdx: {sbsRingIdxInfo}")
    return None


def _getRingMapNum(smiles):
    mol = Chem.MolFromSmiles(smiles)
    ringInfo = mol.GetRingInfo()
    idxToNum = {atm.GetIdx(): atm.GetAtomMapNum() for atm in mol.GetAtoms()}
    ringAtomMapNums = [frozenset([idxToNum[atmid] for atmid in ring]) for ring in ringInfo.AtomRings()]
    return ringAtomMapNums, idxToNum


def _calcStoichScale(rxsmiles, mappedSmiles, args):
    if args.stoichScale in {'N', 'c100'}:
        return None
    elif args.stoichScale == 'simple':
        numsbs = len(rxsmiles.split('>>')[0].split('.'))
        if numsbs > 1:
            return numsbs
    elif args.stoichScale == 'simple2':
        numsbs = len(rxsmiles.split('>>')[0].split('.'))
        if numsbs > 2:
            print("RXSMII", rxsmiles)
            return numsbs - 1
    elif args.stoichScale == 'uniq':
        sbses = rxsmiles.split('>>')[0].split('.')
        counts = [sbses.count(s) for s in set(sbses)]
        if max(counts) > 1:
            print("UNIQ", rxsmiles)
            return max(counts)
    else:
        raise NotImplementedError
    return None


def _correctMappedSmiles(mappedSmiles, rxsmiles):
    sbs, prod = rxsmiles.split('>>')
    sbsP, prodP = sbs.count('+'), prod.count('+')
    mappedSbs, mappedProd = mappedSmiles.split('>>')
    if prodP != sbsP:
        mappedSbsP, mappedProdP = mappedSbs.count('+'), mappedProd.count('+')
        if sbsP != mappedSbsP or prodP != mappedProdP:
            if sbsP == mappedProdP and prodP == mappedSbsP:
                return mappedProd + '>>' + mappedSbs
            else:
                print("PROBLEM NOT FIXED", rxsmiles)
        return mappedSmiles
    uniqSbs = [getUniqSmiles(Chem.MolFromSmiles(smi)) for smi in mappedSbs.split('.')]
    sbs = sbs.split('.')
    prod = prod.split('.')
    if set(uniqSbs) == set(sbs):
        return mappedSmiles
    # print("REVERSE::", rxsmiles, "MAPPEDBEFOR", mappedSmiles)
    return mappedProd + '>>' + mappedSbs


def loadPickle(dicts, rxInfo, args):
    dataDict = {args.basename: dict()}
    for gen, data in enumerate(dicts):
        dataDict[args.basename][gen] = {'data': data, }

    wildToMain = dict()  # smiles : smiles
    mainToWild = dict()  # smiles : {smiles1, smiles2}
    wildFormType = dict()
    noMainForm = set()
    name = args.basename
    fixed = _detectFixed(dataDict[name][0]['data']['rxToFullDict'])
    additives = _detectAdditives(dataDict[name])

    mainPathInfo = getMainPathInfo(dataDict[name][0], rxInfo, args)

    allowedCond, mainPathInfo = detectAllowedCond(dataDict[name][0], rxInfo, mainPathInfo, args)
    if args.verbose > 100:
        print(f"ALLOWEDCOND::76: {allowedCond}")
    theCond = getWorkingConditions(allowedCond)

    allFixed = set()
    for rx in fixed:
        _ = [allFixed.add(cmd) for cmd in fixed[rx]]
    allInit = mainPathInfo['initSubstrates'].union(allFixed)
    isOnepot = dicts[0]['conditions']['isOnepot']
    dataDict[name][0]['common'] = {'fixed': fixed, 'allowedCond': allowedCond, 'mainPathInfo': mainPathInfo,
                                   'workingCond': theCond, 'allInitSbs': allInit, 'isOnepot': isOnepot}
    if args.catAddTime == 'init':
        pass
        # no acction required
    elif args.catAddTime in {'jitEarly', 'jitLate'}:
        rxEnablePoint = loaderCat.getRxEnablePoint(allowedCond, mainPathInfo, args)
        if rxEnablePoint:
            dataDict[name][0]['common']['rxEnablePoint'] = rxEnablePoint
    elif args.catAddTime == 'peak':
        rxEnableCmd = loaderCat.getRxEnableCmd(mainPathInfo, args)
        if rxEnableCmd:
            dataDict[name][0]['common']['rxEnableCmd'] = rxEnableCmd
    else:
        raise NotImplementedError
    dataDict[name][0]['common']['onePot'] = loaderCat.getRxEnableCmdOnepot(mainPathInfo, args)
    # print("THEEE", theCond, len(allowedCond['full'][0]))
    for gen in dataDict[name]:
        try:
            dataDict[name][gen]['data']['idToFullRx'] = _mapIdToFullRx(dataDict[name][gen]['data'], args)
            wildToMain, noMainForm, mainToWild, wildFormType, allWilds = formWildInfo(dicts[gen]['conditions']['formAB'])
            # print(dicts
            dataDict[name][0]['common']['acidBase'] = {'noMainForm': noMainForm, 'wildToMain': wildToMain, 'mainToWild': mainToWild, 'wildFormType': wildFormType, 'wildForms': allWilds}
            if args.scaleSideInOldCondition:
                dataDict[name][gen]['data']['idToFullRx'] = _addScalingForSideInOldConditions(dataDict[name][gen]['data'],
                                                                                              args, dataDict[name][0]['common'])
            # dataDict[name][gen]['data']['idToFullRx'] = _addShortcutInfo(dataDict[name][gen]['data']['idToFullRx'], shortcutToBanSet, args)
            # dataDict[name][gen]['data']['idToFullRx'] = _addDeprotonationInfo(dataDict[name][gen]['data']['idToFullRx'], rxInfo, args)
            dataDict[name][gen]['data']['rxAvail'] = loaderCat.addRxAvailInfo(dataDict[name][gen]['data'], dataDict[name][0]['common'], args)
            # below also ['data']['idToFullRx'] may be changed
            networkxGraph, dataDict[name][gen]['data'] = _buildGraph(dataDict[name][gen]['data']['graph'], additives,
                                                                     dataDict[name][gen]['data'], rxInfo, dataDict[name][0]['common'], args)
            # dataDict[name][gen]['data'] is updated dataDict[name][gen]['data']['idToFullRx']
            rxPozToRm = set()
            if rxPozToRm:
                networkxGraph, removedRxPoz = removeRxFromGraph(networkxGraph, rxPozToRm, allInit, args)
                if args.debug > 8:
                    tmpdata = dataDict[name][gen]['data']['idToFullRx']
                    inforxtorm = [(tmpdata[rxpoz]['smiles'], tmpdata[rxpoz]['rxid']) for rxpoz in rxPozToRm]
                    if args.debug:
                        print(f"REMOVED_RX_FROM_GR:91 {rxPozToRm} RM: {inforxtorm}")
                dataDict[name][gen]['data'] = removeRxFromData(dataDict[name][gen]['data'], removedRxPoz)
            printRx = False
            if printRx:
                for rxnum in dataDict[name][gen]['data']['idToFullRx']:
                    dct = dataDict[name][gen]['data']['idToFullRx'][rxnum]
                    print("F", dct['rxid'], dct['smiles'])
            dataDict[name][gen]['graph'] = networkxGraph
        except:
            print("PROBLEM WITH", name, gen)
            raise
    return dataDict


def formWildInfo(abinfo):
    # dicts[gen]['formAB'])
    # print("===ab===", abinfo)
    wildToMain, mainToWild, = dict(), dict()
    noMainForm, wildFormType = dict(), dict()
    allWilds = set()
    for wildForm in abinfo:
        formType, mainForm = abinfo[wildForm]
        wildToMain[wildForm] = mainForm
        if mainForm in mainToWild:
            raise ValueError
        mainToWild[mainForm] = wildForm
        wildFormType[wildForm] = wildFormType
        allWilds.add(wildForm)
    return wildToMain, noMainForm, mainToWild, wildFormType, allWilds


def _getLateSideToRm(graph, data, commonData, args):
    # assign generation
    smiToGen, toAdd = _assignSmilesGeneration(graph, commonData)
    # print(data['allSmiles'].keys())
    notfound = []
    for smi in toAdd:
        if smi not in data['allSmiles']:
            notfound.append(smi)
            continue
        # print(smi, data['allSmiles'][smi]['generation'])
        smiToGen[smi] = data['allSmiles'][smi]
    if notfound:
        print('missed  generation for smiles', notfound)
        raise

    rxToRm = []
    for reasbs in commonData['mainPathInfo']['reactiveSubstrates']:
        # print("SN", reasbs, data['allSmiles'][reasbs]['generation'], data['idToFullRx'].keys())
        node = [no for no in graph.nodes() if no == reasbs][0]
        nextRxes = tuple(graph.successors(node))
        if len(nextRxes) == 1:
            # is involved in one reaction so it must by on main path - no late side
            continue
        mainRxes = []
        mainRxProdToGen = dict()
        sideRxProdToGenAndNum = dict()
        debug = False
        for rxnum in nextRxes:
            rxinfo = data['idToFullRx'][rxnum]
            rxhash = tuple([rxinfo['smiles'], tuple(rxinfo['rxid'])])
            prd = rxinfo['smiles'].split('>>')[1]
            if prd not in smiToGen:
                print("SMITOGEN", smiToGen)
                print("REMOVED FROM GRAPH??", prd, rxhash)
                raise
            gen = smiToGen[prd]
            if rxhash in commonData['mainPathInfo']['rxes']:
                mainRxes.append(rxnum)
                mainRxProdToGen[prd] = gen
            else:
                sideRxProdToGenAndNum[prd] = (gen, rxnum)
        if debug:
            print("MAIN", mainRxes, mainRxProdToGen, "SIDE", sideRxProdToGenAndNum)
        mainProdsGen = set(mainRxProdToGen.values())
        if not mainProdsGen:
            if args.debug:
                print("REA__", reasbs, "NEXT", nextRxes, "REAS", commonData['mainPathInfo']['reactiveSubstrates'])
            continue
        # assert len(mainProdsGen) == 1
        allowedGen = list(range(min(mainProdsGen)))
        allowedGen.extend(mainProdsGen)
        if args.removeLateSides == 'rm1':
            # nothing to add
            pass
        elif args.removeLateSides == 'rm2':
            _ = [allowedGen.append(g+1) for g in mainProdsGen]
        else:
            raise NotImplementedError
        allowedGen = set(allowedGen)
        for prd in sideRxProdToGenAndNum:
            gen, rxnum = sideRxProdToGenAndNum[prd]
            if gen in allowedGen:
                continue
            if args.debug:
                print("lateSIde: RM", gen, rxnum, "NOT in", allowedGen, data['idToFullRx'][rxnum]['smiles'])
            rxToRm.append(rxnum)
    return rxToRm


def _getRedoxToRm(graph, data, commonData, rxInfo, args):
    # args.redox in {'rm', 'rmRedOxWater'}
    redoxInfo = commonData['mainPathInfo']['redox']
    mainPathRxes = set()
    rxPozToRm = []
    for smarx in commonData['mainPathInfo']['redox']:
        _ = [mainPathRxes.add(rx) for rx in smarx[1]]
    for rxpoz in data['idToFullRx']:
        rxnums = data['idToFullRx'][rxpoz]['rxid']
        if isinstance(rxnums, (int, str)):
            rxnums = (rxnums, )
        for rxnum in rxnums:
            if rxnum == 1:
                continue
            if rxnum in mainPathRxes:
                continue
            rxredox = rxInfo['rxes'][rxnum]['rxoxred']
            if not rxredox or rxredox == 'NA':
                continue
            if not redoxInfo:
                if args.redox == 'rm':
                    rxPozToRm.append(rxpoz)
                    break
                elif args.redox == 'rmRedOxWater':
                    print("RXRED:137", rxredox)
                    if rxredox == 'RED':
                        rxPozToRm.append(rxpoz)
                        break
                else:
                    raise NotImplementedError
            if rxredox != 'OX':
                print("RX", rxpoz, "RXRED:", rxredox, rxnum)
                print("COMMON", redoxInfo)
                # raise
    return rxPozToRm


def removeRxFromGraph(graph, rxPozToRm, initSubstrates, args):
    # print("ALLEDGESbefore:", graph.edges())
    removedRxes = set()
    problems = []
    removedNoProblem = []
    itr = 0
    if args.verbose > 5:
        print("RX TORM0", rxPozToRm, "INIT", initSubstrates)
    while rxPozToRm:
        removedRxesThis = set()
        for rxPoz in set(rxPozToRm):
            try:
                prods = tuple(graph.successors(rxPoz))
                # sbses = tuple(graph.predecessors(rxPoz))
                # print(f"LONGCUT: {rxPoz} RM RX", '.'.join(sbses) + '>>'+'.'.join(prods), f"Is sbs: {tuple(graph.successors(prods[0]))}")
                graph.remove_node(rxPoz)
                removedRxes.add(rxPoz)
                graph, nextRxPozToRm = _removeProdsAndReturnNextRxToRm(graph, prods, initSubstrates, args)
                if args.debug > 20:
                    print(f"RMpart {rxPoz} next: {nextRxPozToRm}")
                removedRxesThis.update(nextRxPozToRm)
                removedNoProblem.append(rxPoz)
            except networkx.exception.NetworkXError:
                problems.append((rxPoz, itr))
                continue
        rxPozToRm = removedRxesThis
        itr += 1
        if args.verbose > 5:
            print(f"RX TO RM ITR {itr} :: {rxPozToRm}")
    if problems:
        print("non existing nodes:", problems, "REMOVED:", removedNoProblem, "RXP", rxPozToRm)
    # print("ALLEDGESafter:", graph.edges())
    if args.verbose > 0:
        print("REMOVED RX from graph", removedRxes, "PROBLEMS?", problems, "RMED", removedNoProblem)
    return graph, removedRxes


def _removeProdsAndReturnNextRxToRm(graph, prods, initSubstrates, args):
    rxToRm = set()
    if args.debug > 123:
        print(f"toremoveal {prods} INIT: {initSubstrates}")
    for prod in prods:
        if prod in initSubstrates:
            if args.debug > 10:
                print(f"DONTrm:: {prod} is in initial substrate ")
            continue
        rxesToTheProd = tuple(graph.predecessors(prod))
        if not rxesToTheProd:
            # during mass initialization nodes without predecessor are treated as starting (mass is added to them)
            # therefore it is important to remove products (of cource if there is no other reaction(s)) together with reaction
            for nextrx in graph.successors(prod):
                rxToRm.add(nextrx)
            graph.remove_node(prod)
            if args.debug > 5:
                print(f"YES_REMOVE_NODE {prod}")
        else:
            if args.debug > 5:
                print(f"NOT_REMOVE_NODE {prod} {rxesToTheProd}")
    return graph, rxToRm


def _assignSmilesGeneration(graph, commonData):
    smiToGen = dict()
    for cmd in commonData['mainPathInfo']['smilesToGen']:
        smiToGen[cmd] = commonData['mainPathInfo']['smilesToGen'][cmd]
    toAdd = set([smiOrInt for smiOrInt in graph.nodes() - set(smiToGen.keys()) if isinstance(smiOrInt, str)])
    mainByprods = commonData['mainPathInfo']['byprods']['main']
    # alterByprods = commonData['mainPathInfo']['byprods']['alter']
    # prevToAdd = None
    while toAdd:
        toRm = set()
        for cmd in toAdd:
            for rx in graph.pred[cmd]:
                if not all(sbs in smiToGen for sbs in graph.pred[rx]):
                    continue
                gen = 1 + max(smiToGen[sbs] for sbs in graph.pred[rx])
                toRm.add(cmd)
                if cmd not in smiToGen:
                    smiToGen[cmd] = gen
                else:
                    if smiToGen[cmd] > gen:
                        smiToGen[cmd] = gen
        if not toRm:
            for smi in toAdd:
                if smi not in mainByprods:
                    continue
                gen = _getByprodsGen(graph, commonData, mainByprods[smi], smiToGen)
                if gen is None:
                    continue
                smiToGen[smi] = gen
                toRm.add(smi)
        if not toRm:
            break
        toAdd -= toRm
    return smiToGen, toAdd


def _getByprodsGen(graph, commonData, byprodrx, smiToGen):
    susprx = [smarxid for smarxid in commonData['mainPathInfo']['rxes'] if smarxid[1][0] in byprodrx]
    if not susprx:
        return None
    therx = None
    if len(susprx) == 1:
        therx = susprx[0]
    else:
        raise NotImplementedError
    sbs, prod = therx[0].split('>>')
    sbs = sbs.split('.')
    sbsgen = [smiToGen[s] for s in sbs]
    prodgen = smiToGen[prod]
    print("GENS", sbsgen, "PG", prodgen)
    if max(sbsgen) + 1 == prodgen:
        return prodgen

    print("XX", [smarxid[1][0] for smarxid in commonData['mainPathInfo']['rxes']])
    print("BP", byprodrx)
    print("S", susprx)
    print("FIX", commonData['fixed'], "commonK", commonData.keys())
    raise


def _buildGraph(edgesList, additives, data, rxInfo, commonData, args):
    graph = networkx.DiGraph()
    for edge in edgesList:
        graph.add_edge(*edge)
        if args.verbose > 100:
            print("ADDED EDGE", edge)

    rxPozToRm = []
    if args.includeConditions != 'no':
        if args.includeConditions in {'minimal', 'minimalMild'}:
            rxPozToRmDueToCond = _getRxPozToRmDueCond(data, commonData, rxInfo, args)
            rxPozToRm.extend(rxPozToRmDueToCond)
            if args.verbose > 0:
                print(f"RM due to:conditions: {rxPozToRm}")
        else:
            raise NotImplementedError
    # if args.rmRxMode != 'N':
    #    rxToRm = _getRxPozToRmDueToWildLA(data, commonData, rxInfo, args)
    #    if args.verbose > 0:
    #        print("RXTORM", rxToRm)
    #    rxPozToRm.extend(rxToRm)
    if args.removeLateSides != 'N':
        rxToRmDueToReactiveSbs = _getLateSideToRm(graph, data, commonData, args)
        if args.verbose > 0:
            print(f"rm rx due to ractive sbs {rxToRmDueToReactiveSbs}")
        rxPozToRm.extend(rxToRmDueToReactiveSbs)
    if args.redox in {'rm', 'rmRedOxWater'}:
        rxToRmDueToRedox = _getRedoxToRm(graph, data, commonData, rxInfo, args)
        if args.verbose > 0:
            print("RX TO REDOX", rxToRmDueToRedox, "XX", rxPozToRm)
        rxPozToRm.extend(rxToRmDueToRedox)
    if rxPozToRm:
        initSubstrates = commonData['allInitSbs']
        graph, removedRxPoz = removeRxFromGraph(graph, rxPozToRm, initSubstrates, args)
        if args.verbose > 0:
            print("RX TO rm", removedRxPoz, "XX", rxPozToRm)
        if args.verbose > 8:
            print(f"REMOVED_RX_FROM_GR:50 {rxPozToRm} RM: {[(data['idToFullRx'][rxpoz]['smiles'], data['idToFullRx'][rxpoz]['rxid']) for rxpoz in rxPozToRm]}")
            print("RXIDS", [data['idToFullRx'][rxpoz]['rxid'] for rxpoz in rxPozToRm])
        data = removeRxFromData(data, removedRxPoz)
    # check add missing side reagents
    newDataIdToFullRx = dict()
    for rxnum in data['idToFullRx']:
        rxids = data['idToFullRx'][rxnum]['rxid']
        if isinstance(rxids, (str, int)):
            rxids = (rxids, )
        try:
            sbses = set(graph.predecessors(rxnum))
            prods = set(graph.successors(rxnum))
        except networkx.exception.NetworkXError:
            if args.verbose > 0:
                print(f"rxpoz {rxnum} NOT IN GRAPH")
            continue
        rsbses, rprods = data['idToFullRx'][rxnum]['fullRxSmiles'].split('>>')
        isReversible = False
        if args.addReverse == 'Y':
            isReversible = set(data['idToFullRx'][rxnum]['rxid']).intersection(rxInfo['reversibleRx'])
        elif args.addReverse == 'YnoTAB':
            isReversible = set(data['idToFullRx'][rxnum]['rxid']).intersection(rxInfo['reversibleRx'])
            if _isWildFormingRx(data['idToFullRx'][rxnum]):
                if args.verbose > 100:
                    print("is wild forming", rxnum, data['idToFullRx'][rxnum])
                isReversible = False
            if any(rxInfo[rxid]['tautomerisation'] for rxid in rxids):
                if args.verbose > 100:
                    print("is tautomerisation", rxnum, data['idToFullRx'][rxnum]['rxid'], data['idToFullRx'][rxnum]['smiles'])
                isReversible = False
        elif args.addReverse == 'wildToMainOnly':
            isReversible = _isWildFormingRx(data['idToFullRx'][rxnum])
        if isReversible:
            if _hasGraphReverseRx(graph, data['idToFullRx'][rxnum], data['idToFullRx']):
                if args.verbose > 100:
                    print("has reversible", rxids, data['idToFullRx'][rxnum])
                isReversible = False
            if args.filterReverse == 'arom':
                if _isAromFormingRx(data['idToFullRx'][rxnum]):
                    isReversible = False
        rprods = set(rprods.split('.'))
        rsbses = set(rsbses.split('.'))
        if isReversible:
            rxnumRev = f'{rxnum}_REV'
            # if rxidRev == '38_REV':
            #    pri nt("::::::::", rxidRev, "S", rsbses, "P", rprods, "EL", len(edgesList))
            for sbs in rsbses:
                graph.add_edge(rxnumRev, sbs)
            for prod in rprods:
                graph.add_edge(prod, rxnumRev)
            newDataIdToFullRx[rxnumRev] = _formReverseRx(data['idToFullRx'][rxnum])
            print("ADD REVERSIBLE", rxnumRev, newDataIdToFullRx[rxnumRev])
        if sbses == rsbses and prods == rprods:
            if args.debug > 100:
                print(f">== ignore (reaction with exactly the same prod and sbs are in graph):: {rxnum} rsbses: {rsbses} {sbses} prods:: {rprods} {prods}")
                print(">== FULL", data['idToFullRx'][rxnum])
            continue
        for sbs in rsbses - sbses:
            if args.debug > 10:
                print("ADD edge to side from sbs", sbs, rxnum)
            graph.add_edge(sbs, rxnum)
        for prod in rprods - prods:
            if args.debug > 10:
                print("ADD edge to side to prod", rxnum, prod)
            graph.add_edge(rxnum, prod)
        if args.debug:
            print("==>", rxnum, sbses, prods, data['idToFullRx'][rxnum]['fullRxSmiles'])
    print("NEW DATA TO FULLRX KEYS::", newDataIdToFullRx.keys())
    if newDataIdToFullRx:
        data['idToFullRx'].update(newDataIdToFullRx)
    if args.verbose > 100:
        print(f"FULLGRAPH:: {graph.edges}")
    print("DATA IDTOFULLRX", data['idToFullRx'].keys())
    return graph, data


def _isWildFormingRx(rxinfo, ):
    print("IS WILD", rxinfo)
    #pp
    raise


def _formReverseRx(data):
    newDict = {k: data[k] for k in data}
    newDict['smiles'] = '>>'.join(reversed(data['smiles'].split('>>')))
    newDict['fullRxSmiles'] = '>>'.join(reversed(data['fullRxSmiles'].split('>>')))
    newDict['origin'] = 'hand-made reverse rx'
    # print("DATA", data)
    if 'mappedRx' in newDict:
        newDict['mappedRx'] = '>>'.join(reversed(data['mappedRx'].split('>>')))
    else:
        newDict['additionalInfo']['mappedRx'] = '>>'.join(reversed(data['additionalInfo']['mappedRx'].split('>>')))
    # newDict['mappedRx'] = '>>'.join(reversed(data['mappedRx'].split('>>')))
    return newDict


def _hasGraphReverseRx(graph, rxinfo, allrxinfo):
    sbses, prod = rxinfo['smiles'].split('>>')
    sbses = sbses.split('.')
    if '.' in prod:
        raise NotImplementedError
    # paths = dict()
    sbsesSet = set(sbses)
    for sbs in sbses:
        path = tuple(networkx.all_simple_paths(graph, source=prod, target=sbs, cutoff=2))
        if not path:
            continue
        revRxCandidateNums = [p[1] for p in path]
        for revRxCandNum in revRxCandidateNums:
            try:
                revCandSbses, revCandProds = allrxinfo[revRxCandNum]['fullRxSmiles'].split('>>')
            except:
                print(f"REV PROBLEM:: {revRxCandNum} ALLNUM:: {allrxinfo.keys()} ::PATH:: {path} GRAPH:: {graph.edges}")
                # GRAPH:: [(1, 'CCCC(C(C(=O)O)C(=O)O)N1CCCCC1'), (1, 'C[O-]'), (1, '[K+]'), (1, '[Na+]'), ('CCCC(C(C(=O)O)C(=O)O)N1CCCCC1', 7),
                # ('CCCC(C(C(=O)O)C(=O)O)N1CCCCC1', '1_REV'), ('C[O-]', 6), ('C[O-]', '1_REV'), ('CCCC(C(C(=O)OC)C(=O)OC)N1CCCCC1', 1), ('CCCC(C(C(=O)OC)C(=O)OC)N1CCCCC1', 4),
                # ('CCCC(C(C(=O)OC)C(=O)OC)N1CCCCC1', 9), (4, 'C[O-]'), (9, 'CCCC(C(C(=O)O)C(=O)OC)N1CCCCC1'), (7, 'CCCC(C(C(=O)[O-])C(=O)O)N1CCCCC1'), (6, 'CO'), ('[OH-]', 1),
                #  ('[OH-]', 5), (5, 'O'), ('1_REV', 'CCCC(C(C(=O)OC)C(=O)OC)N1CCCCC1'), ('1_REV', '[OH-]')
                # raise
                continue
            revCandProdsSet = set(revCandProds.split('.'))
            if revCandProdsSet == sbsesSet:
                return True
    return False


def _minimizeMinimalCond(minAB, minTemp, topoOrder, rxGraph):
    bestAB = _getMinimalMild(minAB, mode='AB')
    bestT = _getMinimalMild(minTemp, mode='T')
    mintopo = list((frozenset((ab, )), frozenset((temp, ))) for ab, temp in zip(bestAB, bestT))
    return mintopo


def _getMinimalMild(minAB, mode):
    allAB = [len(elem) for elem in minAB]
    if allAB[0] == 1 and len(set(allAB)) == 1:
        # no conditions range so return whan we have at this moment
        return [tuple(elem)[0] for elem in minAB]
    uniqAB = [minAB[0], ]
    for cond in minAB[1:]:
        if cond == uniqAB[-1]:
            continue
        uniqAB.append(cond)
    if mode == 'AB':
        uniqABint = [[condABtoInt[cnd] for cnd in cond] for cond in uniqAB]
        bestABint = _getBestCondCombinations(uniqABint)
        # empty frozenset if fo reaction in which only LA is spec. in conditionsAB
        bestAB = [condIntToAB.get(cint, frozenset()) if cint is not None else set() for cint in bestABint]
    elif mode == 'T':
        uniqABint = [[condTempToInt[cnd] for cnd in cond] for cond in uniqAB]
        bestABint = _getBestCondCombinations(uniqABint)
        bestAB = [condIntToTemp[cint] for cint in bestABint]
    else:
        raise NotImplementedError
    # print("UNIQ", uniqABint, "BEST", bestABint, "B", bestAB)
    bestCond = []
    uniqPoz = 0
    prevMin = None
    for cond in minAB:
        if prevMin and prevMin != cond:
            uniqPoz += 1
        bestCond.append(bestAB[uniqPoz])
        prevMin = cond
    return bestCond


def _getBestCondCombinations(uniqCondInt):
    bestAB = None
    bestABval = None
    uniqCondIntMod = [elem for elem in uniqCondInt if elem]
    emptyPoz = [poz for poz, elem in enumerate(uniqCondInt) if not elem]
    for comb in itertools.product(*uniqCondIntMod):
        diffAB = max(comb) - min(comb)
        if bestABval is None or diffAB < bestABval:
            bestAB = comb
            bestABval = diffAB
        elif diffAB == bestABval:
            # consitions closer to N
            sumAB = sum(comb)
            bestSum = sum(bestAB)
            if sumAB < bestSum:
                bestAB = comb
                bestABval = diffAB
    if emptyPoz:
        newbestAB = []
        bpoz = 0
        for poz, elem in enumerate(uniqCondInt):
            if poz in emptyPoz:
                newbestAB.append(elem)
            else:
                newbestAB.append(bestAB[bpoz])
                bpoz += 1
    if not bestAB or emptyPoz:
        print("U", uniqCondInt, "BEST", bestAB, "N", newbestAB)
    return bestAB


def _isAromFormingRx(rxdata):
    uniqSbs, uniqProd = rxdata['smiles'].split('>>')
    if 'c' not in uniqProd:
        return False
    if 'c' in uniqProd and 'c' not in uniqSbs:
        return True
    if 'additionalInfo' not in rxdata:
        mappedSbs, mappedProd = rxdata['mappedRx'].split('>>')
    else:
        mappedSbs, mappedProd = rxdata['additionalInfo']['mappedRx'].split('>>')
    cSbs = re.findall(r'\[c:\d+\]', mappedSbs)
    cHsbs = re.findall(r'\[cH:\d+\]', mappedSbs)
    cProd = re.findall(r'\[c:\d+\]', mappedProd)
    cHprod = re.findall(r'\[cH:\d+\]', mappedProd)
    aromSbs = set()
    for gr in (cSbs, cHsbs):
        _ = [aromSbs.add(int(atm.split(':')[1][:-1])) for atm in gr]
    aromProd = set()
    for gr in (cProd, cHprod):
        _ = [aromProd.add(int(atm.split(':')[1][:-1])) for atm in gr]
    newAromInProd = aromProd - aromSbs
    if newAromInProd:
        return True
    return False


def _findProperRxtuple(sbs, prod, rxid, availTuples, args):
    suitableTuples = [tpl for tpl in availTuples if tpl[2] == rxid and prod == tpl[1]]
    if not suitableTuples and '.' in prod:
        for prd in set(prod.split('.')):
            suitableTuples = [tpl for tpl in availTuples if tpl[2] == rxid and prd == tpl[1]]
            if suitableTuples:
                break
    if args.debug:
        print("SBS", sbs, prod, rxid)
        print("SUI", suitableTuples, "AA", availTuples)
    if len(suitableTuples) == 1:
        return suitableTuples[0]
    print("==============")
    suitableTuples2 = []
    for tpl in suitableTuples:
        if set(tpl[0]).intersection(set(sbs)):
            suitableTuples2.append(tpl)
        # print("=====", tpl[0], type(tpl[0]), "SBS", type(sbs), sbs)
    if len(suitableTuples2) == 1:
        return suitableTuples2[0]
    print(f"?? 11: {sbs}  22: {prod} 33: {rxid} 444: {len(availTuples)}       555:  {suitableTuples} 222: {suitableTuples2}")
    raise


def _getRxPozToRmDueCond(data, commonData, rxInfo, args):
    acidBaseSet = set()
    tempSet = set()
    for variant in commonData['allowedCond']['minimalVariants']:
        for step in variant:
            ab, temp = step
            acidBaseSet.update(ab)
            tempSet.update(temp)
    if args.verbose > 5:
        print("MINIM", acidBaseSet, tempSet, "FROM", commonData['allowedCond'])
    rxToRm = []
    for rxpoz in data['idToFullRx']:
        rxids = data['idToFullRx'][rxpoz]['rxid']
        if isinstance(rxids, (str, int)):
            rxids = (rxids, )
        for rxid in rxids:
            if rxid in (1, ):
                continue
            try:
                rxab = set(rxInfo['rxes'][rxid]['rxacidbase'].replace('LA', '').replace('+', '').split('.')) - {''}
            except KeyError:
                print("RXIDS", rxid in rxInfo['rxes'], rxInfo['rxes'].keys())
                raise
            rxtemp = set(rxInfo['rxes'][rxid]['rxtemp'].upper().split('.'))
            if rxab:
                isABok = rxab.intersection(acidBaseSet)
            else:
                isABok = True
            if isABok and rxtemp.intersection(tempSet):
                if args.verbose > 15:
                    print("OK RXID", rxid, rxab, rxtemp)
            else:
                if args.verbose > 5:
                    print(f"RM RX due cond rx/path:: {rxab} {acidBaseSet} T: {rxtemp} {tempSet}")
                rxToRm.append(rxpoz)
    return rxToRm


def removeRxFromData(data, rxPozToRm):
    rxPozToRmSet = set(rxPozToRm)
    for rx in rxPozToRmSet:
        # print(f"Rmed:: {data['idToFullRx'][rx]}")
        del data['idToFullRx'][rx]
        # print("RXK1", tuple(data['idToFullRx'].keys()))
    # print("DATA", data.keys())
    rxtupleToRm = [rxt for rxt in data['rxtupleToId'] if data['rxtupleToId'][rxt] in rxPozToRmSet]
    # print("TO RM", rxtupleToRm)
    for rxt in rxtupleToRm:
        del data['rxtupleToId'][rxt]
    # print(f"toid: {data['rxtupleToId'].keys()}")
    # DATA  'rxToFullDict', 'graph', 'rxtupleToId'
    return data


def getWorkingConditions(allowedCond):
    variants = []
    for condOption in allowedCond['minimalVariants']:
        stepsList = [(0, condOption[0]), ]
        stepPrev = condOption[0]
        for poz, step in enumerate(condOption[1:], 1):
            if step == stepPrev:
                continue
            stepsList.append((poz, step))
            stepPrev = step
        variants.append(tuple(stepsList))
    return variants


def _addScalingForSideInOldConditions(data, args, commonData):
    mainPathRxSmiles = set([rxsmiAndRxid[0] for rxsmiAndRxid in commonData['allowedCond']['stepsVariants'][0]])
    if commonData['allowedCond']['minimalPozBeforeChangeToTarget']:
        forbidSbs = set(commonData['allowedCond']['minimalPozBeforeChangeToTarget'][0].values())
    else:
        forbidSbs = set()

    for graphId in data['idToFullRx']:
        rxSmi = data['idToFullRx'][graphId]['smiles']
        if rxSmi in mainPathRxSmiles:
            continue
        sbs = set(rxSmi.split('>>')[0].split('.'))
        if sbs.intersection(forbidSbs):
            data['idToFullRx'][graphId]['speedScale'] = args.scaleSideInOldCondition
    return data['idToFullRx']


def detectAllowedCond(data, rxinfo, mainPathInfo, args):
    topos = {'full': [], 'minimalVariants': [], 'stepsVariants': [], 'catUsage': [], 'conditionsOcc': []}
    topoOrders, rxgraph = _getAllowedTopoOrder(data)
    # doPrint = False
    if not topoOrders or (len(topoOrders) == 1 and not topoOrders[0]):
        # single step path
        #  [('CCCC(C(C(=O)OC)C(=O)OC)N1CCCCC1>>CCCC(C(C(=O)O)C(=O)O)N1CCCCC1', ('AW6431',))]
        print("MAINPATINF", mainPathInfo)
        rxid = mainPathInfo['rxes'][0][1][0]
        if args.verbose > 4:
            print("RXID", rxid)
            print("COND",  rxinfo['rxes'][rxid])
        # {'rxacidbase': 'B', 'rxspec': 'NaOH;;;KOH', 'alterCond': '', 'rxoxred': 'NA', 'rxfreewater': 'W',  'rxspeed': 'S', 'rxprotonic': 'P.AP', 'rxpolar': 'P.Non-P',
        # 'solvent': 'water', 'rxtemp': 'rt.H.VH', 'module': 'mech', 'rxidOryg': 'AW6431', 'rxidNew': 18090, 'fixedSb': '[OH-];;;2', 'whenUseFixed': '', 'intermol': 'Y',
        # 'rxclass': '2.1.5.ED', 'multistep': 'B', 'additionalReagents': '1..17.007', 'additionalReagentsForAlter': '1..17.007', 'rx_reactive': 'NA', 'needAlterSolvCond': ''}
        topos['stepsVariants'] = mainPathInfo['rxes']
        rxcondABoneRx = frozenset(rxinfo['rxes'][rxid]['rxacidbase'].split('.'))
        rxcondToneRx = frozenset([t.upper() for t in rxinfo['rxes'][rxid]['rxtemp'].split('.')])
        minVariant = [(rxcondABoneRx, rxcondToneRx), ]
        topos['minimalVariants'] = [minVariant, ]
        topos['minimalPozBeforeChangeToTarget'] = [{0: mainPathInfo['rxes'][0][0].split('>>')[1]}]
        topos['rxToGenSbs'] = {mainPathInfo['rxes'][0]: 1}
        topos['rxToGenProd'] = {mainPathInfo['rxes'][0]: 1}
        topos['conditionsOcc'] = [{'temp': {temp: [0, ] for temp in rxcondToneRx}, 'AB': {ab: [0, ] for ab in rxcondABoneRx}}, ]
        return topos, mainPathInfo
    for poz, topoOrder in enumerate(topoOrders):
        order = _getConditions(topoOrder, rxinfo)
        if not order:
            print(f"NO ORDER {order} FROM:: {topoOrders} Is not empty::: {bool(topoOrders)}")
        catalUsage = [poz for poz, t in enumerate(topoOrder) if t in mainPathInfo['hasCat']]
        minSet, minAB, minTemp = _getMinimalCondSet(order)
        condPopul = {'temp': dict(), 'AB': dict()}
        for poz, abset in enumerate(minAB):
            for ab in abset:
                if ab not in condPopul['AB']:
                    condPopul['AB'][ab] = []
                condPopul['AB'][ab].append(poz)
        for poz, tempset in enumerate(minTemp):
            for temp in tempset:
                if temp not in condPopul['temp']:
                    condPopul['temp'][temp] = []
                condPopul['temp'][temp].append(poz)
        topos['conditionsOcc'].append(condPopul)
        if args.includeConditions == 'minimalMild':
            minSet = _minimizeMinimalCond(minAB, minTemp, topoOrder, rxgraph)
        topos['full'].append(order)
        topos['catUsage'].append(catalUsage)
        topos['minimalVariants'].append(minSet)
        topos['stepsVariants'].append(topoOrder)
    topos['minimalPozBeforeChangeToTarget'] = [_getConditionChangesPos(topos['minimalVariants'][poz], topos['stepsVariants'][poz])
                                               for poz in range(len(topos['stepsVariants']))]
    if args.catAddTime == 'jitLate' and topos['catUsage'][0]:
        latest = [min(usage) for usage in topos['catUsage']]
        maxLatest = max(latest)
        latestCondPoz = [poz for poz, ltst in enumerate(latest) if ltst == maxLatest]
        for key in topos:
            topos[key] = [topoo for poz, topoo in enumerate(topos[key]) if poz in latestCondPoz]
        if args.debug:
            print("TOPOS", topos)
    topos['rxToGenSbs'], topos['rxToGenProd'] = _getGenForRxFromRxgraph(rxgraph)
    print("TOPOS", topos)
    return topos, mainPathInfo


def _getGenForRxFromRxgraph(rxgraph):
    initRxes = set([node for node in rxgraph.nodes if rxgraph.in_degree(node) == 0])
    target = [node for node in rxgraph.nodes if rxgraph.out_degree(node) == 0]
    if len(target) != 1:
        print("TARGETS:452", target)
        raise NotImplementedError
    target = target[0]
    rxToGenSbsDict = dict()
    rxToGenProdDict = dict()
    for rx in rxgraph.nodes:
        # rxToGenProdDict generation
        if rx == target:
            rxToGenProdDict[rx] = 0
        else:
            try:
                dist = networkx.shortest_path_length(rxgraph, source=rx, target=target)
            except:
                print("NO PATH FROM", rx, "TO", target, "IN NODES", rxgraph.nodes)
                raise
            # print("DIST", dist)
            rxToGenProdDict[rx] = dist

        # rxToGenSbsDict generation
        if rx in initRxes:
            rxToGenSbsDict[rx] = 1
        else:
            dists = []
            for initRx in initRxes:
                try:
                    dist = networkx.shortest_path_length(rxgraph, source=initRx, target=rx)
                    dists.append(dist)
                except networkx.NetworkXNoPath:
                    pass
            if not dists:
                raise NotImplementedError
            rxToGenSbsDict[rx] = 1 + min(dists)
    return rxToGenSbsDict, rxToGenProdDict


def _mapIdToFullRx(dataDict, args):
    retDict = dict()
    added = []
    for rxsma, rxid in dataDict['rxToFullDict']:
        sbs, prod = rxsma.split('>>')
        sbs = sbs.split('.')
        rxtuple = (tuple(sbs), prod, rxid)
        if rxtuple not in dataDict['rxtupleToId']:
            print("\n--- NF ", rxtuple, "---  NOT IN  --- ", dataDict['rxtupleToId'].keys())
            raise
        rxnodeID = dataDict['rxtupleToId'][rxtuple]
        retDict[rxnodeID] = dataDict['rxToFullDict'][(rxsma, rxid)]
        added.append(rxnodeID)
    for key in ['med', 'mech']:
        for rxlist in dataDict['competing'][key]:
            for rxinfo in rxlist:
                rxsma = rxinfo['fullRxSmiles']
                sbs, prod = rxsma.split('>>')
                sbs = sorted(sbs.split('.'))
                rxid = tuple(rxinfo['rxid'])
                rxtuple = (tuple(sbs), prod.split('.')[0], rxid)
                rxtuple1, rxtuple2 = rxtuple, None
                if rxtuple not in dataDict['rxtupleToId']:
                    rxtuple = _findProperRxtuple(sbs, prod, rxid, tuple(dataDict['rxtupleToId'].keys()), args)
                    # continue
                    rxtuple2 = rxtuple
                rxnodeID = dataDict['rxtupleToId'][rxtuple]
                if rxnodeID == 16:
                    print(f"1616161616  {rxinfo} 11: {rxtuple1} 222 {rxtuple2}")
                added.append(rxnodeID)
                retDict[rxnodeID] = rxinfo
    if len(retDict) != len(dataDict['rxtupleToId']):
        k1 = set(retDict.keys())
        k2 = set(dataDict['rxtupleToId'].values())
        print("SET DIFF", k1-k2, k2-k1, "ADDED", added)
        print("WRONG1:", len(retDict), len(dataDict['rxToFullDict']), retDict.keys())
        print("WRONG2:", len(dataDict['rxtupleToId']), dataDict['rxtupleToId'].values())
        print([(key, retDict[key]['fullRxSmiles'], retDict[key]['rxid']) for key in retDict])
        print("===", [(key, dataDict['rxtupleToId'][key]) for key in dataDict['rxtupleToId']])
        print("DDD", dataDict.keys())
        raise
    return retDict


def _getConditionChangesPos(conds, rxes):
    changesPos = dict()
    for poz, cnd in enumerate(conds[1:], 1):
        if cnd != conds[poz-1]:
            rxsmi, _ = rxes[poz-1]
            _, target = rxsmi.split('>>')
            if '.' in target:
                raise NotImplementedError
            changesPos[poz-1] = target
    return changesPos


def _getAllowedTopoOrder(data):
    digraph = networkx.DiGraph()
    for rxsmaRxid in data['data']['rxToFullDict']:
        rxsma, _ = rxsmaRxid
        sbses, prods = rxsma.split('>>')
        for sbs in sbses.split('.'):
            digraph.add_edge(sbs, rxsmaRxid)
        for prd in prods.split('.'):
            digraph.add_edge(rxsmaRxid, prd)
    finalrx = tuple(digraph.predecessors(data['data']['targetSmiles']))
    assert len(finalrx) == 1
    finalrx = finalrx[0]
    rxAndCmds = {finalrx: tuple(digraph.predecessors(finalrx))}
    rxgraph = networkx.DiGraph()
    while rxAndCmds:
        nextRxAndCmds = dict()
        for nextrx in rxAndCmds:
            for prod in rxAndCmds[nextrx]:
                prevrx = tuple(digraph.predecessors(prod))
                if not prevrx:
                    continue
                if len(prevrx) != 1:
                    print("PROD", prod, "PREB", prevrx)
                    raise
                prevrx = prevrx[0]
                sbs = tuple(digraph.predecessors(prevrx))
                rxgraph.add_edge(prevrx, nextrx)
                nextRxAndCmds[prevrx] = sbs
        rxAndCmds = nextRxAndCmds
    topord = tuple(networkx.all_topological_sorts(rxgraph))
    # reversed
    # topord = tuple([tuple(reversed(sort)) for sort in topord])
    return topord, rxgraph


def _getMinimalCondSet(order):
    ab = [_formAB(cnd[0]) for cnd in order]
    temp = [_formTemp(cnd[1]) for cnd in order]
    if not ab:
        print("CANNOT FORM AB FROM", order)
        print(f"AB:: {ab} temp:: {temp}")
        raise
    if not temp:
        print("CANNOT FORM TEMP FROM", order)
        print(f"AB:: {ab} temp:: {temp}")
        raise
    minab = _getMinimSet(ab)
    mintemp = _getMinimSet(temp)
    mintopo = list((frozenset(ab), frozenset(temp)) for ab, temp in zip(minab, mintemp))
    return mintopo, minab, mintemp


def _formAB(cond):
    cond = cond.replace('LA', '').replace('+', '')
    cond = set([cnd for cnd in cond.split('.') if cnd])
    return cond


def _formTemp(cond):
    cnds = set([cnd.upper() for cnd in cond.split('.')])
    return cnds


def _getMinimSet(cond):
    intersec = set.intersection(*cond)
    if intersec:
        return [intersec for _ in cond]
    # greedy set
    resList = []
    curSet = cond[0]
    length = 1
    for poz, cset in enumerate(cond[1:]):
        intersec = curSet.intersection(cset)
        if intersec:
            curSet = intersec
            length += 1
        else:
            for _ in range(length):
                resList.append(curSet)
            curSet = cset
            length = 1
    if curSet:
        for _ in range(length):
            resList.append(curSet)
    if len(resList) != len(cond):
        print("COND", cond)
        print("RES", resList)
        raise
    return resList


def _getConditions(topoOrder, rxinfo):
    conds = []
    for rxsmaRxid in topoOrder:
        _, rxid = rxsmaRxid
        assert len(rxid) == 1
        rxid = rxid[0]
        if rxid != 1:
            rxinf = rxinfo['rxes'][rxid]
            conds.append((rxinf['rxacidbase'], rxinf['rxtemp']))
    return conds


def getMainPathInfo(data, rxInfo, args):
    # rxinfo {
    mainPathInfo = {'isWaterSbs': False, 'isWaterCondNeeded': False, 'reactiveSubstrates': dict(), 'rxes': [], 'byprods': rxInfo['byprods'],
                    'redox': set(), 'hasCat': dict(), 'initSubstrates': set(), 'smilesToGen': dict(), 'graph': None}
    # checkedCmd = set()
    print("DATA", data)
    # dg = data['graph']
    allSbses = set()
    allProds = set()
    graph = networkx.DiGraph()
    for rxsmaRxid in data['data']['rxToFullDict']:
        rxsma, rxid = rxsmaRxid
        fullSbses, fullProds = data['data']['rxToFullDict'][rxsmaRxid]['additionalInfo']['mappedRx'].split('>>')
        fullSbses = [getUniqSmiles(smi) for smi in fullSbses.split('.')]
        fullProds = [getUniqSmiles(smi) for smi in fullProds.split('.')]
        for sbs in fullSbses:
            graph.add_edge(sbs, rxsmaRxid)
        for prd in fullProds:
            graph.add_edge(rxsmaRxid, prd)
        sbses, prods = rxsma.split('>>')
        sbses = sbses.split('.')
        sbsesSet = set(sbses)
        allSbses.update(sbsesSet)
        allProds.update(set(prods.split('.')))
        mainPathInfo['rxes'].append(rxsmaRxid)
        if 'O' in sbses:
            mainPathInfo['isWaterSbs'] = True
        assert len(rxid) == 1
        # rxids = tuple(rxid)
        rxid = rxid[0]
        if rxid == 1:
            continue
        rxinfRedox = rxInfo['rxes'][rxid]["rxoxred"].upper().split('.')
        if 'RED' in rxinfRedox:
            mainPathInfo['redox'].add('RED')
            print("REDOX", rxinfRedox)
        if 'OX' in rxinfRedox:
            mainPathInfo['redox'].add('OX')
            print("REDOX", rxinfRedox)
        rxinfWF = rxInfo['rxes'][rxid]['rxfreewater']
        if rxinfWF == 'W':
            mainPathInfo['isWaterCondNeeded'] = True

        initSbs = allSbses - allProds
        mainPathInfo['initSubstrates'] = initSbs
        mainPathInfo['graph'] = graph

        if data['data']['conditions'].get('catalysts', []):
            mainPathInfo['hasCat'] = tuple(data['data']['conditions'])
    # sbsAndCat = initSbs.union(cats)
    # mainPathInfo['smilesToGen'] = loaderGraph.assignSmilesToGenInGraph(graph, sbsAndCat)
    return mainPathInfo


def _detectFixed(rxToFullDict):
    # dataDict[n][0]['data']['rxToFullDict']
    fixedDict = dict()
    for rx in rxToFullDict:
        rxinfo = rxToFullDict[rx]
        sbs, prod = rxinfo['smiles'].split('>>')
        print("RX", rxinfo)
        fullsbs, fullprods = rxinfo['additionalInfo']['mappedRx'].split('>>')
        sbses = sbs.split('.')
        fullsbses = [getUniqSmiles(sbs) for sbs in fullsbs.split('.')]
        if len(sbses) == len(fullsbses):
            continue
        fixed = [s for s in fullsbses if s not in sbs]
        fixedDict[rx] = fixed
    return fixedDict


def _detectAdditives(data):
    maxgen = max(data.keys())
    retDict = dict()
    for rx in data[maxgen]['data']['rxToFullDict']:
        addits = set()
        # if data[maxgen]['data']['rxToFullDict'][rx].get('byproducts'):
        #     pri nt("BYPRODS", data[maxgen]['data']['rxToFullDict'][rx]['byproducts'])
        try:
            additives = data[maxgen]['data']['rxToFullDict'][rx]['byproducts']['additives']
            for additeList in additives:
                for addit in additeList:
                    if addit and addit != '-':
                        canon = Chem.CanonSmiles(addit)
                        addits.add(canon)
        except:
            print("PROBLEM WITH gen:", maxgen,  "rx: ", rx)
            print('missed byprods entry in json in rx',  data[maxgen]['data']['rxToFullDict'][rx])
            raise
        if addits:
            # ('C=O.CC(=O)C[N+](=O)[O-]>>C=C(C(C)=O)[N+](=O)[O-]', ('AW7761',))
            sbs, prod = rx[0].split('>>')
            rxtuple = (tuple(sorted(sbs.split('.'))), prod, rx[1])
            rxidgraph = data[maxgen]['data']['rxtupleToId'][rxtuple]
            retDict[rxidgraph] = addits
    return retDict


# #####################


def getUniqSmiles(molOrSmiles):
    if isinstance(molOrSmiles, str):
        smi = molOrSmiles
    else:
        smi = Chem.MolToSmiles(molOrSmiles)

    nonums = re.sub(':[0-9]*]', ']', smi)
    try:
        smi2 = Chem.CanonSmiles(nonums)
    except:
        print("WRONG SMI", smi)
        raise
    return smi2


# ###################################################


def buildDicts(args):
    jsonData = json.load(open(args.input))
    baseDict = {'competing': {'med': [], 'mech': []},
                'conditions':  {'catalysts': jsonData.get('catalysts', []),
                                'isOnepot': jsonData.get('isOnepot', False),
                                'formAB': jsonData.get('formAB', {}),
                                'initC': jsonData.get('initC', {})
                                },
                'graph': [tuple(edge) for edge in jsonData['graph']],
                'targetSmiles': jsonData['targetSmiles']}
    rxToFullDict, rxtupleToId = formRxinfo(jsonData)
    baseDict['rxToFullDict'] = rxToFullDict
    baseDict['rxtupleToId'] = rxtupleToId
    baseDict['allSmiles'], graph = getAllSmiles(jsonData['graph'])
    dictL2 = copy.deepcopy(baseDict)
    dictL3 = copy.deepcopy(baseDict)
    dictL4 = copy.deepcopy(baseDict)
    if "graphL3" in jsonData:
        dictL3, graph = addL(dictL3, jsonData['graphL3'], jsonData['rxesL3'], jsonData.get('formABL3', dict()), graph)
        if "graphL4" in jsonData:
            toCopy = ('graph', 'rxToFullDict', 'rxtupleToId')
            for elem in toCopy:
                dictL4[elem] = copy.deepcopy(dictL3[elem])
            dictL4, graph = addL(dictL4, jsonData['graphL4'], jsonData['rxesL4'], jsonData.get('formABL4', dict()), graph)
        else:
            dictL4 = copy.deepcopy(dictL3)
    return dictL2, dictL3, dictL4


def addL(dct, edges, rxes, formAB, pathGraph):
    for edge in edges:
        dct['graph'].append(tuple(edge))
    for rxnum in rxes:
        rxdict = rxes[rxnum]
        rxid = rxdict['rxid']
        rxhash = (rxdict['smiles'], (rxid, ))
        mapSbs, mapProd = rxdict['mappedRx'].split('>>')
        noMapSbs = [getUniqSmiles(smi) for smi in mapSbs.split('.')]
        noMapProds = [getUniqSmiles(smi) for smi in mapProd.split('.')]
        fullNoMap = '.'.join(noMapSbs) + '>>' + '.'.join(noMapProds)
        dct['rxToFullDict'][rxhash] = {'smiles': rxdict['smiles'], 'rxid': rxid, 'additionalInfo': {'mappedRx': rxdict['mappedRx'], },
                                       'byproducts': rxdict['byproducts'], 'fullRxSmiles': fullNoMap}
        if 'DE' in rxdict:
            dct['rxToFullDict'][rxhash]['DE'] = rxdict['DE']
        if 'rate' in rxdict:
            dct['rxToFullDict'][rxhash]['rate'] = rxdict['rate']
        # (('FC(F)(F)c1[c-]cc2c(c1)OCO2', 'IC#Cc1ccccc1'), 'FC(F)(F)c1cc2c(cc1I)OCO2', ('AW7188',)): 1,
        sbs, prod = rxdict['smiles'].split('>>')
        rxtuple = (tuple(sbs.split('.')), prod, (rxid, ))
        dct['rxtupleToId'][rxtuple] = rxnum
    dct['allSmiles'], pathGraph = appendSmiles(dct['allSmiles'], pathGraph, edges)
    if formAB:
        dct['formAB'].update(formAB)
    return dct, pathGraph


def appendSmiles(allSmiles, pathGraph, edges):
    for edge in edges:
        pathGraph.add_edge(*edge)
    pathRxes = [n for n in pathGraph.nodes if isinstance(n, int)]
    pathSmiles = [n for n in pathGraph.nodes if not isinstance(n, int)]
    seenrx = set()
    while len(allSmiles) < len(pathSmiles):
        for rx in pathRxes:
            if rx in seenrx:
                continue
            sbses = [smi for smi in pathGraph.pred[rx]]
            if all(cmd in allSmiles for cmd in sbses):
                gen = 1 + max([allSmiles[smi] for smi in sbses])
                for cmd in pathGraph.succ[rx]:
                    if cmd not in allSmiles:
                        allSmiles[cmd] = gen
                seenrx.add(rx)
    return allSmiles, pathGraph
# ########


def formRxinfo(jsonData):
    rxToFullDict = dict()
    rxtupleToId = dict()
    for rxnum in jsonData['rxes']:
        rxdict = jsonData['rxes'][rxnum]
        rxid = rxdict['rxid']
        rxhash = (rxdict['smiles'], (rxid, ))
        mapSbs, mapProd = rxdict['mappedRx'].split('>>')
        noMapSbs = [getUniqSmiles(smi) for smi in mapSbs.split('.')]
        noMapProds = [getUniqSmiles(smi) for smi in mapProd.split('.')]
        fullNoMap = '.'.join(noMapSbs) + '>>' + '.'.join(noMapProds)
        rxToFullDict[rxhash] = {'smiles': rxdict['smiles'], 'rxid': rxid,
                                'additionalInfo': {'mappedRx': rxdict['mappedRx'], },
                                'byproducts': rxdict['byproducts'], 'fullRxSmiles': fullNoMap}
        # (('FC(F)(F)c1[c-]cc2c(c1)OCO2', 'IC#Cc1ccccc1'), 'FC(F)(F)c1cc2c(cc1I)OCO2', ('AW7188',)): 1,
        if 'DE' in rxdict:
            rxToFullDict[rxhash]['DE'] = rxdict['DE']
        if 'rate' in rxdict:
            rxToFullDict[rxhash]['rate'] = rxdict['rate']
        sbs, prod = rxdict['smiles'].split('>>')
        rxtuple = (tuple(sbs.split('.')), prod, (rxid, ))
        rxtupleToId[rxtuple] = rxnum
    return rxToFullDict, rxtupleToId


def getAllSmiles(edges):
    graph = networkx.DiGraph()
    for edge in edges:
        graph.add_edge(edge[0], edge[1])
    allSmiles = [n for n in graph.nodes if not isinstance(n, int)]
    allRxes = [n for n in graph.nodes if isinstance(n, int)]
    allGens = dict()
    for smi in allSmiles:
        if not graph.pred[smi]:
            allGens[smi] = 0
    seenrx = set()
    while len(allGens) < len(allSmiles):
        for rx in allRxes:
            if rx in seenrx:
                continue
            sbses = [smi for smi in graph.pred[rx]]
            if all(cmd in allGens for cmd in sbses):
                gen = 1 + max([allGens[smi] for smi in sbses])
                for cmd in graph.succ[rx]:
                    if cmd not in allGens:
                        allGens[cmd] = gen
                seenrx.add(rx)
    return allGens, graph


if __name__ == "__main__":
    args = parseArgs()
    main(args)
