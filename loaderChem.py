from rdkit import Chem
import helpers #, mayr


def addReactionInfos(dataDict, mayrData, adbData, mutexDict, args):
    cache = dict()
    # hasProblems = False
    # values = []
    for name in dataDict:
        commonData = dataDict[name][0]['common']
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
        if args.debug:
            print("MUTEX:3409:", name, mutexToRxOnMainPath, "MP", mainPathRxid)
        for gen in dataDict[name]:
            for rxid in dataDict[name][gen]['data']['idToFullRx']:
                rxsmiles = dataDict[name][gen]['data']['idToFullRx'][rxid]['fullRxSmiles']
                if rxsmiles not in cache:
                    try:
                        mappedSmiles = dataDict[name][gen]['data']['idToFullRx'][rxid]['additionalInfo']['mappedRx']
                    except:
                        mappedSmiles = dataDict[name][gen]['data']['idToFullRx'][rxid]['mappedRx']
                    # fix problem with reaction with reverse order: products>>substrates
                    mappedSmiles = _correctMappedSmiles(mappedSmiles, rxsmiles)
                    rxidOryg = dataDict[name][gen]['data']['idToFullRx'][rxid]['rxid']
                    assert len(set(rxidOryg)) == 1
                    if rxidOryg[0] == 1:
                        # condition changes
                        continue
                    rxInfo = adbData['rxes'][rxidOryg[0]]
                    tryCalcMayr = True
                    if is1c(rxInfo):
                        if args.debug > 100:
                            print(f"Mayr: no rate for 1c: {rxsmiles} {rxidOryg}")
                        tryCalcMayr = False
                    if not _isAddSub(rxInfo['rxclass']):
                        if args.debug > 100:
                            print(f"Mayr: no rate for non AD/SUB: {rxsmiles} {rxidOryg}")
                        tryCalcMayr = False
                    try:
                        mayrRate = None
                        if tryCalcMayr:
                            mayrRate = mayr.calcMayrRate(rxsmiles, mappedSmiles, mayrData, rxInfo, args)
                        ringBoostFactor = _calcRingBoost(rxsmiles, mappedSmiles, rxInfo, args)
                        fullrx = dataDict[name][gen]['data']['idToFullRx'][rxid]['fullRxSmiles']
                        fullsbs, fullprod = fullrx.split('>>')
                        # if fullsbs.count('.') > 1:
                        # print("FFF", rxidOryg, fullsbs.upper().count('C'), fullprod.upper().count('C'), fullrx, rxsmiles)
                        stoichScale = _calcStoichScale(rxsmiles, mappedSmiles, args)

                        isMutex = set(rxidOryg).intersection(mutexToRxOnMainPath)
                        cache[rxsmiles] = {'mayr': mayrRate, 'boostRing': ringBoostFactor, 'stoichScale': stoichScale, 'isMutex': isMutex}
                    except:
                        print(f"\n RXID {rxid} {rxsmiles} \n {dataDict[name][gen]['data']['idToFullRx'][rxid]}")
                        raise
                    # values.append(cache[rxsmiles])
                    if args.debug > 100:
                        print(f"Mayr RXVAL:: {cache[rxsmiles]} {rxsmiles} {rxidOryg}")
                if cache[rxsmiles]['mayr'] is not None:
                    dataDict[name][gen]['data']['idToFullRx'][rxid]['rate'] = cache[rxsmiles]['mayr']
                if cache[rxsmiles]['boostRing'] is not None:
                    dataDict[name][gen]['data']['idToFullRx'][rxid]['boostRing'] = cache[rxsmiles]['boostRing']
                if cache[rxsmiles]['stoichScale'] is not None:
                    dataDict[name][gen]['data']['idToFullRx'][rxid]['stoichScale'] = cache[rxsmiles]['stoichScale']
                if cache[rxsmiles]['isMutex']:
                    if args.debug:
                        print(f"MUTEX:3455:: {name} {gen} {rxid} {dataDict[name][gen]['data']['idToFullRx'][rxid]['rxid']}")
                    dataDict[name][gen]['data']['idToFullRx'][rxid]['isMutex'] = cache[rxsmiles]['isMutex']
                # print(f"RATEmayr {name} {gen}  {rxsmiles} ::  {cache[rxsmiles]}")
    return dataDict


def _correctMappedSmiles(mappedSmiles, rxsmiles):
    sbs, prod = rxsmiles.split('>>')
    sbsP, prodP = sbs.count('+'), prod.count('+')
    mappedSbs, mappedProd = mappedSmiles.split('>>')
    if prodP != sbsP:
        mappedSbsP, mappedProdP = mappedSbs.count('+'), mappedProd.count('+')
        if sbsP != mappedSbsP or prodP != mappedProdP:
            # print("WE HAVE PROBLEM", rxsmiles, mappedSmiles)
            if sbsP == mappedProdP and prodP == mappedSbsP:
                # print("REVERSE::", rxsmiles, "MAPPEDBEFOR", mappedSmiles)
                return mappedProd + '>>' + mappedSbs
            else:
                print("PROBLEM NOT FIXED", rxsmiles)
        return mappedSmiles
    uniqSbs = [helpers.getUniqSmiles(Chem.MolFromSmiles(smi)) for smi in mappedSbs.split('.')]
    sbs = sbs.split('.')
    prod = prod.split('.')
    if set(uniqSbs) == set(sbs):
        return mappedSmiles
    # print("REVERSE::", rxsmiles, "MAPPEDBEFOR", mappedSmiles)
    return mappedProd + '>>' + mappedSbs


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


def _isAddSub(rxtypes):
    for rxtype in rxtypes.split(';;;'):
        if rxtype.endswith('.C') or rxtype.endswith('.E') or rxtype.endswith('.ED'):
            return True
    return False


def is1c(rxInfo):
    sbs, _ = rxInfo['rxSmarts'].split('>>')
    if len(rxInfo['reactingAtomsIdx']) == 1 or '.' not in sbs:
        return True
    if len(rxInfo['reactingAtomsIdx']) > 1 and not (rxInfo['same'] or rxInfo['fixedSb']):
        return False
    return False


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
