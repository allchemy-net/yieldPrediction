import re
import networkx
import loaderHelper


def _buildGraph(edgesList, additives, data, rxInfo, longShortCut, commonData, args):
    graph = networkx.DiGraph()
    for edge in edgesList:
        if not args.includeAdditives and edge[0] in additives and edge[1] in additives[edge[0]]:
            if args.debug:
                print("IGNORE", edge, additives)
            continue
        graph.add_edge(*edge)
        if args.debug > 100:
            print("ADDED EDGE", edge)

    rxPozToRm = []
    if args.disableLongCut > 0:
        rxPozToRmDueToLongcut = _getRxPozOfLongCutToRm(data, commonData, longShortCut, args)
        rxPozToRm.extend(rxPozToRmDueToLongcut)
        if args.debug:
            print(f"RM due to:longcut: {rxPozToRm}")
    if args.removeDuplicates != 'N':
        if args.removeDuplicates == 'sbsesProd':
            rxPozToRmDueDuplicates = _getRxPozDuplicatesToRm(data, commonData, rxInfo, args)
            rxPozToRm.extend(rxPozToRmDueDuplicates)
            if args.debug:
                print(f"RM due to:dupicate: {rxPozToRm}")
        else:
            raise NotImplementedError
    if args.includeConditions != 'no':
        if args.includeConditions in {'minimal', 'minimalMild'}:
            rxPozToRmDueToCond = _getRxPozToRmDueCond(data, commonData, rxInfo, args)
            rxPozToRm.extend(rxPozToRmDueToCond)
            if args.debug:
                print(f"RM due to:conditions: {rxPozToRm}")
        else:
            raise NotImplementedError
    if args.rmRxMode != 'N':
        rxToRm = _getRxPozToRmDueToWildLA(data, commonData, rxInfo, args)
        if args.debug:
            print("RXTORM", rxToRm)
        rxPozToRm.extend(rxToRm)
    if args.removeLateSides != 'N':
        rxToRmDueToReactiveSbs = _getLateSideToRm(graph, data, commonData, args)
        if args.debug:
            print(f"rm rx due to ractive sbs {rxToRmDueToReactiveSbs}")
        rxPozToRm.extend(rxToRmDueToReactiveSbs)
    if args.redox in {'rm', 'rmRedOxWater'}:
        rxToRmDueToRedox = _getRedoxToRm(graph, data, commonData, rxInfo, args)
        if args.debug:
            print("RX TO REDOX", rxToRmDueToRedox, "XX", rxPozToRm)
        rxPozToRm.extend(rxToRmDueToRedox)
    if rxPozToRm:
        initSubstrates = commonData['allInitSbs']
        graph, removedRxPoz = loaderHelper.removeRxFromGraph(graph, rxPozToRm, initSubstrates, args)
        if args.debug:
            print("RX TO rm", removedRxPoz, "XX", rxPozToRm)
        if args.debug > 8:
            print(f"REMOVED_RX_FROM_GR:50 {rxPozToRm} RM: {[(data['idToFullRx'][rxpoz]['smiles'], data['idToFullRx'][rxpoz]['rxid']) for rxpoz in rxPozToRm]}")
            print("RXIDS", [data['idToFullRx'][rxpoz]['rxid'] for rxpoz in rxPozToRm])
        data = loaderHelper.removeRxFromData(data, removedRxPoz)
    # check add missing side reagents
    newDataIdToFullRx = dict()
    for rxid in data['idToFullRx']:
        try:
            sbses = set(graph.predecessors(rxid))
            prods = set(graph.successors(rxid))
        except networkx.exception.NetworkXError:
            if args.disableLongCut == 0:
                raise NotImplementedError
            if args.debug:
                print(f"rxpoz {rxid} NOT IN GRAPH")
            continue
        rsbses, rprods = data['idToFullRx'][rxid]['fullRxSmiles'].split('>>')
        isReversible = False
        if args.addReverse == 'Y':
            isReversible = set(data['idToFullRx'][rxid]['rxid']).intersection(rxInfo['reversibleRx'])
        elif args.addReverse == 'YnoTAB':
            isReversible = set(data['idToFullRx'][rxid]['rxid']).intersection(rxInfo['reversibleRx'])
            if loaderHelper._isWildFormingRx(data['idToFullRx'][rxid]):
                if args.debug > 100:
                    print("is wild forming", rxid, data['idToFullRx'][rxid])
                isReversible = False
            if loaderHelper._isTautomerisationRx(data['idToFullRx'][rxid]):
                if args.debug > 100:
                    print("is tautomerisation", rxid, data['idToFullRx'][rxid]['rxid'], data['idToFullRx'][rxid]['smiles'])
                isReversible = False
        elif args.addReverse == 'wildToMainOnly':
            isReversible = loaderHelper._isWildFormingRx(data['idToFullRx'][rxid])
        if isReversible:
            if _hasGraphReverseRx(graph, data['idToFullRx'][rxid], data['idToFullRx']):
                if args.debug > 100:
                    print("has reversible", rxid, data['idToFullRx'][rxid])
                isReversible = False
            if args.filterReverse == 'arom':
                if _isAromFormingRx(data['idToFullRx'][rxid]):
                    isReversible = False
        rprods = set(rprods.split('.'))
        rsbses = set(rsbses.split('.'))
        if isReversible:
            rxidRev = f'{rxid}_REV'
            # if rxidRev == '38_REV':
            #    pri nt("::::::::", rxidRev, "S", rsbses, "P", rprods, "EL", len(edgesList))
            for sbs in rsbses:
                graph.add_edge(rxidRev, sbs)
            for prod in rprods:
                graph.add_edge(prod, rxidRev)
            newDataIdToFullRx[rxidRev] = _formReverseRx(data['idToFullRx'][rxid])
            print("ADD REVERSIBLE", rxidRev, newDataIdToFullRx[rxidRev])
        if sbses == rsbses and prods == rprods:
            if args.debug > 100:
                print(f">== ignore (reaction with exactly the same prod and sbs are in graph):: {rxid} rsbses: {rsbses} {sbses} prods:: {rprods} {prods}")
                print(">== FULL", data['idToFullRx'][rxid])
            continue
        for sbs in rsbses - sbses:
            if args.debug > 10:
                print("ADD edge to side from sbs", sbs, rxid)
            graph.add_edge(sbs, rxid)
        for prod in rprods - prods:
            if args.debug > 10:
                print("ADD edge to side to prod", rxid, prod)
            graph.add_edge(rxid, prod)
        if args.debug:
            print("==>", rxid, sbses, prods, data['idToFullRx'][rxid]['fullRxSmiles'])
    print("NEW DATA TO FULLRX KEYS::", newDataIdToFullRx.keys())
    if newDataIdToFullRx:
        data['idToFullRx'].update(newDataIdToFullRx)
    if args.debug > 100:
        print(f"FULLGRAPH:: {graph.edges}")
    print("DATA IDTOFULLRX", data['idToFullRx'].keys())
    return graph, data


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


def _getRedoxToRm(graph, data, commonData, rxInfo, args):
    # args.redox in {'rm', 'rmRedOxWater'}
    redoxInfo = commonData['mainPathInfo']['redox']
    mainPathRxes = set()
    rxPozToRm = []
    for smarx in commonData['mainPathInfo']['redox']:
        _ = [mainPathRxes.add(rx) for rx in smarx[1]]
    for rxpoz in data['idToFullRx']:
        for rxnum in data['idToFullRx'][rxpoz]['rxid']:
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


def _getRxPozOfLongCutToRm(data, commonData, longShortCut, args):
    # longShortCut
    longcutSmiles = dict()
    shortcutSmiles = dict()
    for rxpoz in data['idToFullRx']:
        rxnums = data['idToFullRx'][rxpoz]['rxid']
        if len(rxnums) != 1:
            raise NotImplementedError
        rxnum = rxnums[0]
        if rxnum in longShortCut['longRxesDict']:
            longShortGroupIds = longShortCut['longRxesDict'][rxnum]
            sbses = data['idToFullRx'][rxpoz]['smiles'].split('>>')[0].split('.')
            if args.debug > 5:
                print("LONGCUT", rxnum, sbses, data['idToFullRx'][rxpoz]['smiles'])
            longcutSmiles = _addSbsesToDictSmiles(sbses, longcutSmiles, longShortGroupIds, rxnum, rxpoz)
        if rxnum in longShortCut['shortRxesDict']:
            longShortGroupIds = longShortCut['shortRxesDict'][rxnum]
            sbses = data['idToFullRx'][rxpoz]['smiles'].split('>>')[0].split('.')
            if args.debug > 5:
                print("Short", rxnum, sbses, data['idToFullRx'][rxpoz]['smiles'])
            shortcutSmiles = _addSbsesToDictSmiles(sbses, shortcutSmiles, longShortGroupIds, rxnum, rxpoz)
    # analyze smiles which are in both longcut and shortcut
    rxPozToRm = []
    for smiles in set(longcutSmiles.keys()).intersection(set(shortcutSmiles.keys())):
        if args.debug:
            print("LONGSHORT", smiles, "LONG", longcutSmiles[smiles], "SHORT:", shortcutSmiles[smiles])
            print("COMMON DATA", commonData['mainPathInfo'].keys())
        commonGroups = set(longcutSmiles[smiles].keys()).intersection(set(shortcutSmiles[smiles].keys()))
        mainPathRxesSet = set([rxid[0] for smirx, rxid in commonData['mainPathInfo']['rxes']])
        for group in commonGroups:
            for rxinfo in longcutSmiles[smiles][group]:
                rxpoz = rxinfo[0]
                sbs = set(rxinfo[2])
                matchShort = []
                for shortInfo in shortcutSmiles[smiles][group]:
                    shortSbs = set(shortInfo[2])
                    if sbs.issubset(shortSbs):
                        matchShort.append(shortInfo)
                if matchShort:
                    doRemove = True
                    if rxinfo[1] in mainPathRxesSet:
                        suspRxes = [rx for rx in commonData['mainPathInfo']['rxes'] if rxinfo[1] in rx[1]]
                        for suspRx in suspRxes:
                            rxsma, _ = suspRx
                            mainRxSbs, _ = rxsma.split('>>')
                            mainRxSbs = set(mainRxSbs.split('.'))
                            thisSbses = set(rxinfo[2])
                            if thisSbses == mainRxSbs:
                                doRemove = False
                                if args.debug > 10:
                                    print("dont remove", rxinfo, "MAINRX", commonData['mainPathInfo']['rxes'])
                                break
                    if doRemove:
                        rxPozToRm.append(rxpoz)
                        if args.debug:
                            print("REMOVE NODE", rxpoz, smiles, "LONG", rxinfo[1:3], "SHORT",  matchShort)
    if args.debug:
        print("REMOVE RX to rm due to longcut", rxPozToRm)
    return rxPozToRm


def _getRxPozToRmDueCond(data, commonData, rxInfo, args):
    acidBaseSet = set()
    tempSet = set()
    for variant in commonData['allowedCond']['minimalVariants']:
        for step in variant:
            ab, temp = step
            acidBaseSet.update(ab)
            tempSet.update(temp)
    if args.debug > 5:
        print("MINIM", acidBaseSet, tempSet, "FROM", commonData['allowedCond'])
    rxToRm = []
    for rxpoz in data['idToFullRx']:
        for rxid in data['idToFullRx'][rxpoz]['rxid']:
            if rxid in (1, ):
                continue
            try:
                rxab = set(rxInfo['rxes'][rxid]['rxacidbase'].replace('LA', '').replace('+', '').split('.')) - {''}
            except KeyError:
                print("RXIDS", rxInfo['rxes'].keys())
                raise
            rxtemp = set(rxInfo['rxes'][rxid]['rxtemp'].upper().split('.'))
            if rxab:
                isABok = rxab.intersection(acidBaseSet)
            else:
                isABok = True
            if isABok and rxtemp.intersection(tempSet):
                if args.debug > 15:
                    print("OK RXID", rxid, rxab, rxtemp)
            else:
                if args.debug > 5:
                    print(f"RM RX due cond rx/path:: {rxab} {acidBaseSet} T: {rxtemp} {tempSet}")
                rxToRm.append(rxpoz)
    return rxToRm

def _getSbsFSandProd(rxsma):
    sbs, prod = rxsma.split('>>')
    sbs = frozenset(sbs.split('.'))
    prod = prod.split('.')
    assert len(prod) == 1
    prod = prod[0]
    return (sbs, prod)


def _getRxPozDuplicatesToRm(data, commonData, adbData, args):
    rxToRm = []
    seenRxes = dict()
    for rxpoz in data['idToFullRx']:
        rxsbs, rxprod = data['idToFullRx'][rxpoz]['smiles'].split('>>')
        if '.' in rxprod:
            raise NotImplementedError
        rxsbs = frozenset(rxsbs.split('.'))
        rxhash = (rxsbs, rxprod)
        if rxhash not in seenRxes:
            seenRxes[rxhash] = list()
        seenRxes[rxhash].append((rxpoz, data['idToFullRx'][rxpoz]))
            # rxToRm.append(rxpoz)
    # print("CCC", commonData['allowedCond']['conditionsOcc'], commonData.keys())
    mainPath = {_getSbsFSandProd(rxinfo[0]): rxinfo[1] for rxinfo in commonData['mainPathInfo']['rxes']}
    for rx in seenRxes:
        if len(seenRxes[rx]) == 1:
            continue
        duplicateList = seenRxes[rx]
        if rx in mainPath:
            mainPathRx, otherRx = _getMainPathAndOtherRxpoz(duplicateList, mainPath[rx])
            if mainPathRx:
                rxToRm.extend(otherRx)
                continue
        shortCuts = [rxinfo[1]['isShortcut'] for rxinfo in duplicateList]
        if any(shortCuts):
            # print("DIFF", rx, seenRxes[rx], "\nCOMMON:", mainPath)
            newDuplicateList = []
            for poz, isShortCut in enumerate(shortCuts):
                if not isShortCut:
                    nodeToRm = duplicateList[poz][0]
                    if args.debug > 2:
                        print("node to rm:", nodeToRm)
                    rxToRm.append(nodeToRm)
                else:
                    newDuplicateList.append(duplicateList[poz])
            if len(newDuplicateList) > 1:
                duplicateList = newDuplicateList
            else:
                continue
        rxids = []
        for rxinfo in duplicateList:
            rxes = rxinfo[1]['rxid']
            assert len(rxes) == 1
            rxids.append(rxes[0])
        rxtemp = [_getInfoForRxid(adbData, rxid, 'rxtemp') for rxid in rxids]
        rxABcond = [_getInfoForRxid(adbData, rxid, 'rxacidbase') for rxid in rxids]
        sameTemp = [temp for temp in rxtemp]
        bestMatching = []
        for poz in range(len(rxtemp)):
            thisMatch = [_getCoveredSteps(rxtemp[poz], rxABcond[poz], optCond) for optCond in commonData['allowedCond']['conditionsOcc']]
            thisMatch = sorted([(len(mtch['both']), len(mtch['AB']), len(mtch['temp'])) for mtch in thisMatch], reverse=True)
            bestMatching.append((thisMatch[0], duplicateList[poz][0]))
        bestMatching = sorted(bestMatching, reverse=True)
        if args.debug > 10:
            print("BET", bestMatching)
        for mtch in bestMatching[1:]:
            if args.debug > 10:
                print("remove duplicate:327::", mtch, data['idToFullRx'][mtch[-1]], "FF", seenRxes[rx] )
            rxToRm.append(mtch[-1])
    return rxToRm


def _getMainPathAndOtherRxpoz(seenRx, mainPathRx):
    mainPathRxNodeNum = []
    otherRxNodeNum = []
    assert len(mainPathRx) == 1
    mainPathRx = mainPathRx[0]
    for rx in seenRx:
        rxNodeNum, rxinfo = rx
        rxid = rxinfo['rxid']
        assert len(rxid) == 1
        rxid = rxid[0]
        if rxid == mainPathRx:
            mainPathRxNodeNum.append(rxNodeNum)
        else:
            otherRxNodeNum.append(rxNodeNum)
    return mainPathRxNodeNum, otherRxNodeNum


def _getInfoForRxid(adbData, rxid, whatData):
    if rxid == 1:
        if whatData == 'rxacidbase':
            return ['B', 'WB', 'N', 'WA', 'A']  # no SA and SB
        elif whatData == 'rxtemp':
            return ['VL', 'L', 'RT','H', 'VH']
        else:
            raise NotImplementedError
    lst = adbData['rxes'][rxid][whatData].split('.')
    if whatData == 'rxtemp':
        lst = [t.upper() for t in lst]
    return lst


def _getCoveredSteps(stepTemp, stepAB, optimalCond):
    matchTemp = set()
    for temp in stepTemp:
        if temp in optimalCond['temp']:
            matchTemp.update(set(optimalCond['temp'][temp]))
    matchAB = set()
    for ab in stepAB:
        if ab in optimalCond['AB']:
            matchAB.update(set(optimalCond['AB'][ab]))
    return {'both': matchTemp.intersection(matchAB), 'temp': matchTemp, 'AB': matchAB}


def _getRxPozToRmDueToWildLA(data, commonData, rxInfo, args):
    pozToRm = set()
    if args.rmRxMode in {'mainAwild', 'mainAwildALA'}:
        pozToRmDueToWild = _getRxPozToRmDueToWildMain(data, commonData, rxInfo, args)
        pozToRm.update(set(pozToRmDueToWild))
    if args.rmRxMode in {'mainAwildALA', }:
        pozToRmDueToLA = _getRxPozToRmDueToAandLA(data, commonData, rxInfo, args)
        pozToRm.update(set(pozToRmDueToLA))
    return list(sorted(pozToRm))


def _getRxPozToRmDueToWildMain(data, commonData, rxInfo, args):
    # print("COMMON", commonData['acidBase'])
    # print("DATA", data['acidBase'])
    rxToRm = []
    seenRxes = set()
    for rxpoz in data['idToFullRx']:
        rxsbs, rxprod = data['idToFullRx'][rxpoz]['smiles'].split('>>')
        if '.' in rxprod:
            raise NotImplementedError
        rxsbs = rxsbs.split('.')
        if len(rxsbs) < 2:
            continue
        wildForms = [smi for smi in rxsbs if smi in data['acidBase']['wildForms']]
        if not wildForms:
            continue
        for smi in wildForms:
            mainForm = commonData['acidBase']['wildToMain'][smi]
            if mainForm in rxsbs:
                rxToRm.append(rxpoz)
                # pass
    if rxToRm:
        print("RM", rxToRm, [data['idToFullRx'][rx] for rx in rxToRm])
        # raise
    return rxToRm


def _getRxPozToRmDueToAandLA(data, commonData, rxInfo, args):
    condsAB = [set(cnd[0].split('.')) for cnd in commonData['allowedCond']['full'][0]]
    hasMainLA = [cnd for cnd in condsAB if 'LA' in cnd]
    hasMainAorSA = [cnd for cnd in condsAB if 'A' in cnd or 'SA' in cnd]
    if not(hasMainLA or hasMainAorSA):
        return []
    else:
        print("LALALA", hasMainLA, "C", condsAB)
        #print("DATA", data['acidBase'])
        # raise
    rxToRm = []
    hasMainLAonly = [cnd for cnd in hasMainLA if len(cnd) == 1]
    if not hasMainLAonly:
        return []
    mainPathSet = set(commonData['mainPathInfo']['rxes'])
    hasSideLA = []
    hasSideAorSA = []
    for rxpoz in data['idToFullRx']:
        therx = data['idToFullRx'][rxpoz]
        rxhash = (therx['smiles'], tuple(therx['rxid']))
        rxid = therx['rxid'][0]
        if rxid == 1:
            continue
        if rxhash in mainPathSet:
            continue
        rxcond = rxInfo['rxes'][rxid]['rxacidbase'].split('.')
        if 'LA' in rxcond:
            hasSideLA.append(rxcond)
        if 'A' in rxcond or 'SA' in rxcond:
            hasSideAorSA.append(rxcond)
    if hasSideAorSA and hasMainAorSA:
        print("ISLA", condsAB, "YY", hasSideAorSA)
    if hasSideLA and hasMainLA:
        print("ISLA", condsAB, "XX", hasSideLA, "YY")
    if rxToRm:
        print("RM", rxToRm, [data['idToFullRx'][rx] for rx in rxToRm])
        # raise
    return rxToRm


def _getLateSideToRm(graph, data, commonData, args):
    # assign generation
    smiToGen, toAdd = _assignSmilesGeneration(graph, commonData)
    #print(data['allSmiles'].keys())
    notfound = []
    for smi in toAdd:
        if smi not in data['allSmiles']:
            notfound.append(smi)
            continue
        # print(smi, data['allSmiles'][smi]['generation'])
        smiToGen[smi] = data['allSmiles'][smi]['generation']
    if notfound:
        smiToGen, toAdd = _appendSmilesToSmiToGen(graph, commonData, smiToGen, notfound)


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
        # for rx in data['
        # print("DD", data.keys())
        # for reasbs in
        # raise
    return rxToRm


def assignSmilesToGenInGraph(graph, initsbs):
    smiToGen = dict()
    initNodes = [node for node in graph.nodes() if graph.in_degree(node) == 0]
    for smi in initNodes:
        smiToGen[smi] = 0
    for smi in initsbs:
        if not isinstance(smi, str):
            raise
        print("==SMI", smi)
        smiToGen[smi] = 0
    front = set(initNodes)
    toRm = set()
    itr = 0
    prevFront = None
    while front:
        newToFront = set()
        for node in front:
            for rx in graph.succ[node]:
                if not all(sbs in smiToGen for sbs in graph.pred[rx]):
                    print("RX", rx, [s for s in graph.pred[rx] if s not in smiToGen])
                    continue
                gen = 1 + max(smiToGen[sbs] for sbs in graph.pred[rx])
                for prd in graph.succ[rx]:
                    if graph.out_degree(prd) > 0:
                        newToFront.add(prd)
                    if prd in smiToGen:
                        if smiToGen[prd] > gen:
                            smiToGen[prd] = gen
                    else:
                        smiToGen[prd] = gen
                toRm.add(node)
        front.update(newToFront)
        front -= toRm
        itr += 1
        if itr > 10:
            if not prevFront:
                prevFront = front.copy()
            else:
                # print(f"SMITOGEN {len(smiToGen)}  ndes: {len(graph.nodes())},GR: {graph.nodes()}")
                if front == prevFront and all(smi in smiToGen for smi in front) and len(graph.nodes()) == len(smiToGen):
                    break
            print("ITR", itr, "FR", front, "FRinfo:", {f: graph.succ[f] for f in front})
            for cmd in front:
                for rx in graph.succ[cmd]:
                    print("   ???", cmd, rx, "PRED", {c: smiToGen.get(c, 'N') for c in graph.pred[rx]})
                    for smi in graph.pred[rx]:
                        if smi in smiToGen:
                            print("   !", smi, smiToGen[smi])
                        else:
                            print("    ??", smi, graph.pred[smi])
        if itr > 15:
            print("GR", smiToGen, "FR", front)
            raise
    return smiToGen


def _assignSmilesGeneration(graph, commonData):
    smiToGen = dict()
    for cmd in commonData['mainPathInfo']['smilesToGen']:
        smiToGen[cmd] = commonData['mainPathInfo']['smilesToGen'][cmd]
    toAdd = set([smiOrInt for smiOrInt in graph.nodes() - set(smiToGen.keys()) if isinstance(smiOrInt, str)])
    mainByprods = commonData['mainPathInfo']['byprods']['main']
    alterByprods = commonData['mainPathInfo']['byprods']['alter']
    prevToAdd = None
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


def _appendSmilesToSmiToGen(graph, commonData, smiToGen, toAdd):
    toAdd = set(toAdd)
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


def _getLateSideToRmBasedOnGen(graph, data, commonData, args):
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
            if prd not in data['allSmiles']:
                debug = True
                continue
            gen = data['allSmiles'][prd]['generation']
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
        # for rx in data['
        # print("DD", data.keys())
        # for reasbs in
        # raise
    return rxToRm


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
                raise
                # continue
            revCandProdsSet = set(revCandProds.split('.'))
            if revCandProdsSet == sbsesSet:
                return True
    return False


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


def _addSbsesToDictSmiles(sbses, dictdict, longShortGroupIds, rxnum, rxpoz):
    for sbs in sbses:
        if sbs not in dictdict:
            dictdict[sbs] = dict()
        for longShortGroup in longShortGroupIds:
            if longShortGroup not in dictdict[sbs]:
                dictdict[sbs][longShortGroup] = set()
            dictdict[sbs][longShortGroup].add((rxpoz, rxnum, tuple(sorted(sbses))))
    return dictdict
