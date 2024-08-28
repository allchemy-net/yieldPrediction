import networkx
import consts, helpers


def _isTautomerisationRx(rxinfo):
    usedRx = set(rxinfo['rxid'])
    return usedRx.intersection(consts.tautomerRxid)


def _isWildFormingRx(rxinfo):
    if not ('+' in rxinfo['smiles'] or '-' in rxinfo['smiles']):
        return False
    sbs, prod = rxinfo['smiles'].split('>>')
    sbsmin = sbs.count('-')
    sbsplu = sbs.count('+')
    prodmin = prod.count('-')
    prodplu = prod.count('+')
    if not (prodmin or prodplu):
        return False
    if sbsmin == prodmin and sbsplu == prodplu:
        return False
    if not ((prodplu and prodplu == sbsplu + 1) or (prodmin and prodmin == sbsmin + 1)):
        return False
    prod = prod.split('.')
    mainForm, _, _, formTypeAB, formStability = helpers.getMainForm(prod[0])
    if mainForm:
        sbs = set(sbs.split('.'))
        if mainForm in sbs:
            return True
    return False


def removeRxFromGraph(graph, rxPozToRm, initSubstrates, args):
    # print("ALLEDGESbefore:", graph.edges())
    removedRxes = set()
    problems = []
    removedNoProblem = []
    itr = 0
    if args.debug > 5:
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
        if args.debug > 5:
            print(f"RX TO RM ITR {itr} :: {rxPozToRm}")
    if problems:
        print("non existing nodes:", problems, "REMOVED:", removedNoProblem, "RXP", rxPozToRm)
    # print("ALLEDGESafter:", graph.edges())
    if args.debug:
        print("REMOVED RX from graph", removedRxes, "PROBLEMS?", problems, "RMED", removedNoProblem)
    return graph, removedRxes


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
        # print("RMed", rxt)
        # rxsma = '.'.join(rxt[0]) + '>>' + rxt[1]
        # rxk = (rxsma, rxt[2])
        # for i in data['rxToFullDict']:
        #    print(type(i), i, rxk)
        # del data['rxToFullDict'][rxk]
    # print("GRAP", [type(el) for el in data['graph']])
    # print("TOF", data['rxToFullDict'].keys())
    # print(f"toid: {data['rxtupleToId'].keys()}")
    # DATA  'rxToFullDict', 'graph', 'rxtupleToId'
    return data


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


# #####


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
                rxhash = (thisRx['fullRxSmiles'], tuple(thisRx['rxid']))
                if rxhash in mainPathSet:
                    continue
                sbs, _ = thisRx['smiles'].split('>>')
                sbs = set(sbs.split('.'))
                if allTargets.intersection(sbs):
                    thisRx['isSideFromTarget'] = True
    return dataDict
