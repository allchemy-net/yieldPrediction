# #################0#####################################
# rxEnablePoint = _getRxEnablePoint(allowedCond, mainPathInfo, args)
# rxEnableCmd = _getRxEnableCmd(mainPathInfo, args)
# dataDict[name][gen]['data']['rxAvail'] = _addRxAvailInfo(dataDict[name][gen]['data'], dataDict[name][0]['common'], args)


def addRxAvailInfo(data, common, args):
    usageRange = _getUsageRange(data, common, args)
    activatingCmd = _getActivatingCmd(data, common, args)
    retDict = {'gen': usageRange, 'cmd': activatingCmd}
    if common['isOnepot']:
        if args.debug > 10:
            print(f"onePot option: {common['onePot']}")
            print(f"RETDD {retDict} CMD {common.keys()}")
        onepotActToRx = _getActivatingCmdOnepot(data, common, args)
        if args.debug:
            print("ONEpotInfo", onepot)
        for activator in onepotActToRx:
            if activator in activatingCmd:
                print("ACTIVATOR PROBLEM", activator, activatingCmd)
                # ACTIVATOR PROBLEM C=CCC1CCCCC12OC(=O)C2CC=C {'C=CCC1CCCCC12OC(=O)C2CC=C': {3}}
                raise NotImplementedError
            activatingCmd[activator] = onepotActToRx[activator]
        retDict['cmd'] = activatingCmd
    #if args.debug > 10:
    if True:
        print("RXavailInfo", retDict)
    return retDict


def _getActivatingCmdOnepot(data, common, args):
    if 'onePot' not in common:
        raise NotImplementedError
    # cmdToRx = dict()  # cmd which enable given reactions
    sbsToActivator = dict()
    for cmd in common['onePot']:
        for rxdict in common['onePot'][cmd]:
            sbs = rxdict['subs']
            if sbs not in sbsToActivator:
                sbsToActivator[sbs] = set()
            sbsToActivator[sbs].add(cmd)
    activatorToRx = dict()
    sbsToActivatorSet = set(sbsToActivator)
    print("SBS TO ACTIVATOR", sbsToActivator, "\nCMDONEPOT", common['onePot'])
    print("GENS", {smi: data['allSmiles'][smi]['generation'] for smi in data['allSmiles']})
    for rxnum in data['idToFullRx']:
        rxSbses = set(data['idToFullRx'][rxnum]['fullRxSmiles'].split('>>')[0].split('.'))
        usedSbs = sbsToActivatorSet.intersection(rxSbses)
        if not usedSbs:
            continue
        if len(usedSbs) != 1:
            print(f"some reaction has more than one activator {usedSbs} RXXX: {data['idToFullRx'][rxnum]}")
            # raise NotImplementedError
        for sbs in usedSbs:
            #if data['allSmiles'][sbs]['generation'] == 0:
            #    continue
            activator = sbsToActivator[sbs]
            if len(activator) > 1:
                activator = [act for act in activator if act not in rxSbses and data['allSmiles'][act]['generation'] >= 1]
            if not activator:
                continue
            if len(activator) != 1:
                print(f"sbs: {sbs} infogen:  {data['allSmiles'][sbs]['generation']}")
                print(f"MULTIPLE ACTIVATORS {sbs} {activator} in RX: {data['idToFullRx'][rxnum]}")
                raise NotImplementedError
            activator = tuple(activator)[0]
            if activator not in activatorToRx:
                activatorToRx[activator] = set()
            activatorToRx[activator].add(rxnum)
            if data['allSmiles'][sbs]['generation'] == 0:
                print(f"notIGNORED:: {activator} {rxnum} :ALL: {activatorToRx[activator]}")
    print("FULLACTIVATOR", activatorToRx)
    return activatorToRx


def _getActivatingCmd(data, common, args):
    cmdToBlockedRx = dict()
    if 'rxEnableCmd' not in common:
        return dict()
    catToCmd = dict()
    if args.debug > 100:
        print("RXENABLE", common['rxEnableCmd'])
    for cmd in common['rxEnableCmd']:
        for rxdict in common['rxEnableCmd'][cmd]:
            for cat in rxdict['cat']:
                if cat not in catToCmd:
                    catToCmd[cat] = {cmd, }
                else:
                    catToCmd[cat].add(cmd)
    if args.debug > 10:
        for cmd in common['rxEnableCmd']:
            print("rxEnableCmd:: cmd::", cmd, common['rxEnableCmd'][cmd])
        print("CATTOCMD", catToCmd)
    catSet = set(catToCmd.keys())
    cmdToRx = dict()
    for rxnum in data['idToFullRx']:
        rxSbses = set(data['idToFullRx'][rxnum]['fullRxSmiles'].split('>>')[0].split('.'))
        usedCat = rxSbses.intersection(catSet)
        if args.debug > 120:
            print(f"RX_act_cmd {rxnum} :: {rxSbses} {catSet}")
        for cat in usedCat:
            for cmd in catToCmd[cat]:
                if cmd not in cmdToRx:
                    cmdToRx[cmd] = {rxnum, }
                else:
                    cmdToRx[cmd].add(rxnum)
    if args.debug > 10:
        print("CMD TO RX", cmdToRx)
        for cmd in cmdToRx:
            for rx in cmdToRx[cmd]:
                print(f"CMDtoRX {cmd}  {rx} {data['idToFullRx'][rx]['rxid']}  {data['idToFullRx'][rx]['fullRxSmiles']} ")
    #raise
    return cmdToRx


def _getUsageRange(data, common, args):
    """ get usage range for for each graph; arguments:
        - data - dataDict for one particular graph for given level of additional calculation
        - common - dataDict with common data for path (same for 0th 1st and 2nd level of additional calculation)
    """
    usageRange = dict()
    catUsage = dict()
    detectedCats = set()
    if 'rxEnablePoint' not in common:
        # usageRange[rxnum] = (None, None)
        return usageRange
    if 'rxEnablePoints' in common:
        for cmdlist in common['rxEnablePoints']:
            print("USAGE", cmdlist)
    for rxnum in data['idToFullRx']:
        rxSbses = set(data['idToFullRx'][rxnum]['fullRxSmiles'].split('>>')[0].split('.'))
        begin = None  # if int 1 or more
        end = None
        usedCat = rxSbses.intersection(set(common['rxEnablePoint'].keys()))
        for cat in usedCat:
            minGen = min(common['rxEnablePoint'][cat])
            if cat not in catUsage:
                catUsage[cat] = minGen
            else:
                if catUsage[cat] > minGen:
                    catUsage[cat] = minGen
            if not begin or begin > minGen:
                begin = minGen
        if not begin:
            continue
        usageRange[rxnum] = {'min': begin, 'max': end}
    # check reactions which use catalysts as normal reagents
    # such reaction also need to be blocked
    catSet = set(catUsage.keys())
    for rxnum in data['idToFullRx']:
        rxSbses = set(data['idToFullRx'][rxnum]['fullRxSmiles'].split('>>')[0].split('.'))
        usedCat = rxSbses.intersection(catSet)
        if not usedCat:
            continue
        minGen = min([catUsage[cat] for cat in usedCat])
        if rxnum in usageRange:
            if minGen < usageRange[rxnum]['min']:
                usageRange[rxnum]['min'] = minGen
            continue
        # no in usageRange
        usageRange[rxnum] = {'min': minGen, 'max': None}
    if args.debug > 40:
        print("USAGE", usageRange, "ENABLEpoint::", common.get('rxEnablePoint', 'NON'))
        print("USAGE", [(rxnum, data['idToFullRx'][rxnum]['rxid'], data['idToFullRx'][rxnum]['smiles']) for rxnum in usageRange])
    return usageRange


def getRxEnableCmdOnepot(mainPathInfo, args):
    gen0 = [smi for smi in  mainPathInfo['smilesToGen'] if mainPathInfo['smilesToGen'][smi] == 0]
    mainGraph = mainPathInfo['graph']
    # gen0nodes = [n for n in mainGraph.nodes() if n in set(gen0)]
    depCmdCat = dict()
    for smi in gen0:
        for rx in mainGraph.succ[smi]:
            pgens = {cmd: mainPathInfo['smilesToGen'][cmd] for cmd in mainGraph.pred[rx] if mainPathInfo['smilesToGen'][cmd] > 0}
            if not pgens:
                continue
            if len(pgens) != 1:
                raise NotImplementedError
            sbs, _ = pgens.popitem()
            if sbs not in depCmdCat:
                depCmdCat[sbs] = []
            depCmdCat[sbs].append({'rx': rx, 'subs': smi})
    return depCmdCat


def getRxEnableCmd(mainPathInfo, args):
    if not mainPathInfo.get('hasCat', False):
        return None
    if args.catAddTime != 'peak':
        raise NotImplementedError
    depCmdCat = dict()
    mainPathProds = {smarx[0].split('>>')[1]: smarx for smarx in mainPathInfo['rxes']}
    for rxWithCat in mainPathInfo['hasCat']:
        allrea = rxWithCat[0].split('>>')[0].split('.')
        reactants = [rea for rea in allrea if rea in mainPathProds]
        # reactant from reviousgeneration which need to accomulate before the rx starts
        if not reactants:
            # reactants = [rea for rea in allrea if rea in mainPathProds or rea in mainPathInfo['initSubstrates']]
            continue
        if len(reactants) == 2 and '[Pd]' in reactants:
            reactants.remove('[Pd]')
        if len(reactants) != 1:
            print("ALLREA", allrea, "MP", mainPathInfo['initSubstrates'])
            print("RXWCAT", rxWithCat, "EACTANTS", reactants, "MAINPATH", mainPathProds.keys() )
            print("MA", mainPathInfo['hasCat'])
            raise NotImplementedError
        if reactants[0] not in depCmdCat:
            depCmdCat[reactants[0]] = []
        depCmdCat[reactants[0]].append({'rx': rxWithCat, 'cat': mainPathInfo['hasCat'][rxWithCat]})
    return depCmdCat


def getRxEnablePoint(allowedCond, mainPathInfo, args):
    retdict = dict()
    if not mainPathInfo.get('hasCat', False):
        return None
    if args.catAddTime in {'jitEarly', 'jitLate'}:
        res = _getRxEnableGen(allowedCond, mainPathInfo, args)
        retdict['gen'] = res
    elif args.catAddTime == 'peak':
        res = _getRxEnableCmd(allowedCond, mainPathInfo, args)
        retdict['cmd'] = res
    else:
        raise NotImplementedError
    return res


def _getRxEnableGen(allowedCond, mainPathInfo, args):
    lateCats = dict()
    lenpath = len(allowedCond['minimalVariants'][0])
    for rxWithCat in mainPathInfo['hasCat']:
        if args.catAddTime == 'jitEarly':
            gen = allowedCond['rxToGenSbs'][rxWithCat]
        elif args.catAddTime == 'jitLate':
            gen = lenpath - allowedCond['rxToGenProd'][rxWithCat]
            # print("RXWC", allowedCond['rxToGenProd'], "===",  rxWithCat, "GEN", gen)
        else:
            print("MODE", args.catAddTime)
            raise NotImplementedError
        if gen == 1:
            continue
        cats = mainPathInfo['hasCat'][rxWithCat]
        for cat in cats:
            if cat not in lateCats:
                lateCats[cat] = []
            # print("USAGE", cat)
            lateCats[cat].append(gen)
    if not lateCats:
        return None
    # print("ALLOWEDCON", allowedCond['rxToGenSbs'], "==P==", allowedCond['rxToGenProd'])
    # print("MAIN", mainPathInfo['hasCat'], "LATE", lateCats, "ENN", lenpath)
    return lateCats


def _getRxEnableCmd(allowedCond, mainPathInfo, args):
    pass