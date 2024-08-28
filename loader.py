import pickle, json, itertools
import networkx
from networkx.drawing.nx_pydot import write_dot
from rdkit import Chem
import helpers, consts, loaderCat, loaderGraph, loaderHelper


def loadOnepotCmds(fn):
    onepots = set()
    if not fn:
        return onepots
    fh = open(fn)
    for line in fh:
        line = line.strip()
        if line.startswith('#'):
            continue
        onepots.add(line)
    fh.close()
    return onepots


# ######################################################
# load pickle
# ##


def loadPickle(fnames, rxInfo, args):
    dataDict = dict()
    longShortCut = None
    if args.disableLongCut > 0:
        longShortCut = parseShortCutFile(args.longCutFile)
    reactiveSmarts = None
    shortcutToBanSet = _getShortcutToBan(longShortCut)
    if args.reactiveSmarts:
        reactiveSmarts = loadReactiveSmarts(args.reactiveSmarts, args)
    baseUsage = loadBaseUsageInfo(args.baseUsage)
    onePotSet = loadOnepotCmds(args.onepotfile)
    for fn in fnames:
        try:
            basename, rest = fn.split('.calcedGen_')
            basename = basename.split(r'/')[-1]
        except:
            print("PROBLEM WITH FN", fn)
            raise
        gen, _ = rest.split('.')
        gen = int(gen)
        if basename not in dataDict:
            dataDict[basename] = dict()
        fh = open(fn, 'rb')
        data = pickle.load(fh)
        fh.close()
        dataDict[basename][gen] = {'data': data, }
        inplaceAddFullRxSmi(dataDict[basename][gen])

    for name in dataDict:
        if args.debug:
            print("NAME", name)
        wildToMain = dict()  # smiles : smiles
        mainToWild = dict()  # smiles : {smiles1, smiles2}
        wildFormType = dict()
        noMainForm = set()

        fixed = _detectFixed(dataDict[name][0]['data']['rxToFullDict'])
        if args.debug > 100:
            print(f"FIXED:67 {fixed}")
        additives = _detectAdditives(dataDict[name])
        if args.debug > 100:
            print(f"ADDIT:72: {additives}")
        mainPathInfo = getMainPathInfo(dataDict[name][0], rxInfo, reactiveSmarts, args)
        allowedCond, mainPathInfo = detectAllowedCond(dataDict[name][0], rxInfo, mainPathInfo, args)
        if args.debug > 100:
            print(f"ALLOWEDCOND::76: {allowedCond}")
        theCond = getWorkingConditions(allowedCond)
        allFixed = set()
        for rx in fixed:
            _ = [allFixed.add(cmd) for cmd in fixed[rx]]
        allInit = mainPathInfo['initSubstrates'].union(allFixed)
        isOnepot = name.split('_')[0] in onePotSet
        print(f"ONEPOT? {name} ??? {isOnepot} ALLL: {onePotSet}")
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
                tmptuple = _addWildMainToInfoForDatadictNameGen(dataDict[name][gen], wildToMain, mainToWild, noMainForm, wildFormType)
                dataDict[name][gen], wildToMain, mainToWild, noMainForm, wildFormType = tmptuple
                # uptoupdate
                dataDict[name][0]['common']['acidBase'] = {'noMainForm': noMainForm, 'wildToMain': wildToMain, 'mainToWild': mainToWild, 'wildFormType': wildFormType}
                if args.scaleSideInOldCondition:
                    dataDict[name][gen]['data']['idToFullRx'] = _addScalingForSideInOldConditions(dataDict[name][gen]['data'],
                                                                                                  args, dataDict[name][0]['common'])
                dataDict[name][gen]['data']['idToFullRx'] = _addShortcutInfo(dataDict[name][gen]['data']['idToFullRx'], shortcutToBanSet, args)
                dataDict[name][gen]['data']['idToFullRx'] = _addDeprotonationInfo(dataDict[name][gen]['data']['idToFullRx'], args)
                dataDict[name][gen]['data']['rxAvail'] = loaderCat.addRxAvailInfo(dataDict[name][gen]['data'], dataDict[name][0]['common'], args)
                # below also ['data']['idToFullRx'] may be changed
                networkxGraph, dataDict[name][gen]['data'] = loaderGraph._buildGraph(dataDict[name][gen]['data']['graph'], additives,
                                                                                     dataDict[name][gen]['data'], rxInfo, longShortCut,
                                                                                     dataDict[name][0]['common'], args)
                # dataDict[name][gen]['data'] is updated dataDict[name][gen]['data']['idToFullRx']
                rxPozToRm = set()
                if args.includepka != 'N':
                    rxPozToRm.update(deprotonation.findCompetingDeprotToRm(networkxGraph, dataDict[name][gen]['data']['idToFullRx'], name, baseUsage, args))
                    if args.debug:
                        print("REMOVE DUE TO PKA:", name, gen, rxPozToRm)
                if rxPozToRm:
                    networkxGraph, removedRxPoz = loaderHelper.removeRxFromGraph(networkxGraph, rxPozToRm, allInit, args)
                    if args.debug > 8:
                        tmpdata = dataDict[name][gen]['data']['idToFullRx']
                        inforxtorm = [(tmpdata[rxpoz]['smiles'], tmpdata[rxpoz]['rxid']) for rxpoz in rxPozToRm]
                        if args.debug:
                            print(f"REMOVED_RX_FROM_GR:91 {rxPozToRm} RM: {inforxtorm}")
                    dataDict[name][gen]['data'] = loaderHelper.removeRxFromData(dataDict[name][gen]['data'], removedRxPoz)
                printRx = False
                if printRx:
                    for rxnum in dataDict[name][gen]['data']['idToFullRx']:
                        dct = dataDict[name][gen]['data']['idToFullRx'][rxnum]
                        print("F", dct['rxid'], dct['smiles'])
                dataDict[name][gen]['graph'] = networkxGraph
                if args.exportGraph:
                    graphName = f"{name}_gen{gen}"
                    fn = f'exportedGraph_{args.exportGraph}_{graphName}'
                    fh = open(fn, 'w')
                    json.dump(tuple(networkxGraph.edges()), fh)
                    fh.close()
                    fname = f'exportedGraph_{args.exportGraph}_{graphName}.dot'
                    write_dot(networkxGraph, fname)
            except:
                print("PROBLEM WITH", name, gen)
                raise
    return dataDict


def _addWildMainToInfoForDatadictNameGen(dataDictNameGen, wildToMain, mainToWild, noMainForm, wildFormType):
    # dataDict[name][gen]['data'] = _addTautomerInfo(dataDict[name][gen]['data'], args)
    thisWilds = set()
    for idrx in dataDictNameGen['data']['idToFullRx']:
        if dataDictNameGen['data']['idToFullRx'][idrx].get('origin', ''):
            # ignore reaction reversed by "hand"
            continue
        sbses, prods = dataDictNameGen['data']['idToFullRx'][idrx]['smiles'].split('>>')
        sbses = sbses.split('.')
        wildToMain, mainToWild, noMainForm, wildFormType = _addWildMainToInfo(sbses, wildToMain, mainToWild, noMainForm, wildFormType)
        prods = prods.split('.')
        if len(prods) != 1:
            print("XXX", prods, "RX", dataDictNameGen['data']['idToFullRx'][idrx])
            raise NotImplementedError
        wildToMain, mainToWild, noMainForm, wildFormType = _addWildMainToInfo(prods, wildToMain, mainToWild, noMainForm, wildFormType)
        for cmdList in (sbses, prods):
            for smi in cmdList:
                if smi in wildToMain:
                    thisWilds.add(smi)
    dataDictNameGen['data']['acidBase'] = {'wildForms': thisWilds}  # , 'xOfMainForm': thisMains}
    # objects in dict should link to the same wildToMain mainToWild and noMainForm
    # dataDictName[0]['common']['acidBase'] = {'noMainForm': noMainForm, 'wildToMain': wildToMain, 'mainToWild': mainToWild, 'wildFormType': wildFormType}
    return dataDictNameGen, wildToMain, mainToWild, noMainForm, wildFormType


def _addWildMainToInfo(smilesList, wildToMain, mainToWild, noMainForm, wildFormType):
    for smi in smilesList:
        if smi in wildToMain or smi in mainToWild or smi in noMainForm:
            continue
        mainForm, _, _, formTypeAB, formTypeStability = smilesHelper.getMainForm(smi)
        # if mainForm and mainForm != smi and formTypeAB not in {'A', 'B'}:
        #    print(f"MAIN: {mainForm} in: {smi} type: {formTypeAB} stab: {formTypeStability}")
        #    raise
        if mainForm and mainForm != smi:
            wildToMain[smi] = mainForm
            if mainForm not in mainToWild:
                mainToWild[mainForm] = set()
            mainToWild[mainForm].add(smi)
            wildFormType[smi] = formTypeAB
        else:
            noMainForm.add(smi)
    return wildToMain, mainToWild, noMainForm, wildFormType


def _addDeprotonationInfo(fullrxDict, args):
    for rxnum in fullrxDict:
        isDeprot = False
        if set(fullrxDict[rxnum]['rxid']).intersection(consts.shortcutToPkaCheck):
            deprotonInfo = deprotonation.getComplexDeprotInfo(fullrxDict[rxnum])
            isDeprot = ['Complex', deprotonInfo]  # deprotonated form will be added
            if args.debug > 50:
                print("complex deprotonation", fullrxDict[rxnum])
        elif _isDeprotonationRx(fullrxDict[rxnum]):
            if args.debug > 50:
                print("simple deprotonation", fullrxDict[rxnum])
            deprotonInfo = deprotonation.getSimpleDeprotInfo(fullrxDict[rxnum])
            isDeprot = ['Simple', deprotonInfo]
        fullrxDict[rxnum]['isDeprot'] = isDeprot
    return fullrxDict


def _isDeprotonationRx(rxinfo):
    sbs, prod = rxinfo['smiles'].split('>>')
    panion = prod.count('-')
    if not panion:
        return False
    sanion = sbs.count('-')
    if sanion - panion >= 0:
        return False
    return loaderHelper._isWildFormingRx(rxinfo)


def _addShortcutInfo(fullrxDict, shortcutToBanSet, args):
    for rxnum in fullrxDict:
        isShortCut = False
        if set(fullrxDict[rxnum]['rxid']).intersection(shortcutToBanSet):
            if args.debug > 50:
                print("mark rx as shortcut", fullrxDict[rxnum])
            isShortCut = True
        fullrxDict[rxnum]['isShortcut'] = isShortCut
    return fullrxDict


def _getShortcutToBan(longShortCut):
    fullset = set(longShortCut['shortRxesDict'].keys())
    fullset.update(consts.shortcutToPkaCheck)
    return fullset


def parseShortCutFile(fname):
    fh = open(fname)
    info = []
    for lid, line in enumerate(fh):
        if lid == 0:
            continue
        elems = line.split('\t')
        shortcutName, shortcutRxes, shortcutCond = elems[0:3]
        longcutName, longcutRxes, longcutCond = elems[3:6]
        shortcutRxes = resolveRxRange(shortcutRxes)
        longcutRxes = resolveRxRange(longcutRxes)
        info.append({'shortRxes': shortcutRxes, 'shortName': shortcutName, 'shortCond': shortcutCond,
                     'longRxes': longcutRxes,  'longName': longcutName, 'longCond': longcutCond})
        if len(elems) > 6:  # has Long2
            longcut2Name, longcut2Rxes, longcut2Cond = elems[6:9]
            longcut2Rxes = resolveRxRange(longcut2Rxes)
            info.append({'shortRxes': shortcutRxes, 'shortName': shortcutName, 'shortCond': shortcutCond,
                         'longRxes': longcut2Rxes,  'longName': longcut2Name, 'longCond': longcut2Cond})
            info.append({'shortRxes': longcut2Rxes, 'shortName': longcut2Name, 'shortCond': longcut2Cond,
                         'longRxes': longcutRxes,  'longName': longcutName, 'longCond': longcutCond})
    fh.close()
    longDict = dict()
    shortcutDict = dict()
    for idx, infodict in enumerate(info):
        for rxid in infodict['longRxes']:
            if rxid not in longDict:
                longDict[rxid] = set()
            longDict[rxid].add(idx)
        for rxid in infodict['shortRxes']:
            if rxid not in shortcutDict:
                shortcutDict[rxid] = set()
            shortcutDict[rxid].add(idx)
    retDict = {'fullInfo': info, 'longRxesDict': longDict, 'shortRxesDict': shortcutDict}
    # print("longCUTfull", retDict)
    return retDict


# #####################################################
# load


def loadMutexRxFile(fname):
    fh = open(fname)
    dictfs = dict()
    for line in fh:
        _, rxids1, _, rxids2 = line.split('\t')
        rxids1 = resolveRxRange(rxids1.strip()) #.split(',')
        rxids2 = resolveRxRange(rxids2.strip()) #.split(',')
        dictfs[frozenset(rxids1)] = frozenset(rxids2)
    fh.close()
    # print("MUTEXfromFile:", dictfs)
    return dictfs


def resolveRxRange(rxrange):
    rxes = set()
    for rx in rxrange.split(','):
        if '-' not in rx:
            rxes.add(rx.strip())
            continue
        minrange, maxrange = rx.strip().split('-')
        suffix = ''
        if minrange[-1] == maxrange[-1] and not minrange[-1].isdigit():
            suffix = minrange[-1]
            minrange = minrange[:-1]
            maxrange = maxrange[:-1]
        prefix = ''.join([char for char in minrange if not char.isdigit()])
        minValue = int(minrange.replace(prefix, ''))
        maxValue = int(maxrange.replace(prefix, ''))
        for num in range(minValue, maxValue + 1):
            rxes.add(f'{prefix}{num}{suffix}')
    return rxes


# ###################################
def loadReactiveSmarts(fn, args):
    retDict = {'obj': [], 'smarts': []}
    fh = open(fn)
    added = set()
    for line in fh:
        line = line.split('.')
        for sma in line:
            sma = sma.strip()
            if not sma:
                continue
            if sma in added:
                if args.debug:
                    print("INFO: duplicated smarts", sma)
                continue
            mol = Chem.MolFromSmarts(sma)
            if not mol:
                print("incorrect smarts", sma)
                raise consts.WrongData
            retDict['obj'].append(mol)
            retDict['smarts'].append(sma)
            added.add(sma)
    fh.close()
    return retDict


def inplaceAddFullRxSmi(data):
    for rx in data['data']['rxToFullDict']:
        rxinfo = data['data']['rxToFullDict'][rx]
        mapped = rxinfo['additionalInfo']['mappedRx']
        sbs, prod = mapped.split('>>')
        sbs = [helpers.demapSmiles(s) for s in sbs.split('.')]
        prod = [helpers.demapSmiles(s) for s in prod.split('.')]
        demaped = '.'.join(sbs) + '>>' + '.'.join(prod)
        data['data']['rxToFullDict'][rx]['fullRxSmiles'] = demaped
    for rxtype in ['med', 'mech']:
        for poz, rxlist in enumerate(data['data']['competing'][rxtype]):
            for rxinfo in rxlist:
                mapped = rxinfo['mappedRx']
                sbs, prod = mapped.split('>>')
                sbs = [helpers.demapSmiles(s) for s in sbs.split('.')]
                prod = [helpers.demapSmiles(s) for s in prod.split('.')]
                demaped = '.'.join(sbs) + '>>' + '.'.join(prod)
                rxinfo['fullRxSmiles'] = demaped


# #####################################


def parseAdb(fn, args):
    fh = open(fn)
    rxidToData = {'rxes': dict(), 'fastRxSet': set(), 'slowRxSet': set(), 'verySlowRxSet': set(), 'reversibleRx': set(),
                  'redoxred': set(), 'redoxox': set(), 'byprods': {'main': dict(), 'alter': dict()}, }
    fw = open('reakcje_tautomeryzacje.tsv', 'w')
    for line in fh:
        if not line.startswith('RX\t'):
            continue
        if 'mech' not in line:
            if 'none' in line:
                pass
            else:
                continue
        rxinfo = json.loads(line[3:])
        rxid = rxinfo['rxidOryg']
        if rxinfo["byproducts"]:
            for smi in rxinfo['byproducts'].split('.'):
                smi = smi.strip()
                if not smi or smi == '-':
                    continue
                if smi not in rxidToData['byprods']['main']:
                    rxidToData['byprods']['main'][smi] = set()
                rxidToData['byprods']['main'][smi].add(rxid)
        if rxinfo['alterByproducts']:
            for smi in rxinfo['byproducts'].split('.'):
                smi = smi.strip()
                if not smi or smi == '-':
                    continue
                if smi not in rxidToData['byprods']['alter']:
                    rxidToData['byprods']['alter'][smi] = set()
                rxidToData['byprods']['alter'][smi].add(rxid)
        rxidToData['rxes'][rxid] = rxinfo
        if 'automer' in rxinfo['rxname'] and 'romati' not in rxinfo['rxname'] and 'isomer' not in rxinfo['rxname']:
            print(line[:-1], file=fw)
        if rxinfo["rxspeed"].strip().upper() == 'F':
            rxidToData['fastRxSet'].add(rxid)
        elif rxinfo["rxspeed"].strip().upper() == 'S':
            rxidToData['slowRxSet'].add(rxid)
        elif rxinfo["rxspeed"].strip().upper() == 'VS':
            rxidToData['verySlowRxSet'].add(rxid)
        #
        if rxinfo['irreversible'] == 'N':
            rxidToData['reversibleRx'].add(rxid)
        rxredox = rxinfo["rxoxred"].upper().split('.')
        isRed = 'RED' in rxredox
        isOx = 'OX' in rxredox
        if isRed and isOx:
            # print("RX", rxinfo)
            # raise NotImplementedError
            pass
        elif isRed:
            rxidToData['redoxred'].add(rxid)
        elif isOx:
            rxidToData['redoxox'].add(rxid)
    fh.close()
    return rxidToData


def _detectFixed(rxToFullDict):
    # dataDict[n][0]['data']['rxToFullDict']
    fixedDict = dict()
    for rx in rxToFullDict:
        rxinfo = rxToFullDict[rx]
        sbs, prod = rxinfo['smiles'].split('>>')
        fullsbs, fullprods = rxinfo['fullRxSmiles'].split('>>')
        sbses = sbs.split('.')
        fullsbses = fullsbs.split('.')
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
        try:
            additives = data[maxgen]['data']['rxToFullDict'][rx]['byproducts']['additives']
            for additeList in additives:
                for addit in additeList:
                    if addit and addit != '-':
                        canon = Chem.CanonSmiles(addit)
                        addits.add(canon)
        except:
            pass
        if addits:
            sbs, prod = rx[0].split('>>')
            rxtuple = (tuple(sorted(sbs.split('.'))), prod, rx[1])
            rxidgraph = data[maxgen]['data']['rxtupleToId'][rxtuple]
            retDict[rxidgraph] = addits
    return retDict


# ###############################################
# condition AB and temperature


def getMainPathInfo(data, rxInfo, reactiveSmarts, args):
    mainPathInfo = {'isWaterSbs': False, 'isWaterCondNeeded': False, 'reactiveSubstrates': dict(), 'rxes': [], 'byprods': rxInfo['byprods'],
                    'redox': set(), 'hasCat': dict(), 'initSubstrates': set(), 'smilesToGen': dict(), 'graph':None}
    checkedCmd = set()
    #dg = data['graph']
    allSbses = set()
    allProds = set()
    graph = networkx.DiGraph()
    for rxsmaRxid in data['data']['rxToFullDict']:
        rxsma, rxid = rxsmaRxid
        fullsbs, fullprods = data['data']['rxToFullDict'][rxsmaRxid]['fullRxSmiles'].split('>>')
        for sbs in fullsbs.split('.'):
            graph.add_edge(sbs, rxsmaRxid)
        for prd in fullprods.split('.'):
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
        rxids = tuple(rxid)
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
        if True:  # args.catAddTime != 'init':
            # limit consideration to main path reaction only
            # print("ARGS", args.catDef, fullsbs)
            if args.catDef == 'Pd':
                if '[Pd]' in fullsbs:
                    #print("PDTEST", rxsmaRxid, "FF:", data['data']['rxToFullDict'][rxsmaRxid])
                    mainPathInfo['hasCat'][rxsmaRxid] = '[Pd]'
            elif args.catDef == 'fixPd':
                fullrx = data['data']['rxToFullDict'][rxsmaRxid]['fullRxSmiles']
                fullsbs = set(fullrx.split('>>')[0].split('.'))
                cat = fullsbs - sbsesSet
                if '[Pd]' in fullsbs:
                    cat.add('[Pd]')
                if cat:
                    mainPathInfo['hasCat'][rxsmaRxid] = tuple(cat)
                    print(f"HAS CAT {rxsmaRxid}")
            else:
                raise NotImplementedError
        if reactiveSmarts:
            for sbs in sbses:
                if sbs in checkedCmd:
                    continue
                mol = Chem.MolFromSmiles(sbs)
                checkedCmd.add(sbs)
                for poz, patt in enumerate(reactiveSmarts['obj']):
                    if mol.HasSubstructMatch(patt):
                        if sbs not in mainPathInfo['reactiveSubstrates']:
                            mainPathInfo['reactiveSubstrates'][sbs] = []
                        mainPathInfo['reactiveSubstrates'][sbs].append(reactiveSmarts['smarts'][poz])
                        break
    initSbs = allSbses - allProds
    mainPathInfo['initSubstrates'] = initSbs
    mainPathInfo['graph'] = graph
    cats = set()
    for smis in mainPathInfo['hasCat'].values():
        if isinstance(smis, (tuple, set, list)):
            # catDef fixPd
            _ = [cats.add(s) for s in smis]
        else:
            # catDef Pd
            cats.add(smis)
    sbsAndCat = initSbs.union(cats)
    mainPathInfo['smilesToGen'] = loaderGraph.assignSmilesToGenInGraph(graph, sbsAndCat)
    return mainPathInfo


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


def detectAllowedCond(data, rxinfo, mainPathInfo, args):
    topos = {'full': [], 'minimalVariants': [], 'stepsVariants': [], 'catUsage': [], 'conditionsOcc': []}
    topoOrders, rxgraph = _getAllowedTopoOrder(data)
    doPrint = False
    if not topoOrders or (len(topoOrders) == 1 and not topoOrders[0]) :
        print("MAINPATINF", mainPathInfo)
        rxid = mainPathInfo['rxes'][0][1][0]
        print("RXID", rxid)
        print("COND",  rxinfo['rxes'][rxid])
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
        uniqABint = [[consts.condABtoInt[cnd] for cnd in cond] for cond in uniqAB]
        bestABint = _getBestCondCombinations(uniqABint)
        # empty frozenset if fo reaction in which only LA is spec. in conditionsAB
        bestAB = [consts.condIntToAB.get(cint, frozenset()) if cint is not None else set() for cint in bestABint]
    elif mode == 'T':
        uniqABint = [[consts.condTempToInt[cnd] for cnd in cond] for cond in uniqAB]
        bestABint = _getBestCondCombinations(uniqABint)
        bestAB = [consts.condIntToTemp[cint] for cint in bestABint]
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
    mintopo =  list((frozenset(ab), frozenset(temp)) for ab, temp in zip(minab, mintemp))
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


def getHeatOfFormation(args):
    data = dict()
    for smiHeat in open(args.heatFormation):
        smiHeat = smiHeat.strip()
        if not smiHeat or smiHeat[0] == '#':
            continue
        try:
            smi, heat = smiHeat.split()
            data[smi] = float(heat)
        except:
            print("LINE DE PROBLEM", smiHeat)
            raise
    print(f"LOADED HoF for {len(data)} smiles from {args.heatFormation}")
    return data


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


def _mapIdToFullRx(dataDict, args):
    retDict = dict()
    added = []
    for rxsma, rxid in dataDict['rxToFullDict']:
        sbs, prod = rxsma.split('>>')
        sbs = sorted(sbs.split('.'))
        rxtuple = (tuple(sbs), prod, rxid)
        if rxtuple not in dataDict['rxtupleToId']:
            print("NF", rxtuple)
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
                rxtuple1, rxtuple2 =  rxtuple, None
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


# tautomer info
def addTautomerInfoToDataset(dataDict, args):
    for name in dataDict:
        for gen in dataDict[name]:
            dataDict[name][gen]['data'] = _addTautomerInfo(dataDict[name][gen]['data'], dataDict[name][0]['common'], args)
    return dataDict


def _addTautomerInfo(data, commonData, args):
    tinfo = {'prod': dict(), 'sbs': dict(), 'full': set(), 'rxnum': dict()}
    isEmpty = True
    if args.fixTautomerRatio == 'DE':
        raise NotImplementedError
    # xProd = args.setTautomerRatio / (1 + args.setTautomerRatio)
    for idrx in data['idToFullRx']:
        rxid = data['idToFullRx'][idrx]['rxid']
        if len(rxid) != 1:
            raise NotImplementedError
        rxid = rxid[0]
        if rxid not in consts.tautomerRxid:
            # tautomer rxid is read from global variable for now in the future will be improved
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
        hasSbsEnolC = sbsMol.HasSubstructMatch(consts.enolateC)
        hasSbsEnolO = sbsMol.HasSubstructMatch(consts.enolateO)
        isTautoEnolate = False
        if hasSbsEnolC and hasSbsEnolO:
            raise NotImplementedError
        elif hasSbsEnolC or hasSbsEnolO:
            prodMol = Chem.MolFromSmiles(prod)
            hasProdEnolC = prodMol.HasSubstructMatch(consts.enolateC)
            hasProdEnolO = prodMol.HasSubstructMatch(consts.enolateO)
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


# #############################


def loadBaseUsageInfo(basefile):
    try:
        fh = open(basefile)
        dct = dict()
        for line in fh:
            striped = line.strip()
            if not striped or striped.startswith('#'):
                continue
            molname, baseEq, baseName = line.split('\t')
            dct[molname] = {'base': baseName, 'eq': float(baseEq)}
        fh.close()
        return dct
    except:
        return dict()


def loadOverwritenInitC(dataDict, fname, args):
    fh = open(fname)
    retDict = dict()
    for line in fh:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        try:
            graphName, cmdSmiles, conc = line.split('\t')
        except:
            print(f"Wrong format in line {line}")
            raise consts.WrongData
        if graphName not in retDict:
            retDict[graphName] = dict()
        canonSmiles = Chem.CanonSmiles(cmdSmiles)
        if not canonSmiles:
            raise consts.WrongData
        if canonSmiles in retDict[graphName]:
            raise consts.WrongData
        try:
            conc = float(conc)
        except:
            raise consts.WrongData
        if conc < 0:
            raise consts.WrongData
        retDict[graphName][canonSmiles] = conc
    fh.close()
    if args.debug:
        print(f"initial concentration overwriten for {retDict}")
    return retDict
