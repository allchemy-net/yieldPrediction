import math
import networkx
# from networkx.drawing.nx_agraph import write_dot
# from numba import jit
# import scipy.stats
from rdkit import Chem
from cachetools import cached
from cachetools.keys import hashkey
import helpers, consts, calcScaling, calcHelper
try:
    import calcCopasi
except ModuleNotFoundError:
    print("copasi mode not available")

# orion is needed for hyperparam opt, you can install it with:
# pip install orion
# if above command failed try to update pip
# pip3.9 install --upgrade pip
# if you have problem like
# [...]
# FileNotFoundError: [Errno 2] No such file or directory: '/usr/bin/pip3.9'
# [....]
# do make symlink:
# ln -s /usr/local/bin/pip3.9 /usr/bin/pip3.9


class WrongData(Exception):
    pass


class WrongOptions(Exception):
    # when combination of used cmd arguments has no sense
    pass


class NegativeMassException(Exception):
    # when after step calculation one of mass is negative
    pass


#################################
def getActiveRxesAndSbs(rxes, graph, mass, addedRx, args):
    actRx = []
    rxSbs = []
    keepCmd = set()
    for rx in rxes:
        if rx in addedRx:
            if args.debug:
                print(f":79::RX {rx} in added IGNORE")
            continue
        sbses = tuple(graph.predecessors(rx))
        if args.debug:
            print(f":83::RX {rx} ANALIZE: {sbses}")
        if not sbses:
            print("---RX", rx, "has no substrates")
            raise
        # seenSbses = [smi for smi in sbses if smi in mass]
        unseenSbses = [smi for smi in sbses if smi not in mass]
        if not unseenSbses:
            actRx.append(rx)
            addedRx.add(rx)
            rxSbs.append(sbses)
            if args.debug:
                print(f"  OK {rxSbs} ACT: {actRx}")
        else:
            if args.debug:
                print(":97:::UNSEEN", unseenSbses, [e for e in graph.edges if set(e).intersection(set(unseenSbses))])
            keepCmd.update(set(unseenSbses))
    if args.debug:
        print(f"ACTRX:66: after itr: {actRx} \n")
    return actRx, rxSbs, keepCmd, addedRx


def getAllActiveRxesAndSbs(graph, mass, rxToSubstrates, args):
    actRx = []
    for rx in rxToSubstrates:
        unseenSbses = [smi for smi in rxToSubstrates[rx] if smi not in mass or mass[smi] < args.epsilon]
        if not unseenSbses:
            actRx.append(rx)
    if args.debug:
        print(f"ACTRX:77: after itr: {actRx} \n")
    return actRx


# @jit(nopython=True)
def makePropagation(data, mass, front, target, adbData, commonData, calcInfo, args):
    prevFront = None
    addedRx = set()
    fullData = []
    # print("DATAK", type(data['graph']), data['data'].keys())
    massInRx = {rx: [] for rx in data['data']['rxtupleToId'].values()}
    rev = {f'{rx}_REV': [] for rx in data['data']['rxtupleToId'].values()}
    massInRx.update(rev)
    if args.propagationMode == 'front':
        if calcInfo:
            raise NotImplementedError
        for itr in range(args.maxiter):
            mass, front, addedRx = makePropStepFront(data['graph'], mass, front, addedRx, data['data'], commonData, itr, target, adbData, args)
            # if args.saveDataForPlot:
            fullData.append(mass.copy())
            if args.debug:
                print(f":102:ITR {itr} RX::: {addedRx} M::: {mass} :::FRONT {front}")
            if front == prevFront:
                print("\n\n=============PROBLEM================")
                args.debug = 10
                mass, front, addedRx = makePropStepFront(data['graph'], mass, front, addedRx, data['data'], commonData, itr, target, adbData, args)
                print("EDH", data['graph'].edges)
                raise
            if not front:
                break
            prevFront = front
        else:
            print("not converged!!!!")
            getState(mass, addedRx, front, data)
            raise
        return mass, fullData, massInRx
    elif args.propagationMode == 'all':
        targetMass = []
        rxToSubstrates = _getSubstrateForRxInGraph(data)
        if args.debug > 10:
            print("RX TO SBS", rxToSubstrates)
        numIterToPerform = args.maxiter
        if args.includeConditions != 'no':
            if args.includeConditions in {'minimal', 'minimalMild'}:
                if args.debug:
                    print(f"minimalCONDITION {len(commonData['allowedCond']['minimalVariants'])} :0: {commonData['allowedCond']['minimalVariants'][0]}")
                numIterToPerform *= len(commonData['allowedCond']['minimalVariants'][0])
            else:
                raise NotImplementedError
        for itr in range(numIterToPerform):
            mass, changedRx = makePropStepAll(data['graph'], mass, data['data'], target, rxToSubstrates, adbData, commonData, itr, calcInfo, args)
            for rx in changedRx:
                massInRx[rx].append(changedRx[rx])
            # if args.saveDataForPlot:
            fullData.append(mass.copy())
            if args.debug > 10:
                print("SUMMASS", sum(mass.values()))
            targetMass.append(mass.get(target, 0))
            calcInfo = calcHelper.updateCalcInfo(calcInfo, fullData, args)
            if args.debug:
                print(f"::141::ITR {itr} RX::: {addedRx} M::: {mass} :::FRONT {front}")
        massDict = {target: targetMass}
        for smi in mass:
            if smi == target:
                continue
            massDict[smi] = mass[smi]
        if args.debug > 123:
            print("MASS IN RX", {key: (max(massInRx[key]) if massInRx[key] else massInRx[key], len(massInRx[key])) for key in massInRx})
        return massDict, fullData, massInRx
    else:
        raise NotImplementedError


def _getSubstrateForRxInGraph(data):
    rxnodeNames = list(data['data']['idToFullRx'].keys())
    retDict = dict()
    for rxnum in rxnodeNames:
        try:
            retDict[rxnum] = tuple(data['graph'].predecessors(rxnum))
            # print(f"DATA_GRAPH: {rxnum} {retDict[rxnum]}, PRD: {tuple(data['graph'].successors(rxnum))}")
        except networkx.exception.NetworkXError:
            # if args.disableLongCut == 0:
            #     raise NotImplementedError
            continue
    return retDict


# @jit(nopython=True)
def makePropStepFront(graph, mass, front, addedRx, data, commonData, itr, target, adbData, args):
    if args.debug > 10:
        print("DEBUG::174 begin propagation step")
    nextFront = set()
    activeRx = []  # list of rx node
    rxSbses = []   # list of substrateslist
    for node in front:
        rxes = graph.successors(node)
        actRxes, actSbs, keepCmd, addedRx = getActiveRxesAndSbs(rxes, graph, mass, addedRx, args)
        if keepCmd:
            nextFront.update(set(keepCmd))
        for poz, rx in enumerate(actRxes):
            activeRx.append(rx)
            rxSbses.append(actSbs[poz])
    # sbsToCount = countSubstratesFromListList(rxSbses)
    keyToFullRx = data['idToFullRx']
    sbsCoef, prodCoef = getSbsAndProdCoef(activeRx, keyToFullRx, args)
    if args.debug:
        print(f"130:SBSC: {sbsCoef} PROD: {prodCoef}")
    rxYields = getRxYields(activeRx, keyToFullRx, target, adbData, commonData, itr, sbsCoef, mass, args)
    sbsUsage = getIzolatedSubstrateUsage(activeRx, keyToFullRx, sbsCoef, rxYields, mass)
    if args.debug > 120:
        print(f"RXYield:: {rxYields}")
        print(f"SBSUSAGE {sbsUsage} M:: {mass}")
    sbsUsage = getNormalizedSubstrateUsage(sbsUsage, sbsCoef, mass, args)
    if args.debug:
        print("inSTEP: sbsUsage:", sbsUsage, "Mass:", mass)
        print("SBSusage", sbsUsage, "YLD", rxYields, "RX", activeRx)
    # mass = doPropagateMass(propagFactor, activeRx, rxSbses, mass, data['idToFullRx'], args)
    mass, changedRx = performMassTransfer(sbsUsage, activeRx, sbsCoef, prodCoef, mass, sbsCoef, args)
    nextCmds = getProductsFromRx(activeRx, graph)
    nextFront.update(set(nextCmds))
    if 'tautomers' in data or args.fixTautomerRatio != 'no':
        raise NotImplementedError
    if args.debug:
        print(f"::132 nextfr {nextFront} init front {front} ADDEDRX {addedRx} MASS:: {mass}")
    return mass, nextFront, addedRx


# @jit(nopython=True)
def makePropStepAll(graph, mass, data, target, rxToSubstrates, adbData, commonData, itr, calcInfo, args):
    # rx -> sbs
    if args.debug > 10:
        print(f"=== DEBUG::174 begin propagation step {itr} ===")
    activeRx = getAllActiveRxesAndSbs(graph, mass, rxToSubstrates, args)
    if data['rxAvail']['gen'] or calcInfo['disableRx']:
        bannedSet = set()
        if data['rxAvail']['gen']:
            gen = 1 + (itr // args.maxiter)
            bannedRx = [rx for rx in data['rxAvail']['gen'] if data['rxAvail']['gen'][rx]['min'] and data['rxAvail']['gen'][rx]['min'] > gen]
            bannedSet.update(set(bannedRx))
        if calcInfo.get('disableRx'):
            bannedSet.update(calcInfo['disableRx'])
        if args.debug > 150:
            print("ACT", activeRx, bannedSet)
        activeRx = [rx for rx in activeRx if rx not in bannedSet]
        if args.debug > 150:
            print("TR", activeRx)

    if args.debug > 20:
        print("STEP", len(activeRx), len(rxToSubstrates))
    keyToFullRx = data['idToFullRx']
    sbsCoef, prodCoef = getSbsAndProdCoef(activeRx, keyToFullRx, args)
    if args.debug:
        print(f"130:SBSC: {sbsCoef} PROD: {prodCoef}")
    rxYields = getRxYields(activeRx, keyToFullRx, target, adbData, commonData, itr, sbsCoef, mass, args)
    sbsUsage = getIzolatedSubstrateUsage(activeRx, keyToFullRx, sbsCoef, rxYields, mass)
    if args.debug > 10:
        print(f"SBSUsage::a:: {itr} USE: {sbsUsage} MASS: {mass} COEF: {sbsCoef}")
        print(f"RXYIELD:: {itr} YY: {rxYields} {activeRx}")
    sbsUsage = getNormalizedSubstrateUsage(sbsUsage, sbsCoef, mass, args)
    if args.debug > 10:
        rxToSmi = {rx: data['idToFullRx'][rx]['smiles'] for rx in data['idToFullRx']}
        print(f"inSTEP: sbsUsage: {sbsUsage} yld: {rxYields} ACTrx: {activeRx} Mass: {mass} F: {rxToSmi}")
        print("SBSusage", sbsUsage, "YLD", rxYields, "RX", activeRx)
    # mass = doPropagateMass(propagFactor, activeRx, rxSbses, mass, data['idToFullRx'], args)
    mass, changedRx = performMassTransfer(sbsUsage, activeRx, sbsCoef, prodCoef, mass, args)
    # nextCmds = getProductsFromRx(activeRx, graph)
    # print(f"ITR {itr} {gen} {bannedRx} {activeRx} RX: {changedRx} {sbsUsage}")
    if 'tautomers' in data:
        mass, changedRx = _setTautomerRatio(mass, changedRx, data['tautomers'], activeRx, args)
    if 'acidBase' in data:
        mass, changedRx = _setAcidBaseRatio(mass, changedRx, data['acidBase'], commonData, itr, args)
    return mass, changedRx


def _setTautomerRatio(mass, changedRx, tautoData, activeRx, args):
    for actRx in activeRx:
        if actRx not in tautoData['rxnum']:
            continue
        sbs, prod, _, _, whichToCalcRatio, isEnolate = tautoData['rxnum'][actRx]
        totalMass = mass[sbs] + mass[prod]
        if args.fixTautomerRatio == 'constEnolate' and isEnolate:
            xProd = args.setTautoEnolateRatio / (1 + args.setTautoEnolateRatio)
        else:
            xProd = args.setTautomerRatio / (1 + args.setTautomerRatio)
        if whichToCalcRatio == 0:
            xProd = 1 - xProd
        prodMass = totalMass * xProd
        sbsMass = totalMass - prodMass
        # diffProd = prodMass - mass[prod]
        # diffSbs = sbsMass - mass[sbs]
        mass[sbs] = sbsMass
        if mass[sbs] < 0:
            raise NegativeMassException
        mass[prod] = prodMass
        if mass[prod] < 0:
            raise NegativeMassException
        if args.debug > 120:
            print(f"TAUTO_RATIO sbs: {sbs} {mass[sbs]}  prod: {prod} {mass[prod]} xProd: {xProd} TDATA: {tautoData['rxnum'][actRx]}")
        # if actRx in changedRx:
        #    print(changedRx[actRx], diffSbs, diffProd, "SBS PR", mass[sbs], mass[prod])
        #    raise NotImplementedError
    return mass, changedRx


# acid base eq
def _getCondAorB(commonData, itr, args):
    condAorB = None
    if args.fixAcidBaseRatio == 'cond' or args.fixAcidBaseRatio == 'condIgnoreWeak':
        condNum = itr // args.maxiter
        conds = commonData['allowedCond']['minimalVariants'][0][condNum][0]
        if args.fixAcidBaseRatio == 'cond':
            if conds.intersection({'B', 'WB', 'SB'}):
                if not conds.intersection({'A', 'SA', 'WA'}):
                    condAorB = 'B'
            if conds.intersection({'A', 'SA', 'WA'}):
                if not conds.intersection({'SB', 'B', 'WB'}):
                    condAorB = 'A'
        elif args.fixAcidBaseRatio == 'condIgnoreWeak':
            if conds.intersection({'B', 'SB'}):
                if not conds.intersection({'A', 'SA'}):
                    condAorB = 'B'
            if conds.intersection({'A', 'SA'}):
                if not conds.intersection({'SB', 'B'}):
                    condAorB = 'A'
        else:
            raise NotImplementedError
    return condAorB


def _setAcidBaseRatio(mass, changedRx, dataAB, commonData, itr, args):
    return _setAcidBaseRatioSimple(mass, changedRx, dataAB, commonData, itr, args)
    # complex mode seems to be more reasonable but better correlation was obtained using simple mode
    # therefore we use simple mode however complex mode will stay for future use/tests
    # return _setAcidBaseRatioComplex(mass, changedRx, dataAB, commonData, itr, args)


def _setAcidBaseRatioSimple(mass, changedRx, dataAB, commonData, itr, args):
    massCmdSet = set(mass.keys())
    changedWild = dataAB['wildForms'].intersection(massCmdSet)
    if not changedWild:
        return mass, changedRx
    condAorB = _getCondAorB(commonData, itr, args)
    for idx, smilesWild in enumerate(changedWild):
        mainForm = commonData['acidBase']['wildToMain'][smilesWild]
        if mainForm not in mass:
            # print(f"MAIN FORM NOT FOUND! WILD: {smilesWild} MAIN: {mainForm} :: ")
            continue
        xOfMainForm = 1 - (args.setAcidBaseRatio / (1 + args.setAcidBaseRatio))
        formType = None
        if condAorB:
            formType = commonData['acidBase']['wildFormType'][smilesWild]
            # print("WILD FORM TYPE", smilesWild, formType)
            if 'La' in smilesWild and 'La' not in mainForm:
                pass
            elif formType != condAorB:  # acid form in base cond or vice versa
                xOfMainForm = 0.9999
        summMass = mass[mainForm] + mass[smilesWild]
        # xOfMainForm = dataAB['xOfMainForm'][smilesWild]
        if args.debug > 30:
            print(f"{idx} ACIDBASE:S M: {mainForm} W: {smilesWild} mass: {mass[mainForm]} {mass[smilesWild]} xmain: {xOfMainForm} TYPE: {formType} AB: {condAorB}")
        mass[mainForm] = xOfMainForm * summMass
        if mass[mainForm] < 0:
            raise NegativeMassException
        mass[smilesWild] = summMass - mass[mainForm]
        if mass[smilesWild] < 0:
            raise NegativeMassException
        # TO DO  do mod changedRX
    return mass, changedRx


def _setAcidBaseRatioComplex(mass, changedRx, dataAB, commonData, itr, args):
    massCmdSet = set(mass.keys())
    changedWild = dataAB['wildForms'].intersection(massCmdSet)
    if not changedWild:
        return mass, changedRx
    # print("DATA AB", dataAB)
    condAorB = _getCondAorB(commonData, itr, args)
    mainToWilds = dict()
    for idx, smilesWild in enumerate(changedWild):
        mainForm = commonData['acidBase']['wildToMain'][smilesWild]
        if mainForm not in mass:
            # print(f"MAIN FORM NOT FOUND! WILD: {smilesWild} MAIN: {mainForm} :: ")
            continue
        if mainForm not in mainToWilds:
            mainToWilds[mainForm] = []
        mainToWilds[mainForm].append(smilesWild)
    #
    xOfMainFormRaw = 1 - (args.setAcidBaseRatio / (1 + args.setAcidBaseRatio))
    for mainForm in mainToWilds:
        summMass = mass[mainForm] + sum([mass[smilesWild] for smilesWild in mainToWilds[mainForm]])
        # xOfMainForm = dataAB['xOfMainForm'][smilesWild]
        bannedWildForm = []
        if condAorB:
            for smilesWild in mainToWilds[mainForm]:
                formType = commonData['acidBase']['wildFormType'][smilesWild]
                if formType != condAorB:  # acid form in base cond or vice versa
                    bannedWildForm.append(smilesWild)
        allowedWildForm = [smi for smi in mainToWilds[mainForm] if smi not in set(bannedWildForm)]
        if args.debug > 30:
            print(f"{idx} ACIDBASE:C {mainForm} {mainToWilds[mainForm]} mass: {mass[mainForm]} {[mass[smi] for smi in mainToWilds[mainForm]]}")
        numWilds = len(mainToWilds[mainForm])
        if bannedWildForm and allowedWildForm:
            xOfMainForm = xOfMainFormRaw
            mass[mainForm] = xOfMainForm * summMass
            summMassOfWild = summMass - mass[mainForm]
            summMassOfBanWild = summMassOfWild / 10_000
            summMassOfWildAllow = summMassOfWild - summMassOfBanWild
            numWildsBanned = len(bannedWildForm)
            numWildsAllowed = len(allowedWildForm)
            for smilesWild in bannedWildForm:
                mass[smilesWild] = summMassOfBanWild / numWildsBanned
            for smilesWild in allowedWildForm:
                mass[smilesWild] = summMassOfWildAllow / numWildsAllowed
        else:  # all wild allowed or all wild banned
            if not allowedWildForm:
                xOfMainForm = 0.9999
            else:
                xOfMainForm = xOfMainFormRaw
            mass[mainForm] = xOfMainForm * summMass
            summMassOfWild = summMass - mass[mainForm]
            for smilesWild in mainToWilds[mainForm]:
                mass[smilesWild] = summMassOfWild / numWilds
        if mass[mainForm] < 0:
            raise NegativeMassException
        for smilesWild in mainToWilds[mainForm]:
            if mass[smilesWild] < 0:
                raise NegativeMassException
        # TO DO  do mod changedRX
    return mass, changedRx


# @jit(nopython=True)
def performMassTransfer(sbsUsage, activeRx, sbsCoef, prodCoef, mass, args):
    # perform mass transfer id est substrat from substrate and add to product all in mass variable
    # print("\n\n11", sbsUsage, "22:", activeRx, "33:", sbsCoef, "44:", prodCoef, "555:", mass)
    changedMass = dict()
    changedRx = dict()
    for poz, num in enumerate(activeRx):
        # check stoichiometry
        stoichMass = [sbsUsage[poz][sbs] / sbsCoef[poz][sbs] for sbs in sbsUsage[poz]]
        if len(set(stoichMass)) != 1:
            if max(stoichMass) - min(stoichMass) > 1e-9:
                print("stoich", stoichMass)
                raise
        changedRx[num] = stoichMass[0]
        # newMass = {prod: stoichMass[0] * prodCoef[poz][prod] for prod in prodCoef[poz]}
        for prod in prodCoef[poz]:
            if prod not in changedMass:
                changedMass[prod] = 0
            changedMass[prod] += stoichMass[0] * prodCoef[poz][prod]
        for sbs in sbsCoef[poz]:
            if sbs not in changedMass:
                changedMass[sbs] = 0
            changedMass[sbs] -= stoichMass[0] * sbsCoef[poz][sbs]
        # for oneprod in newMass:
        #    if oneprod not in mass:
        #        mass[oneprod] = 0
        #    mass[oneprod] += newMass[oneprod]
    # print("--RX", num, sbsUsage[poz], sbsCoef[poz], prodCoef[poz])
    # set equilibrium between taumeric species
    for smi in changedMass:
        if smi not in mass:
            mass[smi] = 0
        mass[smi] += changedMass[smi]
        if mass[smi] < 0:
            if args.paropt != 'N':
                print("NEGATIVE MASS", smi, mass[smi], "ACTRX", activeRx, "SBS", sbsUsage, "change", changedMass[smi])
                print(f"ALLMASS: {mass} ALLCHanges: {changedMass} RXchanged:: {changedRx}")
            raise NegativeMassException
    return mass, changedRx


# @jit(nopython=True)
def getIzolatedSubstrateUsage(activeRx, keyToFullRx, sbsCoef, rxYields, mass):
    izolatedSubstrateUsage = []
    for rxpoz, rxkey in enumerate(activeRx):
        try:
            izolateMaxUsage = {sbs: (mass[sbs] * rxYields[rxpoz]) / sbsCoef[rxpoz][sbs] for sbs in sbsCoef[rxpoz]}
        except KeyError:
            print("\nPROBLEM", rxpoz, rxkey, sbsCoef[rxpoz], rxYields[rxpoz], mass, "ACT RX", activeRx)
            print("38rev", keyToFullRx[rxkey])
            raise
        minValUsage = min(izolateMaxUsage.values())
        singleRxSbsUsage = {sbs: (sbsCoef[rxpoz][sbs] * minValUsage) for sbs in sbsCoef[rxpoz]}
        izolatedSubstrateUsage.append(singleRxSbsUsage)
    return izolatedSubstrateUsage


# @jit(nopython=True)
def getNormalizedSubstrateUsage(sbsUsage, sbsCoef, mass, args):
    if args.debug > 10:
        print("DEBUG:::240 enter to normalize substrate usage")
    if args.debug > 20:
        print("DEBUG::241 sbs:", sbsUsage, sbsCoef, "mass", mass)

    # sbsUsage = [{'CC(C)C=O': 75.0}, {'CC(C)C=O': 75.0}, {'CC(C)C=O': 75.0, 'Nc1ccccc1C(=O)O': 75.0}, {'Nc1ccccc1C(=O)O': 75.0}]
    # sbsCoef = [{'CC(C)C=O': 1}, {'CC(C)C=O': 1}, {'CC(C)C=O': 1, 'Nc1ccccc1C(=O)O': 1}, {'Nc1ccccc1C(=O)O': 1}]
    # mass = {'[C-]#[N+]C1=CCCCC1': 100, 'COc1ccc(CN)cc1': 100, 'CC(C)C=O': 100, 'Nc1ccccc1C(=O)O': 100}

    totalUsage = dict()
    multiUsageSet = set()
    for poz, sbsUseInRx in enumerate(sbsUsage):
        for sbs in sbsUseInRx:
            if sbs not in totalUsage:
                totalUsage[sbs] = {'usage': sbsUseInRx[sbs], 'rxpoz': [poz, ]}
            else:
                totalUsage[sbs]['usage'] += sbsUseInRx[sbs]
                totalUsage[sbs]['rxpoz'].append(poz)
                multiUsageSet.add(sbs)
    if not multiUsageSet:
        return sbsUsage
    itr = 0
    if args.debug > 30:
        print("DEBUG::523:: totalUSAGE", totalUsage, "MULTI", multiUsageSet)
    while True:
        itr += 1
        scaleFactor, includedSbs, involvedRx = _getScalerForRx(sbsUsage, sbsCoef, mass, args.maxConsumption, totalUsage, multiUsageSet, args)
        if args.debug >= 1:
            print("SCALE:528:", scaleFactor, includedSbs, involvedRx)
        if not scaleFactor:
            if args.debug >= 1:
                print("BREAAAk")
            break
        sbsUsage, totalUsage = _scaleReactions(scaleFactor, sbsUsage, includedSbs, involvedRx, totalUsage)
        if args.debug >= 1:
            print(":535:SBS", sbsUsage, "TT", totalUsage)
        if itr > 10:
            print(":537:SBS", sbsUsage, "TT", totalUsage, scaleFactor)
            if itr > 15:
                raise
    return sbsUsage


# @jit(nopython=True)
def _getScalerForRx(sbsUsage, sbsCoef, mass, maxConsumption, totalUsage, multiUsageSet, args, eps=1e-6):
    if args.debug > 30:
        print("===== IN:: _getScalerForRx")
    includedSbs = set()
    involvedRx = set()
    activeRx = set()
    allScaleFactors = []
    for poz, sbsUseInRx in enumerate(sbsUsage):
        scaleFactors = []
        for sbs in sbsUseInRx:
            _ = [activeRx.add(poz) for poz in totalUsage[sbs]['rxpoz']]
            if sbs in includedSbs:
                continue
            if sbs not in multiUsageSet:
                continue
            availSbs = mass[sbs] * maxConsumption
            if availSbs >= totalUsage[sbs]['usage']:
                continue
            scaleFactor = availSbs / (totalUsage[sbs]['usage'] * sbsCoef[poz][sbs])
            scaleFactor -= eps
            includedSbs.add(sbs)
            _ = [involvedRx.add(poz) for poz in totalUsage[sbs]['rxpoz']]
            scaleFactors.append(scaleFactor)
        if scaleFactors:
            if args.debug:
                print("\n======", scaleFactors, poz, sbsUsage, "ALL", allScaleFactors)
            scaleFactor = min(scaleFactors)
            # sbsUsage, totalUsage = _scaleReactions(scaleFactor, sbsUsage, includedSbs, involvedRx)
            if len(scaleFactors) > 1 and args.debug > 5:
                print("need to check if it works correct, e.g. A+B>>, A+C>>  and B+D>> when A and B has different scaling factor")
                # in such case it can be mutually dependent i.e. apply scaling for A may lead to situation where scaling for B is not need
                # and vice versa so we need to scale scalling based on mutual realtion between set of reaction involving A nad B
                # it is better to analyze with some example
                # raise NotImplementedError
            allScaleFactors.append(scaleFactor)
    if allScaleFactors:
        scaleFactor = min(allScaleFactors)
    else:
        scaleFactor = None
    if args.debug > 10:
        print("RETURN", allScaleFactors, scaleFactor)
    if args.propagationMode == 'front':
        return scaleFactor, includedSbs, involvedRx
    elif args.propagationMode == 'all':
        return scaleFactor, includedSbs, activeRx
    else:
        raise NotImplementedError


# @jit(nopython=True)
def _scaleReactions(scaleFactor, sbsUsage, includedSbs, involvedRx, totalUsage):
    # print("OPTS 1:", scaleFactor, "2:", sbsUsage, "3:", includedSbs, "4:", involvedRx, "5::", totalUsage)
    for rxpoz in involvedRx:
        for sbs in sbsUsage[rxpoz]:
            prevUsage = sbsUsage[rxpoz][sbs]
            sbsUsage[rxpoz][sbs] = scaleFactor * prevUsage
            changeInUsage = prevUsage - sbsUsage[rxpoz][sbs]
            if changeInUsage < 0:
                print("CH", prevUsage, scaleFactor, sbsUsage[rxpoz][sbs], "CH", changeInUsage)
                raise
            # print("CH", changeInUsage, "SBS", sbsUsage)
            totalUsage[sbs]['usage'] -= changeInUsage
    return sbsUsage, totalUsage


def getProductsFromRx(activeRx, graph):
    products = set()
    for rx in activeRx:
        prods = graph.successors(rx)
        products.update(set(prods))
    return products


# @jit(nopython=True)
def getSbsAndProdCoef(activeRx, keyToFullRx, args):
    # return list of dict with stoichiometry coef for each substrate in each reaction
    sbsToCoef = []
    prodToCoef = []
    for rxpoz, rxkey in enumerate(activeRx):
        if args.debug > 1450:
             print("KEY:1708:", rxkey, rxpoz, keyToFullRx[rxkey]['fullRxSmiles'], keyToFullRx[rxkey])
        fullSbs, fullProd = keyToFullRx[rxkey]['fullRxSmiles'].split('>>')
        fullSbs = fullSbs.split('.')
        fullProd = fullProd.split('.')
        # sbses = tuple(sorted(fullSbs))
        # factor = propagFactor[rxpoz]
        #
        sbsToCoef.append({sbs: fullSbs.count(sbs) for sbs in fullSbs})
        prodToCoef.append({prd: fullProd.count(prd) for prd in fullProd})
    return sbsToCoef, prodToCoef


def getRxYields(activeRx, keyToFullRx, target, adbData, commonData, itr, sbsCoef, masses, args):
    yields = []
    relFast = None
    # print(activeRx)
    if len(activeRx) > 1 and args.heatFormationBoltzman:
        des = []
        for rxkey in activeRx:
            try:
                des.append(max(keyToFullRx[rxkey]['DE'], 1e-3))
            except:
                print("PROBLEM WTH", rxkey, keyToFullRx[rxkey])
                raise
        relFast = _partitioning(des, args)

    trueFullPoz = []  # when some of normal reaction has yield higher than allowed scale all reaction except this from list
    maxNormalYield = 0  # max normal yield used to determine scaling factor
    for rxpoz, rxkey in enumerate(activeRx):
        rxids = keyToFullRx[rxkey]['rxid']
        if args.useEvansPolanyi == 'Y':
            rxDE = keyToFullRx[rxkey].get('DE', 0)
            if abs(rxDE) > args.maxAllowDE or args.scaleByDE < 1e-10:
                rxDE = 0
        elif args.useEvansPolanyi == 'N':
            rxDE = 0
        else:
            raise NotImplementedError
        rxidsSet = set(rxids)
        if rxidsSet.intersection(adbData['veryFastRx']):
            yields.append(1.0)
            trueFullPoz.append(rxpoz)
            continue
        hasF = rxidsSet.intersection(adbData['fastRxSet'])
        hasS = rxidsSet.intersection(adbData['slowRxSet'])
        hasVS = rxidsSet.intersection(adbData['verySlowRxSet'])
        hasSpeeds = (hasF, hasS, hasVS)
        isSidemech = rxidsSet.intersection(adbData['lowYield'])
        if keyToFullRx[rxkey].get('isMutex', False):
            if args.mutexRxMode == 'VS':
                hasSpeeds = (False, False, True)  # Fast, Slow, VerySlow
            elif args.mutexRxMode == 'sidemech':
                isSidemech = True
                if args.debug > 15:
                    print("MUTEXSDEMECH", rxidsSet)
            else:
                raise NotImplementedError
        if args.redox == 'NAsidemech':
            if not commonData['mainPathInfo']['redox']:  # redox
                if rxidsSet.intersection(adbData['redoxred']) or rxidsSet.intersection(adbData['redoxox']):
                    isSidemech = True
        rate = 0  # ignored when 0 otherwise rate/step is multiply by rate
        if args.scaleByMayr:
            rate = keyToFullRx[rxkey].get('rate', 0)
            if abs(rate) > args.maxAllowMayrNE:
                if rate < 0:
                    rate = -1 * args.maxAllowMayrNE
                else:
                    rate = args.maxAllowMayrNE
            rate = calcScaling.getMayrScaller(rate, args)
        rate = calcScaling.getMayrScaller(rate, args)
        yld = None
        if 1 in rxidsSet:
            yld = args.mainFormYield
        elif isSidemech:
            # print(f"SIDEMECH {rxids} R:: {args.sidemechYield} {rate} {relFast} {hasSpeeds}")
            yld = calcScaling.scaleYield(args.sidemechYield, relFast, rxpoz, hasSpeeds, rate, args)  # args.scaleByADBspeed)
        elif isinstance(rxkey, str) and 'REV' in rxkey:
            if args.reverseSpeedMode == 'simple':
                yld = args.reverseSpeed
            elif args.reverseSpeedMode == 'full':
                # hasSpeedsRev = (False, False, False)
                hasSpeedsRev = hasSpeeds
                yld = calcScaling.scaleYield(args.reverseSpeed, relFast, rxpoz, hasSpeedsRev, rate, args)  # args.scaleByADBspeed)
            else:
                raise NotImplementedError
        else:
            yld = calcScaling.scaleYield(args.defaultYield, relFast, rxpoz, hasSpeeds, rate, args)  # args.scaleByADBspeed)
        initYld = yld
        if keyToFullRx[rxkey]['isSideFromTarget']:
            yld = args.fromTargetYield
        # scale by DE
        scallerDE = calcScaling.getDEscaller(rxDE, args)
        yld *= scallerDE
        # scale by boost
        scallerBoost = calcScaling.getBoostScaler(keyToFullRx[rxkey], args)
        yld *= scallerBoost
        # scale by stoich
        stoichScaller = calcScaling.getStoichScaler(keyToFullRx[rxkey], masses, sbsCoef[rxpoz], args)
        yld *= stoichScaller
        # scale based on water
        if args.scaleWaterRx != 'no' and 'O' in sbsCoef[rxpoz] and not commonData['mainPathInfo']['isWaterCondNeeded']:
            yld /= args.waterRxScaleDesc
        if args.redox == 'rmRedOxWater':
            if not commonData['mainPathInfo']['redox'] and rxidsSet.intersection(adbData['redoxox']):
                yld /= args.waterRxScaleDesc
        # scale by conditions
        if args.includeConditions != 'no':
            if args.includeConditions in {'minimal', 'minimalMild'}:
                condNum = itr // args.maxiter
                conds = commonData['allowedCond']['minimalVariants'][0][condNum]
                # raise
                # if itr % args.maxiter == 0:
                #    print("ITR", itr, rxids, condNum, conds)
                scaler = _getScalerMinimConditions(commonData, rxids, conds, adbData, args)
                # scaler = keyToFullRx[rxkey].get('speedScale', 1)
                # if scaler != 1:
                #    print(f"SCLA {rxpoz} {rxkey} {rxids} ITR: {itr} {condNum} {conds} SCALE: {scaler}")
                if condNum in commonData['allowedCond']['minimalPozBeforeChangeToTarget']:
                    scaler *= _scaleCompetingWhichAreBeforeCondChange(condNum, keyToFullRx[rxkey]['smiles'], commonData)
                yld *= scaler
            else:
                raise NotImplementedError
        # scale by reversibility
        revRxkey = f'{rxkey}_REV'
        if revRxkey in keyToFullRx:
            yld *= args.forwardReversableScaler
        if yld > maxNormalYield:
            maxNormalYield = yld
        yields.append(yld)
        # print(f"YLD {rxids} {rxidsSet} {isSidemech} :Y: {yld} {maxNormalYield} INI {initYld} R {relFast} {rxpoz} {hasSpeeds} RAT {rate} X {xxx}")
    if maxNormalYield > args.maxStepYield:
        scalling = maxNormalYield / args.maxStepYield
        yields = [yld / scalling for yld in yields]
        for poz in trueFullPoz:
            yields[poz] = 1.0
    return yields


@cached(cache={}, key=lambda commonData, rxids, allowConds, adbData, args: hashkey(rxids[0], allowConds))
def _getScalerMinimConditions(commonData, rxids, allowConds, adbData, args):
    if 1 in rxids:
        return 1   # no scaling
    rxinf = adbData['rxes'][rxids[0]]
    rxAB, rxTemp = rxinf['rxacidbase'], rxinf['rxtemp']
    rxAB = set([cnd for cnd in rxAB.replace('LA', '').replace('+', '').split('.') if cnd])
    rxTemp = set([temp.replace('.', '').upper() for temp in rxTemp.split('.')])
    allowAB, allowTemp = allowConds
    if rxAB.intersection(allowAB) and rxTemp.intersection(allowTemp):
        if args.debug > 20:
            print("ALLOW", allowAB, allowTemp, "RX CONDS", rxAB, rxTemp, "11")
        scaller = 1
    else:
        if args.debug > 20:
            print("SCALED TO", allowAB, allowTemp, "RX CONDS", rxAB, rxTemp, args.scaleWrongConditionRx)
        scaller = args.scaleWrongConditionRx
    # if condNum in commonData['minimalPozBeforeChangeToTarget']:
    #    print("CHANGE", rxids, )
    #    raise
    return scaller


@cached(cache={}, key=lambda condNum, rxsmi, commonData: hashkey(condNum, rxsmi))
def _scaleCompetingWhichAreBeforeCondChange(condNum, rxsmi, commonData):
    # pri nt("RXS", rxsmi, type(rxsmi), rxsmi.split('>>'))
    sbs = rxsmi.split('>>')[0].split('.')
    sbs = set(sbs)
    if rxsmi == commonData['allowedCond']['steps'][condNum+1]:
        return 1
    if commonData['allowedCond']['minimalPozBeforeChangeToTarget'][condNum] not in sbs:
        return 1
    print("SCALD", rxsmi)
    return 0.00001


# @jit(nopython=True)
def _partitioning(DEs, args):
    # values in DEs expected in kJ/mol
    k = 8.31446  # J/(mol*K)
    T = args.temperature  # K
    maxde = max(DEs)
    minde = min(DEs)
    if maxde - minde > args.maxAllowDE:
        if args.debug > 10:
            print("DEBUG:519: DE", DEs)
        DEs = [1 for _ in DEs]
    deKTs = [-(de * 1000) / (k * T) for de in DEs]
    exps = [math.exp(deKT) for deKT in deKTs]
    sume = sum(exps)
    return [exp/sume for exp in exps]


# @jit(nopython=True)
def propagateNetwork(data, sbsStoich, adbData, commonData, args):
    if args.propagationMode == 'copasi':
        copasiRes, target = calcCopasi.propagateNetworkCopasi(data, sbsStoich, adbData, commonData, args)
        return calcCopasi.extractCopasiResults(copasiRes, target, args)
    else:
        return propagateNetworkOld(data, sbsStoich, adbData, commonData, args)


def makeABrxesAndHashes(smilesWild, mainForm, formType):
    if formType not in {'A', 'B', 'M'}:
        print(f"SMILES WILD {smilesWild} MAIN: {mainForm} FORM: {formType}")
        raise
    mainhash = f'{smilesWild}>>{mainForm}::MAIN'
    mainstr = f'{smilesWild} -> {mainForm}'
    revhash = f'{mainForm}>>{smilesWild}::{formType}_towild'
    revstr = f'{mainForm} -> {smilesWild}'
    return mainhash, mainstr, revhash, revstr


def propagateNetworkOld(data, sbsStoich, adbData, commonData, args):
    init = [x for x in data['graph'].nodes if data['graph'].in_degree(x) == 0]
    top = [x for x in data['graph'].nodes if data['graph'].out_degree(x) == 0]
    zeroGen = [smi for smi in data['data']['allSmiles'] if data['data']['allSmiles'][smi]['generation'] == 0]
    init = tuple(set(init).union(set(zeroGen)))
    calcInfo = calcHelper.formCalcInfo(data, args)
    try:
        target = data['data']['targetSmiles']
    except:
        print("DATA", data['data'].keys())
        raise
    if args.debug:
        print("GRAPH", data['graph'].edges)
        print(f"INIT {init} 0gen: {zeroGen}", "\nTOP:", top)
    mass = {smi: 100 for smi in init}
    if sbsStoich:
        s1 = set(sbsStoich.keys())
        s2 = set(init)
        if s1 != s2:
            if args.debug:
                print("missed in graph init:", s1 - s2, "missed in given sbs:", s2 - s1)
                print("Mass?", mass, "Stoich", sbsStoich)
                print("TARGET", target, "MISSED", '.'.join(s2 - s1))
            if s1 - s2:
                print("missed cmd", s1 - s2, s1, s2)
                raise NotImplementedError
        mass = {key: sbsStoich[key] for key in sbsStoich}  # make copy mass will be changed
        for smi in s2 - s1:
            # additives, probably was not used in reaction therefore we can set concentration close to zero
            mass[smi] = max(sbsStoich.values()) * args.additivesScale
    front = init
    if args.debug:
        print("INITial concentration", mass)
    mass, fullMass, rxChanged = makePropagation(data, mass, front, target, adbData, commonData, calcInfo, args)
    return mass, fullMass, rxChanged, target


def getState(mass, addedRx, front, data):
    print("D", data['data'].keys())
    notAddedRx = [rx for rx in data['data']['idToFullRx'] if rx not in addedRx]
    print("ID", data['data']['idToFullRx'].keys())
    print("NOT ADD", notAddedRx)
    for rx in notAddedRx:
        print("not added:", rx, data['data']['idToFullRx'][rx])


# ######################################################
# load pickle
# ##


def getReactiveOnMainPath(mainPathInfo, reactiveSmarts):
    print("MP", mainPathInfo)
    raise


def _addScalingForSideInOldConditions(data, args, commonData):
    mainPathRxSmiles = set([rxsmiAndRxid[0] for rxsmiAndRxid in commonData['allowedCond']['stepsVariants'][0]])
    forbidSbs = set(commonData['allowedCond']['minimalPozBeforeChangeToTarget'][0].values())
    for graphId in data['idToFullRx']:
        rxSmi = data['idToFullRx'][graphId]['smiles']
        if rxSmi in mainPathRxSmiles:
            continue
        sbs = set(rxSmi.split('>>')[0].split('.'))
        if sbs.intersection(forbidSbs):
            data['idToFullRx'][graphId]['speedScale'] = args.scaleSideInOldCondition
    return data['idToFullRx']


###############################################################
# d H formation


def _calcDE(rxsmiles, heatFormation):
    sbs, prod = rxsmiles.split('>>')
    sbs = sbs.split('.')
    prod = prod.split('.')
    pE = []
    sE = []
    missed = set()
    for p in prod:
        if p not in heatFormation:
            # print(f"NO heat formation for {p}")
            missed.add(p)
            continue
        pE.append(heatFormation[p])
    for s in sbs:
        if s not in heatFormation:
            print(f"NO heat formation for {s}")
            missed.add(s)
            continue
        sE.append(heatFormation[s])
    if missed:
        print(f"MISSED SMILES heat formation {'.'.join(missed)}")
        raise
    return sum(pE) - sum(sE)


def addReactionDE(dataDict, heatFormation):
    cache = dict()
    hasProblems = False
    values = []
    for name in dataDict:
        for gen in dataDict[name]:
            for rxid in dataDict[name][gen]['data']['idToFullRx']:
                rxsmiles = dataDict[name][gen]['data']['idToFullRx'][rxid]['fullRxSmiles']
                try:
                    if rxsmiles not in cache:
                        cache[rxsmiles] = _calcDE(rxsmiles, heatFormation)
                        values.append(cache[rxsmiles])
                        # print("RXVAL::", cache[rxsmiles], rxsmiles)
                    if abs(cache[rxsmiles]) > 20:
                        isRxOk, reason = checkIsRxOk(rxsmiles)
                        # if not isRxOk:
                        #    print("WARNING", reason, cache[rxsmiles], rxsmiles, dataDict[name][gen]['data']['idToFullRx'][rxid])
                    dataDict[name][gen]['data']['idToFullRx'][rxid]['DE'] = cache[rxsmiles]
                except:
                    dataDict[name][gen]['data']['idToFullRx'][rxid]['DE'] = 0.1
                    hasProblems = True
    if hasProblems:
        # raise
        pass
    return dataDict


def checkIsRxOk(rxsmiles):
    sbses, prods = rxsmiles.split('>>')
    sbses = sbses.split('.')
    prods = prods.split('.')
    sbsDict = _countAtoms(sbses)
    prodDict = _countAtoms(prods)
    if set(sbsDict.keys()) != set(prodDict.keys()):
        return False, 'missed atom type on one side'
    diff = []
    for atm in sbsDict:
        if sbsDict[atm] != prodDict[atm]:
            diff.append((atm, sbsDict[atm], prodDict[atm]))
    if diff:
        return False, diff
    return True, ""


def _countAtoms(smiList):
    atms = dict()
    for sbs in smiList:
        molH = Chem.AddHs(Chem.MolFromSmiles(sbs))
        for atm in molH.GetAtoms():
            sym = atm.GetSymbol()
            if sym not in atms:
                atms[sym] = 0
            atms[sym] += 1
    return atms


# ####################################
# tests
# ####


def testNormalizeSbsUsage():
    sbsUsage = [{'CC(C)C=O': 75.0}, {'CC(C)C=O': 75.0}, {'CC(C)C=O': 75.0, 'Nc1ccccc1C(=O)O': 75.0}, {'Nc1ccccc1C(=O)O': 75.0}]
    sbsCoef = [{'CC(C)C=O': 1}, {'CC(C)C=O': 1}, {'CC(C)C=O': 1, 'Nc1ccccc1C(=O)O': 1}, {'Nc1ccccc1C(=O)O': 1}]
    mass = {'[C-]#[N+]C1=CCCCC1': 100, 'COc1ccc(CN)cc1': 100, 'CC(C)C=O': 100, 'Nc1ccccc1C(=O)O': 100}
    print(sbsUsage, sbsCoef, mass)

##############################
# opt


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


############################################


def _getAllSmiles(dataDict):
    allSmiles = set()
    for calcName in dataDict:
        for gen in dataDict[calcName]:
            smiles = set(dataDict[calcName][gen]['data']['allSmiles'].keys())
            allSmiles.update(smiles)
    return allSmiles


def calcRxRate(allSmiles, mayrList):
    stats = dict()
    for smi in allSmiles:
        mol = Chem.MolFromSmiles(smi)
        matched = []
        for mayrDct in mayrList:
            # print(type(mayrDct), mayrDct)
            if any(mol.HasSubstructMatch(sma) for sma in mayrDct['smarts']):
                matched.append(mayrDct)
        if len(matched) > 10:
            smaset = set([m['smartsString'] for m in matched])
            if len(smaset) > 1:
                print(smi, smaset)
        if len(matched) not in stats:
            stats[len(matched)] = 0
        stats[len(matched)] += 1
        if not matched:
            print("NOMATCH", smi)
    print(stats)
