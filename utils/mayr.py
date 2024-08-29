import re, traceback, argparse
from rdkit import Chem
from rdkit.Chem import AllChem


class WrongData(Exception):
    pass


def parseArgs():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nucleodata", type=str, default='dataMayr/Nucleophiles.tsv', help='file with Mayr parametr for nucleophiles')
    ap.add_argument("--electrodata", type=str, default='dataMayr/Electrophiles.tsv', help='file with Mayr parametr for electrophiles')
    ap.add_argument("--rxsmiles", type=str, help='smiles of reaction', required=True)
    ap.add_argument('--rxsmarts', type=str, help='smarts of reaction', required=True)
    ap.add_argument('--rxpolar', type=str, help='defined allowed solvent polarity, P - means polar, Non-P - non polar, and P.Non-P means both', default="P.Non-P")
    ap.add_argument('--rxprotonic', type=str, help='define is protic (P), aprotic (AP) or both (P.AP) type of solvents allowed', default="P.AP",)
    ap.add_argument('--verbose', type=int, default=0, help='verbosity level')
    ap.add_argument('--constMayrSN', type=float, default=1.0, help='value for Mayr SN parameter')

    return ap.parse_args()


# helpers


def getUniqSmiles(molOrSmiles):
    if isinstance(molOrSmiles, str):
        smi = molOrSmiles
    else:
        smi = Chem.MolToSmiles(molOrSmiles)
    nonums = re.sub(':[0-9]*]', ']', smi)
    smi2 = Chem.CanonSmiles(nonums)
    return smi2


def getMappedSubstratesAndProducts(sbs, mappedSbs, prod, mappedProd):
    sbsMols = [Chem.MolFromSmiles(s) for s in sbs]
    mappedSbsAll = [Chem.MolFromSmiles(s) for s in mappedSbs]
    mappedSbs = _getMatchingToMapped(sbsMols, mappedSbsAll)
    if not mappedSbs:
        print("NO MATCHING", sbsMols, sbs, mappedSbsAll)
        raise NotImplementedError
    prodMols = [Chem.MolFromSmiles(s) for s in prod]
    mappedProdAll = [Chem.MolFromSmiles(s) for s in mappedProd]
    mappedProd = _getMatchingToMapped(prodMols, mappedProdAll)
    return mappedSbs, mappedProd


def _getMatchingToMapped(molList, mappedMolList):
    mapped = []
    usedPoz = set()
    for mol in molList:
        matches = _getMatchingMolFromList(mol, mappedMolList)
        if not matches:
            # print("MATCH", matches, "M::", Chem.MolToSmiles(mol), 'MM::', [Chem.MolToSmiles(m) for m in mappedMolList])
            raise NotImplementedError
        for match in matches:
            poz, mol = match
            if poz in usedPoz:
                continue
            usedPoz.add(poz)
            mapped.append(mol)
            break
        else:
            raise NotImplementedError
    return mapped


def getFragmentBasedSimilarity(mol, coreInMolIdx, coreMolList, opts):
    defaults = {'radius': 2, 'distance':  'tanimoto'}
    for key in defaults:
        if key not in opts:
            opts[key] = defaults[key]
    prodMol = None
    molFP = calcFragmentFingerprintForGivenIdxes(mol, coreInMolIdx, opts['radius'])
    refFP = []
    for smarts, molref in coreMolList:
        fp = calcFragmentFingerprintForGivenFrag(molref, smarts, prodMol, opts)
        refFP.append(fp)
    dsts = [tanimotoDistance(molFP, fp2) for fp2 in refFP]
    return dsts


def calcFragmentFingerprintForGivenFrag(mol, pattList, prodMol, opts):
    """
    Calculate fragment fingerprint (as a dict with counting)
    """
    coreIdxes = getRightFragmentIdxes(mol, pattList, prodMol)
    molFP = calcFragmentFingerprintForGivenIdxes(mol, coreIdxes, opts['radius'])
    return molFP


def getRightFragmentIdxes(mol, pattList, modmol, mode='fast'):
    """
    Find which patterns from pattList undergo reaction. Id est are present in mol but not in modmol.
    Arguments:
    - mol - mol object
    - pattList - list of patterns as mol objects
    - modmol - mol object of molecule without pattern (or None)
    - mode - 'fast' or 'strict'. Fast = when only one pattern match dont check modmol, strict = always check modmol (is not None)
    Return:
    list of indices of atom in mol object
    """
    matchesList = []
    matchesIdx = []
    foundPatt = []
    for idx, patt in enumerate(pattList):
        matches = mol.GetSubstructMatches(patt)
        if not matches:
            continue
        foundPatt.append(patt)
        if len(matches) > 1:
            ranks = tuple(Chem.CanonicalRankAtoms(mol, breakTies=False))
            rankOfMatches = [frozenset([ranks[idx] for idx in match]) for match in matches]
            if len(set(rankOfMatches)) == 1:
                matches = [matches[0], ]
            else:
                if not modmol:
                    print("MULTICORE::50", matches, "FOR MOL", Chem.MolToSmiles(mol), "::", [Chem.MolToSmarts(m) for m in foundPatt], pattList)
                    print("RANKS", rankOfMatches)
                    raise NotImplementedError
                matchesInMod = modmol.GetSubstructMatches(patt)
                if len(matches) - len(matchesInMod) > 1:
                    raise NotImplementedError
                # do compare matches
                raise NotImplementedError
        if len(matches) == 1:
            matchesList.append(matches)
            matchesIdx.append(idx)
        else:
            print("MULTIMATCH", matches, patt)
    if len(matchesIdx) > 1:
        print("MATCHES", matchesIdx, "FOR MOL", Chem.MolToSmiles(mol), "::", [Chem.MolToSmarts(m) for m in foundPatt], pattList)
        print("LIST", matchesList)
        raise NotImplementedError
    return matchesList[0][0]


def calcFragmentFingerprintForGivenIdxes(mol, idxes, radius):
    bi = dict()
    idxes = set(idxes)
    AllChem.GetMorganFingerprint(mol, radius, bitInfo=bi)
    fpc = dict()
    for hsh in bi:
        isin = 0
        for idx in bi[hsh]:
            if idxes.intersection(set(idx)):
                isin += 1
        if isin > 0:
            fpc[hsh] = isin
    return fpc


def _getMatchingMolFromList(mol, mappedMolList):
    matched = []
    for idx, mapped in enumerate(mappedMolList):
        if mol.HasSubstructMatch(mapped) and mapped.HasSubstructMatch(mol):
            matched.append((idx, mapped))
    return matched


def tanimotoDistance(fp1, fp2):
    sfp1 = set(fp1.keys())
    sfp2 = set(fp2.keys())
    sumFP = sfp1.union(sfp2)
    intersectFP = sfp1.intersection(sfp2)
    return len(intersectFP) / len(sumFP)


# ########
# load mayr


def parseMayrData(args):
    nucleoList, dbNucleo = _parseMayrExpFile(args.nucleodata, args)
    electroList, dbElectro = _parseMayrExpFile(args.electrodata, args)
    dictOfSmiles = dict()
    for smi in dbNucleo:
        dictOfSmiles[smi] = {'nucleo': dbNucleo[smi]}
    for smi in dbElectro:
        if smi not in dictOfSmiles:
            dictOfSmiles[smi] = {'electro': dbElectro[smi]}
        else:
            dictOfSmiles[smi]['electro'] = dbElectro[smi]
    return {'nucleo': nucleoList, 'electro': electroList, 'dict': dictOfSmiles}


def _parseMayrExpFile(fn, args):
    fh = open(fn)
    lines = fh.readlines()
    fh.close()
    header = {name: poz for poz, name in enumerate(lines[0].strip().split('\t'))}
    dataList = []
    dbDict = dict()
    smiPoz, solventPoz, smaPoz, reaPoz = header['SMILES'], header['Solvent'], header['SMARTS'], header['Reactivity Parameters']
    for lid, line in enumerate(lines[1:]):
        line = line.split('\t')
        try:
            data = _makeDataEntry(line[smiPoz], line[smaPoz], line[solventPoz], line[reaPoz])
            dataList.append(data)
            if data['smilesString'] not in dbDict:
                dbDict[data['smilesString']] = dict()
            # if data['solvent'] in dbDict[data['smilesString']]:
            #    print("SMI", data['smilesString'], dbDict[data['smilesString']])
            #    print("SMI LIN", line, data)
            #    # raise
            dbDict[data['smilesString']][data['solvent']] = data
        except WrongData:
            if args.verbose > 0:
                print("Mayr WRONG LNE", line)
            # raise
        except:
            print("Mayr PROBLEM WITH LINE", lid)
            raise
    return dataList, dbDict


def _makeDataEntry(smiles, smarts, solvent, reactivity):
    smilesObj = Chem.MolFromSmiles(smiles)
    if not smilesObj:
        if args.verbose > 0:
            print("SMI OBJ", smiles, "CANNOT READ !!!!")
        raise WrongData
    smarts = smarts.strip()
    smaStrList = []
    smaObjList = []
    if smarts:
        for sma in smarts.split('.'):
            smaMol = Chem.MolFromSmarts(sma)
            if not smaMol:
                print("INCORRECT SMARTS", sma)
                raise WrongData
            if not smilesObj.HasSubstructMatch(smaMol):
                sma = Chem.MolToSmarts(smaMol)
                smaMol = Chem.MolFromSmarts(sma)
                smiles = Chem.MolToSmiles(smilesObj)
                smilesObj = Chem.MolFromSmiles(smiles)
            if smilesObj.HasSubstructMatch(smaMol):
                smaStrList.append(sma)
                smaObjList.append(smaMol)
        if not smaStrList:
            # print("!!nor smaStrLisr: SMI ", smiles, orygSmi, " wrongs SMA", smarts)
            raise WrongData
        # if len(smaStrList) > 1:
        #    print("SMI", smiles, "MATCH TO::", smaStrList,)
    else:
        print("NO SMARTS")
        smaStrList = []
        smaObjList = []
    allowedSolvents = {'acetone', 'aq acetone', 'aq MeCN', 'water', 'DMF', 'DMSO', 'MeCN', 'dichloromethane',
                       'EtOH', 'iPrOH', 'nPrOH', 'MeOH', '0', '',
                       'water-MeCN mix', 'MeOH-MeCN mix', 'EtOH-MeCN mix', 'iPrOH-MeCN mix', 'nPrOH-MeCN mix',
                       'THF', 'TFE'}
    if solvent not in allowedSolvents:
        print("SOLV", solvent)
        raise WrongData
    # 'N Parameter: 23.74 s N Parameter: 0.71'
    reactivity = reactivity.strip().split()
    sN = None
    if len(reactivity) == 7:
        NorE = float(reactivity[2])
        sN = float(reactivity[6])
    elif len(reactivity) == 3:
        # E Parameter: 7.52
        NorE = float(reactivity[2])
    else:
        print("Cannot read reactivity information")
        raise WrongData
    retDict = {'smilesString': smiles, 'smiles': smilesObj, 'smartsString': tuple(smaStrList), 'smarts': tuple(smaObjList),
               'solvent': solvent, 'reactivity': NorE}
    if sN:
        retDict['sN'] = sN
    return retDict


# ####################
# calculate


def calcMayrRate(mayrData, args):
    allSbsMapped, allProdMapped = args.rxsmiles.split('>>')
    sbs = [getUniqSmiles(smi) for smi in allSbsMapped.split('.')]
    prod = [getUniqSmiles(smi) for smi in allProdMapped.split('.')]
    allSbsMapped = allSbsMapped.split('.')
    allProdMapped = allProdMapped.split('.')

    try:
        mappedSbs, mappedProd = getMappedSubstratesAndProducts(sbs, allSbsMapped, prod, allProdMapped)
    except:
        print("PROBLEM with RX::  SBS", sbs, allSbsMapped, "\nPRD", prod, allProdMapped)
        raise

    if not mappedSbs:
        print("no mapped sbs: ===>", mappedSbs, "PROD:", mappedProd)
        raise NotImplementedError
    mappedSbs, reactingCores, isIntramol = _getReactingCores(mappedSbs, mappedProd, args)
    if not (reactingCores and reactingCores[0]):
        print("\n====NO REACTING CORE!!!!  for", mappedSbs, mappedProd)

    molNEscores, isFull = _getMatchingMayr(mappedSbs, reactingCores, mayrData, isIntramol, args)
    rate = None
    if isFull:
        # print("RXSMI", rxsmiles)
        try:
            rate = _calcRateFromScores(molNEscores, args)
        except:
            print("RATE CALC FAILED", args)
            print("TRACE::1462", traceback.format_exc())
            rate = None
    else:
        print("Mayr NoFullDataFor RX", *[bool(dct['nucleo']) + bool(dct['electro']) for dct in molNEscores], sep='\t')
        pass
    if args.verbose:
        print("Mayr RATE::", sbs, prod, rate)
    # rate = fragmentBasedFingerprint.getReactionRate(sbs, mappedSbs, prod, mappedProd, mayrData, opts)
    return rate


def _getReactingCores(mappedSbsMolList, mappedProdMolList, args, returnMapNums=False, debug=False):
    # mappedSbsMolList order does not necessary reflect order of substrates in rxInfo['rxSmarts']
    # so
    rxnObj = AllChem.ReactionFromSmarts(args.rxsmarts)
    sbsPattMols = tuple(rxnObj.GetReactants())
    isIntramol = False
    if len(mappedSbsMolList) == 1 and len(sbsPattMols) == 2:
        if args.verbose > 3:
            print("INTRAMOL")
        isIntramol = True
    coresList = []
    for poz, sbsPatt in enumerate(sbsPattMols):
        if isIntramol:
            poz = 0
        sbsMol = mappedSbsMolList[poz]
        cores = sbsMol.GetSubstructMatches(sbsPatt)
        if returnMapNums:
            cores = [_convertIdxToMapNum(sbsMol, idxes) for idxes in cores]
        coresList.append(cores)
    return mappedSbsMolList, coresList, isIntramol


def _orderSbsToRx(sbs, rxInfo):
    if len(sbs) == 1:
        return sbs
    mols = [Chem.MolFromSmiles(smiles) for smiles in sbs]
    rxn = AllChem.ReactionFromSmarts(rxInfo['rxSmarts'])
    pattToMol = []
    hasProblem = False
    for sbsPatt in rxn.GetReactants():
        patts = [poz for poz, mol in enumerate(mols) if mol.HasSubstructMatch(sbsPatt)]
        if not patts:
            hasProblem = True
        pattToMol.append(patts)
    if hasProblem:
        hasProblem = False
        pattToMol = []
        allPatts = tuple(rxn.GetProducts())
        for sbsPatt in allPatts:
            patts = [poz for poz, mol in enumerate(mols) if mol.HasSubstructMatch(sbsPatt)]
            if not patts:
                hasProblem = True
            pattToMol.append(patts)
        if hasProblem:
            hasProblem = False
            prodSma = rxInfo['rxSmarts'].split('>>')[1]
            if prodSma[0] == '(' and prodSma[-1] == ')':
                patts = [Chem.MolFromSmarts(patt) for patt in prodSma[1:-1].split('.')]
                # if len(patts) != len(mols):
                #    print("RXSMA", rxInfo['rxSmarts'], "==", sbs)
                #    raise NotImplementedError
                for sbsPatt in patts:
                    patts = [poz for poz, mol in enumerate(mols) if mol.HasSubstructMatch(sbsPatt)]
                    if not patts:
                        hasProblem = True
                    pattToMol.append(patts)
            else:
                raise NotImplementedError
    if not pattToMol or hasProblem:
        print("HAS PROBLEM???", hasProblem, sbs, "PATT", rxInfo['rxSmarts'])
        raise NotImplementedError
    # check if is normal:
    if all([poz in matches for poz, matches in enumerate(pattToMol)]):
        if not sbs:
            raise NotImplementedError
        return sbs
    newOrder = []
    usedPoz = set()
    for sbsIdxList in pattToMol:
        for sbsIdx in sbsIdxList:
            if sbsIdx in usedPoz:
                continue
            usedPoz.add(sbsIdx)
            newOrder.append(sbs[sbsIdx])
    # print("NEW ORDER", newOrder, "FROM", sbs, "M", pattToMol)
    if not newOrder:
        print(f"???SBB order??? {[Chem.MolToSmiles(m) for m in sbs]}")
        raise NotImplementedError
    return newOrder


def _getMatchingMayr(mappedSbs, coresSbs, mayrData, isIntramol, args):
    bestMayr = []
    isFull = True
    hasN = False
    hasE = False
    uniqSmiles = [getUniqSmiles(m) for m in mappedSbs]
    for idx, core in enumerate(coresSbs):
        if isIntramol:
            idx = 0
        mol = mappedSbs[idx]
        if not core:
            if args.verbose > 0:
                print("!core not detected! Mayr CORE", core, "CORESBS", coresSbs, idx, uniqSmiles)
        try:
            bestMatchNucleo = _getMatchingNorE(mayrData['nucleo'], core, mol)
            bestMatchElectro = _getMatchingNorE(mayrData['electro'], core, mol)
        except:
            bestMatchNucleo, bestMatchElectro = None, None
            print("TRACE::1252", traceback.format_exc())
        if bestMatchNucleo:
            hasN = True
        elif bestMatchElectro:
            hasE = True
        else:
            isFull = False
        # bestMayr.append(bestMatch)
        bestMayr.append({'nucleo': bestMatchNucleo, 'electro': bestMatchElectro})
        # print("MAYR!!", uniqSmiles[idx], 'nucleo', bestMatchNucleo, 'electro', bestMatchElectro, sep=';')
    isTrullyFull = isFull and hasN and hasE
    if args.verbose > 5:
        print(f"MAYR ANION SUM: {isFull} {hasN} {hasE}  {bestMayr} ")
    return bestMayr, isTrullyFull


def _getMatchingNorE(mayrList, coresList, mol):
    matched = []
    debug = False
    for mayrDct in mayrList:
        # pri nt(type(mayrDct), mayrDct)
        # if any(mol.HasSubstructMatch(sma) for sma in mayrDct['smarts']):
        for pattObj in mayrDct['smarts']:
            for match in mol.GetSubstructMatches(pattObj):
                overlaps = [set(match).intersection(core) for core in coresList]
                nonEmptyOverlaps = [st for st in overlaps if st]
                if debug:
                    print("DEUG00", match, Chem.MolToSmiles(mol), Chem.MolToSmarts(pattObj), coresList, "OV", overlaps, nonEmptyOverlaps)
                if nonEmptyOverlaps:
                    maxOverLaps = max([len(idxset) for idxset in nonEmptyOverlaps])
                    if len(match) > 2 and maxOverLaps < 2:
                        continue
                    if debug:
                        print("--MATCH", match, coresList, "MMM", mayrDct)
                    matched.append((coresList, overlaps, match, mayrDct))
    if not matched:
        if debug:
            print("NOMATCHFOR core", Chem.MolToSmiles(mol), "NUM", len(coresList))
        return None
    bestMatch = matched
    if len(matched) > 1:
        bestMatch = _getBestMatch(mol, matched)
    if debug:
        print("  DEUG match", len(matched), len(bestMatch))
    return bestMatch


def _calcRateFromScores(molNEscores, args):
    if len(molNEscores) != 2:
        print("MOLNE:2245:", len(molNEscores), molNEscores)
        raise NotImplementedError
    molNE01, molNE10 = _getBestMatchingSolventsPair(molNEscores, args)
    if molNE01 and molNE10 is None:
        electroRate = molNE01['electro'][3]['reactivity']
        nucleoRate = molNE01['nucleo'][3]['reactivity']
        sn = molNE01['nucleo'][3]['sN']
    elif molNE01 is None and molNE10:
        electroRate = molNE10['electro'][3]['reactivity']
        nucleoRate = molNE10['nucleo'][3]['reactivity']
        sn = molNE10['nucleo'][3]['sN']
    else:
        if args.mayrMatchMode == 'exact':
            print("SCO", molNEscores)
            raise NotImplementedError
        electroRate10 = molNE10['electro'][3]['reactivity']
        nucleoRate10 = molNE10['nucleo'][3]['reactivity']
        rateNE10 = electroRate10 + nucleoRate10
        electroRate01 = molNE01['electro'][3]['reactivity']
        nucleoRate01 = molNE01['nucleo'][3]['reactivity']
        rateNE01 = electroRate01 + nucleoRate01
        if rateNE01 > rateNE10:
            electroRate = electroRate01
            nucleoRate = nucleoRate01
            sn = molNE01['nucleo'][3]['sN']
        else:
            electroRate = electroRate10
            nucleoRate = nucleoRate10
            sn = molNE10['nucleo'][3]['sN']
    if args.constMayrSN != 0.0:
        sn = args.constMayrSN
    rate = sn * (nucleoRate + electroRate)
    return rate


def _getBestMatchingSolventsPair(molNEscore, args):
    # TODO TO DO
    # here add code to filter molNEscore and return only values from the same solvent or the most similar solvents
    # print(molNEscore)
    nNucleo = 0
    nElectro = 0
    for ne in molNEscore:
        if ne['electro']:
            nElectro += 1
        if ne['nucleo']:
            nNucleo += 1
    if args.verbose > 5 and nNucleo > 0 and nElectro > 0:
        print(f"1. BEST MATCH IN:: {molNEscore}")
    mol0Nsolvents = _getAllSolventsFromMayrMatch(molNEscore[0]['nucleo'])
    mol0Esolvents = _getAllSolventsFromMayrMatch(molNEscore[0]['electro'])
    mol1Nsolvents = _getAllSolventsFromMayrMatch(molNEscore[1]['nucleo'])
    mol1Esolvents = _getAllSolventsFromMayrMatch(molNEscore[1]['electro'])
    mol01NEbestSolvents = _getBestSolventMatch(mol0Nsolvents, mol1Esolvents, args)
    mol01ENbestSolvents = _getBestSolventMatch(mol0Esolvents, mol1Nsolvents, args)
    if args.verbose > 5:
        print(f"2. BEST MATCH SOLVENTS: {mol0Nsolvents} .... {mol1Esolvents}  == {mol01NEbestSolvents}",
              f"111: {mol0Esolvents} ... {mol1Nsolvents} == {mol01ENbestSolvents}")
    # molNE1score = {'nucleo': molNEscore[0]['nucleo'] if not molNEscore[0]['nucleo'] else [molNEscore[0]['nucleo'][0], ],
    #               'electro': molNEscore[0]['electro'] if not molNEscore[0]['electro'] else [molNEscore[0]['electro'][0], ]}
    # molNE2score = {'nucleo': molNEscore[1]['nucleo'] if not molNEscore[1]['nucleo'] else [molNEscore[1]['nucleo'][0], ],
    #               'electro': molNEscore[1]['electro'] if not molNEscore[1]['electro'] else [molNEscore[1]['electro'][0], ]}
    ne1N = _getMayrForGivenSolvent(molNEscore[0]['nucleo'], mol01NEbestSolvents[0])
    ne1E = _getMayrForGivenSolvent(molNEscore[1]['electro'], mol01NEbestSolvents[1])
    ne2N = _getMayrForGivenSolvent(molNEscore[1]['nucleo'], mol01ENbestSolvents[1])
    ne2E = _getMayrForGivenSolvent(molNEscore[0]['electro'], mol01ENbestSolvents[0])
    if ne1N and ne1E:
        molNE1score = {'nucleo': ne1N, 'electro': ne1E}
    else:
        molNE1score = None
    if ne2N and ne2E:
        molNE2score = {'nucleo': ne2N, 'electro': ne2E}
    else:
        molNE2score = None
    if molNE2score is None and molNE1score is None:
        print(f"PROBLEM!!!!  {molNEscore} SOLV: {mol01NEbestSolvents} :: {mol01ENbestSolvents}  :::results::  {ne1N} {ne1E} :: {ne2N} {ne2E}")
        raise
    if args.verbose > 5:
        print(f"3. BEST MATCH IN:: {molNEscore} OUT:: {molNE1score} out22: {molNE2score}")
    return molNE1score, molNE2score


def _getMayrForGivenSolvent(dataList, solvent):
    if not dataList:
        return dataList
    entries = [entry for entry in dataList if entry[-1]['solvent'] == solvent]
    if not entries:
        print("SOLVENTs", solvent, [entry[-1]['solvent'] for entry in dataList])
        return None
        # raise ValueError
    if len(entries) > 1:
        print(f"multiple match to the same solvent {solvent}:: {entries}")
    return entries[0]


def _getAllSolventsFromMayrMatch(matchList):
    if not matchList:
        return None
    solvents = []
    for elem in matchList:
        slv = elem[3]['solvent']
        # if not slv:
        #    print("SLV", slv, elem)
        solvents.append(slv)
    return solvents


def _getBestSolventMatch(solventList1, solventList2, args):
    # when full data for more than one solvent available choose the most common one
    solventPopularity = {'water': 221, 'DMSO': 273, 'MeCN': 203, 'dichloromethane': 16, '': 58,
                         'MeOH-MeCN mix': 29, 'DMF': 13, 'nPrOH-MeCN mix': 2, 'THF': 3,
                         'MeOH': 2, 'EtOH': 2, 'acetone': 2, 'aq acetone': 2, 'aq MeCN': 3,
                         }
    # when data for different solvent available match the most similar pair based on groups
    # the most similar combination is between solvents which has more common groups
    groups = {('iPrOH', 'MeOH', 'nPrOH', 'TFE'): 'alcohols',
              ('water', 'iPrOH', 'MeOH', 'nPrOH', 'TFE', 'aq'): 'protic',
              ('acetone', 'DMSO', 'dichloromethane', 'MeCN'): 'aprotic',
              ('benzene', 'toluene'): 'aromatic',
              ('acetone', 'DMF', 'dichloromethane', 'MeCN', 'DMSO', 'water', 'iPrOH', 'MeOH', 'nPrOH', 'TFE', 'aq'): 'polar',
              ('benzene', 'toluene', 'THF'): 'nonpolar',
              }
    solventTypes = {'polar': ('acetone', 'DMF', 'dichloromethane', 'MeCN', 'DMSO', 'water', 'iPrOH', 'MeOH', 'nPrOH', 'TFE', 'aq'),
                    'nonpolar': ('benzene', 'toluene', 'THF'),
                    'protic': ('water', 'iPrOH', 'MeOH', 'nPrOH', 'TFE', 'aq'),
                    'aprotic': ('acetone', 'DMSO', 'dichloromethane', 'MeCN')}
    solv1 = []
    solv2 = []
    allowedProtonic = args.rxprotonic.split('.')
    allowedPolar = args.rxpolar.split('.')
    if len(allowedPolar) == 1 or len(allowedProtonic) == 1:
        solventList1, solventList2 = removeNotAllowedSolvent(solventList1, solventList2, allowedProtonic, allowedPolar, solventTypes)
    # print("RXINFO", rxinfo['rxpolar'], rxinfo['rxprotonic'], "SOLV", , "SPEC", solventList1, solventList2)
    if not (solventList1 or solventList2):
        return None, None
    if not solventList1:
        return None, list(sorted(solventList2, key=lambda x: solventPopularity.get(x, -1), reverse=True))[0]
    if not solventList2:
        return list(sorted(solventList1, key=lambda x: solventPopularity.get(x, -1), reverse=True))[0], None
    set1 = set(solventList1)
    set2 = set(solventList2)
    intsec = set1.intersection(set2)
    if intsec:
        lst = list(sorted(intsec, key=lambda x: solventPopularity.get(x, -1), reverse=True))[0]
        return lst, lst
    else:
        similarityScore = dict()
        for solv1 in set1:
            for solv2 in set2:
                score = _getSolventSimilScore(solv1, solv2, groups)
                if score not in similarityScore:
                    similarityScore[score] = [(solv1, solv2), ]
                else:
                    similarityScore[score].append((solv1, solv2))
        try:
            maxscore = max(similarityScore.keys())
            solv1, solv2 = similarityScore[maxscore][0]
        except:
            print("!!!!!!!!!!PROBLEM", maxscore, "SSS:", similarityScore)
            raise
        return solv1, solv2


def removeNotAllowedSolvent(solventList1, solventList2, allowedProtonic, allowedPolar, solventTypes):
    # return solventList1, solventList2
    # print("CCXXCC", solventList1, solventList2, allowedProtonic, allowedPolar, solventTypes)
    if solventList1:  # and len(solventList1) > 1:
        allowedSolv1 = [s for s in solventList1 if isSolventAllow(s, allowedProtonic, allowedPolar, solventTypes)]
        if allowedSolv1:
            solventList1 = allowedSolv1
    if solventList2:  # and len(solventList2) > 1:
        allowedSolv2 = [s for s in solventList2 if isSolventAllow(s, allowedProtonic, allowedPolar, solventTypes)]
        if allowedSolv2:
            solventList2 = allowedSolv2
    return solventList1, solventList2


def isSolventAllow(name, allowProtonic, allowPolar, solventTypes):
    if len(allowProtonic) == 1:
        assert allowProtonic[0] in {'P', 'AP'}
        if allowProtonic[0] == 'P':
            if name not in solventTypes['protic']:
                return False
        elif allowProtonic[0] == 'AP':
            if name not in solventTypes['aprotic']:
                return False
    if len(allowPolar) == 1:
        assert allowPolar[0] in {'P', 'Non-P'}
        if allowPolar == 'P':
            if name not in solventTypes['polar']:
                return False
        if allowPolar == 'Non-P':
            if name not in solventTypes['onpolar']:
                return False
    return True


def _getSolventSimilScore(solv1, solv2, groups):
    score = 0
    for gr in groups:
        for elem in gr:
            if elem in solv1 and elem in solv2:
                score += 1
    return score

# #############


def _getBestMatch(mol, matches):
    bestMatchScore = (0, 0)  # overlapCore totalOverlap
    bestMatchList = []
    for coresList, overlaps, match, mayrDct in matches:
        # MXX ((6, 5, 8, 4, 9, 10, 11, 0, 1, 2, 3),) [{0, 1, 4}] (1, 0, 4) MM
        bestPoz = get1stArgMax(overlaps)
        thisMatchScore = (len(overlaps[bestPoz]), len(coresList[bestPoz]))
        if thisMatchScore > bestMatchScore:
            bestMatchScore = thisMatchScore
            bestMatchList = [(coresList[bestPoz], overlaps[bestPoz], match, mayrDct), ]
        elif thisMatchScore == bestMatchScore:
            bestMatchList.append((coresList[bestPoz], overlaps[bestPoz], match, mayrDct))
    coreInMol = [tuple(sorted(matchInfo[1])) for matchInfo in bestMatchList]
    if len(set(coreInMol)) != 1:
        molNoMap = Chem.MolFromSmiles(getUniqSmiles(mol))
        ranks = tuple(Chem.CanonicalRankAtoms(molNoMap, breakTies=False))
        rankOfMatches = [frozenset([ranks[idx] for idx in match]) for match in coreInMol]
        if len(set(rankOfMatches)) == 1:
            coreInMol = [coreInMol[0], ]
        else:
            print("take 1CORE IN MOL", coreInMol, Chem.MolToSmiles(mol), matches, "RANK", rankOfMatches)
            # raise NotImplementedError
            coreInMol = [coreInMol[0], ]  # take first it is not strict solution but seems to be reasonable
    coreInfoList = [(matchInfo[3]['smarts'], matchInfo[3]['smiles']) for matchInfo in bestMatchList]
    opts = {'distance': 'tanimoto', 'radius': 2}
    coreSimilarity = getFragmentBasedSimilarity(mol, coreInMol[0], coreInfoList, opts)
    simPoz = [(sim, poz) for poz, sim in enumerate(coreSimilarity)]
    themaxSim = max(simPoz)[0]
    maxSimPozs = [spoz[1] for spoz in simPoz if spoz[0] == themaxSim]
    debug = False
    if debug:
        print("LLISST", bestMatchList)
        print("MOL", Chem.MolToSmiles(mol))
        print("CORE SIM", coreSimilarity)
        sims = tuple(sorted([(sim, bestMatchList[poz][3]['smilesString']) for poz, sim in enumerate(coreSimilarity)]))
        print(getUniqSmiles(mol), sims)
    return [matches[msp] for msp in maxSimPozs]


def _convertIdxToMapNum(mol, idxes):
    mapNums = []
    for idx in idxes:
        num = mol.GetAtomWithIdx(idx).GetAtomMapNum()
        mapNums.append(num)
    return mapNums


def get1stArgMax(lst, mode='len'):
    if mode != 'len':
        raise NotImplementedError
    maxLen = -1
    maxLenPoz = None
    for poz, elem in enumerate(lst):
        thisLen = len(elem)
        if thisLen > maxLen:
            maxLen = thisLen
            maxLenPoz = poz
    return maxLenPoz


if __name__ == "__main__":
    args = parseArgs()
    mayrData = parseMayrData(args)
    mayrRate = calcMayrRate(mayrData, args)
    print("RATE", mayrRate)
