import sys, re, argparse
from functools import cache
from rdkit import Chem
import consts


@cache
def getMainForm(smiles):
    # this function is for cache only which should speed up
    return smilesHelper.getMainForm(smiles)


def getProductOfRx(rxn, mol, maxrx=5):
    prevProd = None
    isPrevSelective = False
    for itr in range(maxrx):
        prods = rxn.RunReactants([mol, ])
        if not prods:
            if not isPrevSelective:
                raise NotImplementedError
            smiles = prevProd
            break
        prods = [p[0] for p in prods]
        uniqSmiles = [getUniqSmiles(p) for p in prods]
        isSelective = len(set(uniqSmiles)) == 1
        smiles = Chem.MolToSmiles(prods[0])
        mol = Chem.MolFromSmiles(smiles)  # kind of sanitization of mol
        if isSelective:
            break
        #if not smiles:
        #    raise NotImplementedError
        prevProd = smiles
        isPrevSelective = isSelective
    else:
        raise NotImplementedError
    return smiles


@cache
def getMainTautomerForm(smiles):
    mol = Chem.MolFromSmiles(smiles)
    for patt in consts.tautomerisationDict:
        if mol.HasSubstructMatch(patt):
            rxn = consts.tautomerisationDict[patt]
            smiles = getProductOfRx(rxn, mol)
            mol = Chem.MolFromSmiles(smiles)
    return smiles


def getMainTautoToSideTautosDict(allCmds):
    mainToSide = dict()
    for cmd in allCmds:
        mainCmd = getMainTautomerForm(cmd)
        if mainCmd == cmd:
            continue
        if mainCmd not in mainToSide:
            mainToSide[mainCmd] = set()
        mainToSide[mainCmd].add(cmd)
    return mainToSide


def getMapNumFromSmiles(smiles):
    nomapsmiles = re.findall(':\d+]', smiles)
    numMaps = [int(rawnum[1:-1]) for rawnum in nomapsmiles]
    return numMaps


def demapSmiles(smiles):
    nomapsmiles = re.sub(':\d+]', ']', smiles)
    if '[1' in nomapsmiles:
        mol = Chem.MolFromSmiles(nomapsmiles)
        for atm in mol.GetAtoms():
            atm.SetIsotope(0)
        return Chem.MolToSmiles(mol)
    return Chem.CanonSmiles(nomapsmiles)


def getUniqSmiles(mol):
    smi = Chem.MolToSmiles(mol)
    nonums = re.sub(':[0-9]*]', ']', smi)
    smi2 = Chem.CanonSmiles(nonums)
    return smi2


def getRxHash(rxinfo):
    rxhash = rxinfo['fullRxSmiles'] + '::' + ':'.join([str(rxid) for rxid in rxinfo['rxid']])
    return rxhash


def getRevRxHashAndStr(rxinfo, mode):
    sbs, prod = rxinfo['fullRxSmiles'].split('>>')
    revrx = '>>'.join([prod, sbs])
    if mode == 'tautomer':
        modeCode = 'T'
    elif mode == 'acidbase':
        modeCode = 'AB'
    else:
        raise NotImplementedError
    rxhash = revrx + f':{modeCode}:' + ':'.join([str(rxid) for rxid in rxinfo['rxid']])
    rxstr = revrx.replace('>>', ' -> ').replace('.', ' + ')
    return rxhash, rxstr


def replaceArgFromDict(args, kwargs):
    newargs = argparse.Namespace(**vars(args))
    multiVals = dict()
    for name in kwargs:
        if name.endswith('_MULTIPART'):
            namemod, pozs = name.split('__PART')
            pozs = pozs.split('_')
            thisPoz = int(pozs[0])
            if namemod not in multiVals:
                thisLen = int(pozs[1])
                multiVals[namemod] = [None for _ in range(thisLen)]
            multiVals[namemod][thisPoz] = kwargs[name]
        else:
            newargs.__dict__[name] = kwargs[name]
    for name in multiVals:
        print("COMPLEX", name, multiVals[name])
        newargs.__dict__[name] = multiVals[name]
    return newargs


#####


def _getMatchingMolFromList(mol, mappedMolList):
    matched = []
    for idx, mapped in enumerate(mappedMolList):
        if mol.HasSubstructMatch(mapped) and mapped.HasSubstructMatch(mol):
            matched.append((idx, mapped))
    return matched


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
