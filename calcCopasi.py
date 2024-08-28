import helpers, consts, calcScaling
import basico


def getAllRxRate(keyToFullRx, target, adbData, commonData, itr, sbsCoef, args):
    yields = dict()
    relFast = None
    # print(activeRx)
    if args.heatFormationBoltzman:
        raise NotImplementedError
    trueFullPoz = []  # when some of normal reaction has yield higher than allowed scale all reaction except this from list
    maxNormalYield = 0  # max normal yield used to determine scaling factor
    for rxpoz, rxkey in enumerate(keyToFullRx):
        rxhash = helpers.getRxHash(keyToFullRx[rxkey])
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
        if rxidsSet.intersection(consts.veryFastRxid):
            yields[rxhash] = args.vfrate
            trueFullPoz.append(rxpoz)
            continue
        hasF = rxidsSet.intersection(adbData['fastRxSet'])
        hasS = rxidsSet.intersection(adbData['slowRxSet'])
        hasVS = rxidsSet.intersection(adbData['verySlowRxSet'])
        hasSpeeds = (hasF, hasS, hasVS)
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
        elif rxidsSet.intersection(consts.sidemechRxid):
            yld = calcScaling.scaleYield(args.sidemechYield, relFast, rxpoz, hasSpeeds, rate, args)  # args.scaleByADBspeed)
        elif isinstance(rxkey, str) and 'REV' in rxkey:
            if args.reverseSpeedMode == 'simple':
                yld = args.reverseSpeed
            elif args.reverseSpeedMode == 'full':
                yld = calcScaling.scaleYield(args.reverseSpeed, relFast, rxpoz, (False, False, False), rate, args)  # args.scaleByADBspeed)
            else:
                raise NotImplementedError
        else:
            sbses = set(keyToFullRx[rxkey]['smiles'].split('>>')[0].split('.'))
            if target in sbses:
                yld = args.fromTargetYield
            else:
                yld = calcScaling.scaleYield(args.defaultYield, relFast, rxpoz, hasSpeeds, rate, args)  # args.scaleByADBspeed)
        # scale by DE
        scallerDE = calcScaling.getDEscaller(rxDE, args)
        yld *= scallerDE
        # scale based on water
        if args.debug > 120:
            print("SVSC", sbsCoef)
        if args.scaleWaterRx != 'no' and 'O' in sbsCoef[rxpoz] and not commonData['mainPathInfo']['isWaterCondNeeded']:
            yld /= args.waterRxScaleDesc
        # scale by conditions
        if args.includeConditions != 'no':
            if args.includeConditions == 'minimal':
                condNum = itr // args.maxiter
                conds = commonData['allowedCond']['minimalVariants'][0][condNum]
                scaler = calcScaling.getScalerMinimConditions(commonData, rxids, conds, adbData, args)
                if condNum in commonData['allowedCond']['minimalPozBeforeChangeToTarget']:
                    scaler *= calcScaling.scaleCompetingWhichAreBeforeCondChange(condNum, keyToFullRx[rxkey]['smiles'], commonData)
                yld *= scaler
            else:
                raise NotImplementedError
        # scale by reversibility
        revRxkey = f'{rxkey}_REV'
        if revRxkey in keyToFullRx:
            yld *= args.forwardReversableScaler
        if yld > maxNormalYield:
            maxNormalYield = yld
        yields[rxhash] = yld
    return yields


######################################################

def extractCopasiResults(copasiRes, target, args):
    mass = {target: copasiRes[target].to_list()}
    fullMass = {cmd: copasiRes[cmd].to_list() for cmd in copasiRes}
    rxChanged = {}
    return mass, fullMass, rxChanged, target


def getRxInCopasiFormat(keyToFullRx, args):
    # return list of dict with stoichiometry coef for each substrate in each reaction
    rxes = []
    for rxpoz, rxkey in enumerate(keyToFullRx):
        if args.debug:
            print("KEY:2075:", rxkey, rxpoz, keyToFullRx[rxkey]['fullRxSmiles'], keyToFullRx[rxkey])
        rxcopasi = keyToFullRx[rxkey]['fullRxSmiles'].replace('>>', ' -> ').replace('.', ' + ')
        rxes.append({'rxstr': rxcopasi, 'rxhash': helpers.getRxHash(keyToFullRx[rxkey])})
    return rxes


def propagateNetworkCopasi(data, sbsStoich, adbData, commonData, args):
    init = [x for x in data['graph'].nodes if data['graph'].in_degree(x) == 0]
    top = [x for x in data['graph'].nodes if data['graph'].out_degree(x) == 0]
    zeroGen = [smi for smi in data['data']['allSmiles'] if data['data']['allSmiles'][smi]['generation'] == 0]
    init = tuple(set(init).union(set(zeroGen)))
    try:
        target = data['data']['targetSmiles']
    except:
        print("DATA", data['data'].keys())
        raise
    if args.debug:
        print("GRAPH", data['graph'].edges)
        print(f"INIT {init} 0gen: {zeroGen}", "\nTOP:", top)
    # mass = {smi: 100 for smi in init}
    m = basico.new_model(name='Simple Model')
    sbsCoef = []
    for rx in getRxInCopasiFormat(data['data']['idToFullRx'], args):
        if args.debug:
            print("rx", rx)
        sbsSmiles = rx['rxstr'].split(' -> ')[0].split(' + ')
        sbsCoef.append({sbs: sbsSmiles.count(sbs) for sbs in set(sbsSmiles)})
        basico.add_reaction(name=rx['rxhash'], scheme=rx['rxstr'])
    if args.debug > 100:
        pars = basico.get_reaction_parameters()
        print("PARS:2149:", type(pars), pars)
    conc = basico.get_species().initial_concentration
    for cmd in conc.index:
        basico.set_species(cmd, initial_concentration=0)
    for cmd in init:
        basico.set_species(cmd, initial_concentration=100)
    conc = basico.get_species().initial_concentration
    #
    if args.includeConditions == 'minimal':
        if args.debug > 20:
            print(f"LEN CON {len(commonData['allowedCond']['minimalVariants'])} :0: {commonData['allowedCond']['minimalVariants'][0]}")
        ratesList = []
        numSteps = []
        pathlen = len(commonData['allowedCond']['minimalVariants'][0])
        prevNumStep = 0
        # print("KKK", commonData['allowedCond'].keys())
        if len(commonData['workingCond']) > 1:
            # raise NotImplementedError
            print("WARNING", commonData, commonData['workingCond'])
        for poz, step in enumerate(commonData['workingCond'][0]):
            itr = step[0] * args.maxiter
            if args.debug:
                print("ITR", step[0])
            rates = getRxConst(data, target, adbData, commonData, itr, sbsCoef, args)
            ratesList.append(rates)
            if args.debug:
                print("RATES", rates)
            if poz > 0:
                thisNumStep = step[0]
                prevNumStep = commonData['workingCond'][0][poz - 1][0]
                numSteps.append((thisNumStep - prevNumStep) * args.maxiter)
        thisNumStep = commonData['workingCond'][0][-1][0]
        # print("THIS", pathlen, thisNumStep)
        numSteps.append((pathlen - 1 - thisNumStep) * args.maxiter)
    else:
        raise NotImplementedError
    # print("CONC", type(conc), conc)
    # print("RATES", ratesList, "\nNUMSTEPS:", numSteps)
    for poz, numStep in enumerate(numSteps):
        rxRates = ratesList[poz]
        for rxhash in rxRates:
            rxrate = rxRates[rxhash]
            # print("RXH", rxhash, rxrate, rxhash in basico.get_reactions().index)
            basico.set_reaction_parameters(name=f'({rxhash}).k1', value=rxrate)
        if args.debug > 120:
            print("+++RATES", basico.get_reactions())
        result = basico.run_time_course(numStep, update_model=True)
        if args.debug > 120:
            print("RES", result[target])
    return result, target


def getRxConst(data, target, adbData, commonData, itr, sbsCoef, args):
    # rx -> sbs
    if args.debug > 10:
        print("DEBUG::174 begin propagation step")
    rxYields = getAllRxRate(data['data']['idToFullRx'], target, adbData, commonData, itr, sbsCoef, args)
    # sbsUsage = getIzolatedSubstrateUsage(activeRx, keyToFullRx, sbsCoef, rxYields, mass)
    # sbsUsage = getNormalizedSubstrateUsage(sbsUsage, sbsCoef, mass, args)
    if itr == 0 and 'tautomers' in data['data']:
        tautoData = data['data']['tautomers']['rxnum']
        for tauto in tautoData:
            sbs, prod, rxid, rxpoz, whichToCalcRatio = tautoData[tauto]
            rxhash = helpers.getRxHash(data['data']['idToFullRx'][rxpoz])
            revrxhash, revrxstr = helpers.getRevRxHashAndStr(data['data']['idToFullRx'][rxpoz], mode='tautomer')
            print(f"\n\n {sbs}>> {prod}", rxhash, rxhash in basico.get_reactions().index)
            # basico.add_reaction(name=revrxhash, scheme=revrxstr)
            xProd = args.setTautomerRatio / (1 + args.setTautomerRatio)
            if whichToCalcRatio == 0:
                xProd = 1 - xProd
            prodSbsRatio = xProd / (1 - xProd)
            if revrxhash not in basico.get_reactions().index:
                basico.add_reaction(name=revrxhash, scheme=revrxstr)
            if prodSbsRatio > 1:  # product dominates increase speed of forward rx
                forwardrate = prodSbsRatio * args.eqrate
                basico.set_reaction_parameters(name=f'({rxhash}).k1', value=forwardrate)
                basico.set_reaction_parameters(name=f'({revrxhash}).k1', value=args.eqrate)
            else:
                revrate = (1 / prodSbsRatio) * args.eqrate
                basico.set_reaction_parameters(name=f'({revrxhash}).k1', value=revrate)
                basico.set_reaction_parameters(name=f'({rxhash}).k1', value=args.eqrate)
    if 'acidBase' in data['data']:
        conds = None
        if args.includeConditions != 'no':
            if args.includeConditions == 'minimal':
                condNum = itr // args.maxiter
                conds = commonData['allowedCond']['minimalVariants'][0][condNum]
            else:
                raise NotImplementedError
        if itr == 0:
            _addAcidBaseRxCopasi(data['data']['acidBase']['wildForms'], commonData['acidBase'], args)
        if conds:
            _changeAcidBaseCopasi(data['data']['acidBase']['wildForms'], commonData['acidBase'], conds, args)
    else:
        print("DATAK", data['data'].keys(), 'acidBase' in data['data'])

    #    mass, changedRx = _setAcidBaseRatio(mass, changedRx, data['acidBase'], commonData, itr, args)
    # get acid base and tautomer fast rate
    return rxYields


def _addAcidBaseRxCopasi(dataWildForms, dataAcidBase, args):
    # dataDict[name][gen]['data']['acidBase'] = {'wildForms': thisWilds}  # , 'xOfMainForm': thisMains}
    # objects in dict should link to the same wildToMain mainToWild and noMainForm
    # dataDict[name][0]['common']['acidBase'] = {'noMainForm': noMainForm, 'wildToMain': wildToMain, 'mainToWild': mainToWild, 'wildFormType': wildFormType}
    for smilesWild in dataWildForms:
        mainForm = dataAcidBase['wildToMain'][smilesWild]
        formType = dataAcidBase['wildFormType'][smilesWild]
        rxhashMain, rxstrMain, rxhashRev, rxstrRev = makeABrxesAndHashes(smilesWild, mainForm, formType)
        xOfWildForm = args.setAcidBaseRatio / (1 + args.setAcidBaseRatio)
        xOfMainForm = 1 - xOfWildForm
        if xOfWildForm > xOfMainForm:
            raise NotImplementedError
        else:
            rateratio = xOfMainForm / xOfWildForm
            if rxhashMain not in basico.get_reactions().index:
                print("RXHASmain", rxhashMain)
                basico.add_reaction(name=rxhashMain, scheme=rxstrMain)
                basico.set_reaction_parameters(name=f'({rxhashMain}).k1', value=args.eqrate * rateratio)
            else:
                print(f"second try of adding same {rxhashMain} {basico.get_reaction_mapping(rxhashMain)} RATE: {args.eqrate * rateratio}")
                raise NotImplementedError
            if rxhashRev not in basico.get_reactions().index:
                basico.add_reaction(name=rxhashRev, scheme=rxstrRev)
                basico.set_reaction_parameters(name=f'({rxhashRev}).k1', value=args.eqrate)
            else:
                print(f"second try of adding same {rxhashRev} {basico.get_reaction_mapping(rxhashRev)} RATE: {args.eqrate}")
                print("second try of adding same rx rev")
                raise NotImplementedError


def _changeAcidBaseCopasi(dataWildForms, dataAcidBase, conds, args):
    print("COND", conds)
    condB = 'SB' in conds[0] or 'B' in conds[0]
    condN = 'W' in conds[0] or 'N' in conds[0]
    condA = 'SA' in conds[0] or 'A' in conds[0]
    cond = 'N'
    if condN:
        cond = 'N'
    elif condB:
        cond = 'B'
    elif condA:
        cond = 'A'
    #
    xOfWildForm = args.setAcidBaseRatio / (1 + args.setAcidBaseRatio)
    xOfMainForm = 1 - xOfWildForm
    if xOfWildForm > xOfMainForm:
        raise NotImplementedError
    else:
        rateratio = xOfMainForm / xOfWildForm
    rate = args.eqrate * rateratio
    for rxhash in basico.get_reactions().index:
        if rxhash.endswith('::A_towild'):
            if cond == 'B':
                basico.set_reaction_parameters(name=f'({rxhash}).k1', value=0)
            else:
                basico.set_reaction_parameters(name=f'({rxhash}).k1', value=rate)
        elif rxhash.endswith('::B_towild'):
            if cond == 'A':
                basico.set_reaction_parameters(name=f'({rxhash}).k1', value=0)
            else:
                basico.set_reaction_parameters(name=f'({rxhash}).k1', value=rate)


def makeABrxesAndHashes(smilesWild, mainForm, formType):
    if formType not in {'A', 'B', 'M'}:
        print(f"SMILES WILD {smilesWild} MAIN: {mainForm} FORM: {formType}")
        raise
    mainhash = f'{smilesWild}>>{mainForm}::MAIN'
    mainstr = f'{smilesWild} -> {mainForm}'
    revhash = f'{mainForm}>>{smilesWild}::{formType}_towild'
    revstr = f'{mainForm} -> {smilesWild}'
    return mainhash, mainstr, revhash, revstr
