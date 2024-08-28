import math
# from numba import jit


def getDEscaller(rxDE, args):
    if args.mayrAndDEMode == 'exp':
        scale = math.exp(-rxDE * args.scaleByDE)
        return scale
    elif args.mayrAndDEMode == 'expsigmoid':
        x = rxDE * args.scaleByDE
        scale = 2 / (1 + math.exp(x))
        return scale
    elif args.mayrAndDEMode == 'expKsigmoid':
        x = rxDE * args.scaleByDE
        scale = 2 / (1 + math.exp(rxDE))
        return scale
    elif args.mayrAndDEMode == 'expLELU':
        if rxDE < 0:
            scale = -rxDE + 1
        else:
            scale = 1 / abs(-rxDE - 1)
        return scale
    elif args.mayrAndDEMode == 'expLOLI':
        # (x < 0 ? 1-x : 1/log(x+2.718)
        if rxDE < 0:
            scale = - rxDE + 1
        else:
            scale = 1 / math.log(rxDE + math.e)
        return scale
    elif args.mayrAndDEMode == 'expTexp':
        rxDEtr = rxDE
        if args.maxAllowDEtr < abs(rxDE):
            if rxDE < 0:
                rxDEtr = - args.maxAllowDEtr
            else:
                rxDEtr = args.maxAllowDEtr
        scale = math.exp(-rxDEtr)
        return scale
    elif args.mayrAndDEMode == 'expabs' or args.mayrAndDEMode == 'linear':
        mode = 'newLinear'
        if args.mayrAndDEMode == 'expabs':
            mode = 'oldStableLinear'
        if rxDE > 0:
            scale = rxDE * args.scaleByDE
            return (1 + scale)
        elif rxDE < 0:
            if mode == 'oldStableLinear':
                scale = rxDE * args.scaleByDE
                return (1 - scale)
            elif mode == 'newLinear':
                # mode linear new
                scale = abs(rxDE * args.scaleByDE)
                return 1 - min(0.99, scale)
            else:
                raise NotImplementedError
        else:  # exact 0
            return 1
    else:
        raise NotImplementedError


def getMayrScaller(rate, args):
    if args.mayrAndDEMode == 'exp':
        rate = math.exp(-rate)
        return rate
    elif args.mayrAndDEMode == 'linear':
        if rate > 0:
            rate = 1 + (rate * args.scaleByMayr)
        elif rate < 0:
            rate = 1 / (abs(rate) * args.scaleByMayr)
        else:
            rate = 0
        if rate < 0:
            print("RATE NEGATIVE", rate)
            raise
        return rate


def getBoostScaler(rxinfo, args):
    if 'boostRing' not in rxinfo:
        scaler = 1
    else:
        scaler = rxinfo['boostRing'] * args.boostIntramol
    if rxinfo['isShortcut']:
        scaler *= args.shortcutspeed
    return scaler


def getStoichScaler(rxinfo, masses, sbsCoef, args):
    if args.stoichScale == 'c100':
        if len(sbsCoef) == 1:
            sbs = list(sbsCoef.keys())[0]
            return sbsCoef[sbs]
        massCoef = [masses[sbs]/sbsCoef[sbs] for sbs in sbsCoef]
        massCoef.sort()
        if args.debug > 10:
            print("c100: coef scaling:", massCoef)
        if args.debug > 150:
            print(f"c100:99: MASS: {masses} COEF: {sbsCoef}")
        return massCoef[1] / args.c100
    elif 'stoichScale' not in rxinfo:
        return 1
    scaler = 1 / rxinfo['stoichScale']
    return scaler


# @jit(nopython=True)
def scaleYield(yld, relFast, rxpoz, hasSpeed, rate, args):  # speedScale):
    # speedScale = args.scaleByADBspeed
    if relFast:
        yld = yld * relFast[rxpoz]
    if any(hasSpeed):
        if hasSpeed[0]:   # fast
            yld = yld * args.scaleByADBspeed[0]
        if hasSpeed[1]:   # slow
            yld = yld * args.scaleByADBspeed[1]
        if hasSpeed[2]:   # very slow
            yld = yld * args.scaleByADBspeed[2]
    if rate:
        # y = yld
        yld *= rate
    # scale by DE
    # if args.scaleByDE:
    #    yld
    return yld
