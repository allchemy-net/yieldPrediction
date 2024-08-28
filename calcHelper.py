def formCalcInfo(data, args):
    infoDict = dict()
    if 'rxAvail' in data['data']:
        infoDict['cmdToDisableRx'] = data['data']['rxAvail']['cmd'].copy()
        infoDict['cmdToDisableMax'] = {cmd: 0 for cmd in data['data']['rxAvail']['cmd']}
        infoDict['disableRx'] = set()
        for cmd in data['data']['rxAvail'].get('cmd', []):
            infoDict['disableRx'].update(data['data']['rxAvail']['cmd'][cmd])
    if args.debug > 10:
        print(f"infoDict {infoDict}")
    return infoDict


def updateCalcInfo(infoDict, masses, args):
    cmdToRm = []
    for cmd in infoDict['cmdToDisableRx']:
        # print("CMD", cmd, len(masses))
        thismass = masses[-1].get(cmd, 0)
        maxmass = infoDict['cmdToDisableMax'][cmd]
        if thismass > 20:
            cmdToRm.append(cmd)
            if args.debug > 0:
                print(f"OLD {infoDict['disableRx']}")
            infoDict['disableRx'] -= infoDict['cmdToDisableRx'][cmd]
            if args.debug > 0:
                print(f"DISABLE CX {cmd} {len(masses)} {thismass} {maxmass} {infoDict['disableRx']}")
            continue
        if len(masses) > 2000:
            cmdToRm.append(cmd)
            if args.debug > 0:
                print(f"OLD {infoDict['disableRx']}")
            infoDict['disableRx'] -= infoDict['cmdToDisableRx'][cmd]
            if args.debug > 0:
                print(f"DISABLE CX {cmd} {len(masses)} {thismass} {maxmass} {infoDict['disableRx']}")
            continue
        if thismass > maxmass:
            infoDict['cmdToDisableMax'][cmd] = thismass
            continue
        if maxmass > 1 and thismass / maxmass < 0.99:
            cmdToRm.append(cmd)
            if args.debug > 0:
                print(f"OLD {infoDict['disableRx']}")
            infoDict['disableRx'] -= infoDict['cmdToDisableRx'][cmd]
            if args.debug > 0:
                print(f"DISABLE CX {cmd} {len(masses)} {thismass} {maxmass} {infoDict['disableRx']}")

    for cmd in cmdToRm:
        del infoDict['cmdToDisableRx'][cmd]
    return infoDict