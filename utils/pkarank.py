import time, pickle, argparse
import statistics
from rdkit import Chem
from rdkit.Chem import AllChem


def parseArgs():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smiles", type=str, required=True,  help='input smiles')
    ap.add_argument('--tables', type=str, default='datapka/pkaAtomicTables.pickle', help='file (pickle) with pka values for atomic fingerprint')
    ap.add_argument('--verbose', type=int, default=0, help='verbosity level')
    ap.add_argument('--radius', type=int, default=3, help='maximal radius used to matching. Should be the same as in data stored in --tables')
    return ap.parse_args()


def getFSfromCH(args):
    mol = Chem.MolFromSmiles(args.smiles)
    bitInfo = dict()
    _ = AllChem.GetMorganFingerprint(mol, args.radius, bitInfo=bitInfo)
    # {787341104: ((2, 1),), 847957139: ((3, 0),), 1173125914: ((1, 1),), 2245384272: ((1, 0), (2, 0)), 2246728737: ((0, 0),), 2542228502: ((1, 2),), 2592785365: ((3, 1),), 3542456614: ((0, 1),)}
    CHs = set([atm.GetIdx() for atm in mol.GetAtoms() if atm.GetTotalNumHs() > 0 and atm.GetSymbol() == 'C'])
    atomFPs = dict()
    for fp in bitInfo:
        for atmFPinfo in bitInfo[fp]:
            atmIdx, fpRad = atmFPinfo
            if atmIdx not in CHs:
                continue
            if fpRad == 0:
                continue
            if atmIdx not in atomFPs:
                atomFPs[atmIdx] = dict()
            atomFPs[atmIdx][fpRad] = fp
    return atomFPs


def loadTables(args):
    t0 = time.time()
    fh = open(args.tables, 'rb')
    data = pickle.load(fh)
    fh.close()
    if args.verbose > 2:
        print(f'tables loaded in {round(time.time() - t0, 2)} s')
    return data


def matchPkas(atomFPs, tables, args):
    atomPkas = dict()
    for atom in atomFPs:
        fps = atomFPs[atom]
        maxRadius = max(fps)
        radius = maxRadius
        while radius > 0:
            if tables[radius] and fps[radius] in tables[radius]:
                pkas = tables[radius][fps[radius]]
                if args.verbose > 1:
                    print(f"atom {atom} radius {radius}: {len(pkas)} entries in range: {min(pkas)} - {max(pkas)}")
                atomPkas[atom] = (pkas, radius)
                break
            radius -= 1
    return atomPkas


def rankPka(pkaData):
    rank = dict()
    for atom in pkaData:
        pkas, radius = pkaData[atom]
        if len(pkas) > 1:
            if radius == 1:
                pka = max(pkas)
            else:
                pka = statistics.mean(pkas)
        else:
            pka = pkas[0]
        rank[atom] = pka
    return sorted(rank, key=lambda x: rank[x]), rank


def main(args):
    atomFPs = getFSfromCH(args)
    if args.verbose > 1:
        print('atomic fingerprint for C-H: ', atomFPs)
    tables = loadTables(args)
    pkas = matchPkas(atomFPs, tables, args)
    ranks, values = rankPka(pkas)
    return ranks, values


if __name__ == "__main__":
    args = parseArgs()
    ranks, values = main(args)
    print(f'position sorted from most acidit: {ranks};  pka Values: {values}')
