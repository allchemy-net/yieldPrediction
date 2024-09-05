import time, pickle, argparse
import statistics
from rdkit import Chem
from rdkit.Chem import AllChem


def parseArgs():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smiles", type=str, required=True,  help='input smiles')
    ap.add_argument('--tables', type=str, nargs='+', default=['datapka/pkaTable_radius_1part_1.pickle', 'datapka/pkaTable_radius_1part_7.pickle',
                    'datapka/pkaTable_radius_2part_4.pickle', 'datapka/pkaTable_radius_3part_1.pickle', 'datapka/pkaTable_radius_3part_7.pickle',
                    'datapka/pkaTable_radius_1part_2.pickle', 'datapka/pkaTable_radius_1part_8.pickle', 'datapka/pkaTable_radius_2part_5.pickle',
                    'datapka/pkaTable_radius_3part_2.pickle', 'datapka/pkaTable_radius_3part_8.pickle', 'datapka/pkaTable_radius_1part_3.pickle',
                    'datapka/pkaTable_radius_1part_9.pickle', 'datapka/pkaTable_radius_2part_6.pickle', 'datapka/pkaTable_radius_3part_3.pickle',
                    'datapka/pkaTable_radius_3part_9.pickle', 'datapka/pkaTable_radius_1part_4.pickle', 'datapka/pkaTable_radius_2part_1.pickle',
                    'datapka/pkaTable_radius_2part_7.pickle', 'datapka/pkaTable_radius_3part_4.pickle', 'datapka/pkaTable_radius_1part_5.pickle',
                    'datapka/pkaTable_radius_2part_2.pickle', 'datapka/pkaTable_radius_2part_8.pickle', 'datapka/pkaTable_radius_3part_5.pickle',
                    'datapka/pkaTable_radius_1part_6.pickle', 'datapka/pkaTable_radius_2part_3.pickle', 'datapka/pkaTable_radius_2part_9.pickle',
                    'datapka/pkaTable_radius_3part_6.pickle'], help='file (pickle) with pka values for atomic fingerprint')
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
    t00 = time.time()
    data = dict()
    for fname in args.tables:
        _, rawradius = fname.split('radius_')
        radius = int(rawradius[0])
        t0 = time.time()
        fh = open(fname, 'rb')
        tmp = pickle.load(fh)
        fh.close()
        if radius not in data:
            data[radius] = dict()
        data[radius].update(tmp)
        if args.verbose > 5:
            print(f'tables loaded in {round(time.time() - t0, 2)} s from {fname} deced radius {radius} SIZES: {[(rad, len(data[rad])) for rad in data]}')
    if args.verbose > 1:
        print(f'table loaded in  {round(time.time() - t00, 2)} s, loaded {len(args.tables)} files, radius: {sorted(data.keys())}')
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
