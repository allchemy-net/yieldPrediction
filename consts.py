from rdkit import Chem
from rdkit.Chem import AllChem


enolateC = Chem.MolFromSmarts('[CX3-][CX3](=[OX1])[#6,O]')
enolateO = Chem.MolFromSmarts('[CX3]=[CX3]([O-])[#6,O]')
enoleToKetone = AllChem.ReactionFromSmarts('[CX3:1]=[CX3:2]([OH:3])[#6,O:4]>>[C:1][C:2](=[OX1:3])[*:4]')
tautomerisationStrList = [('[CX3]=[CX3]([OH])[#6,O]', '[CX3:1]=[CX3:2]([OH:3])[#6,O:4]>>[C:1][C:2](=[OX1:3])[*:4]'), ]
tautomerisationDict = {Chem.MolFromSmarts(lst[0]): AllChem.ReactionFromSmarts(lst[1]) for lst in tautomerisationStrList}

expYields = {"rx29": 100, "rx18": 98, "rx27": 93, "rx20": 89, "rx31": 86, "rx28": 52, "rx15": 78, "rx16": 74, "rx17": 72, "rx12": 69,
             "rx19": 65, "rx30": 62, "rx13": 58, "rx14": 54, "rx22": 53, "rx23": 48, "rx25": 44, "rx21": 54, "rx26": 37, "rx24": 33,
             "mach2": 22, "mach1": 34, "mach4": 97, "mach7": 57, "mach9": 35, "mach8": 31, "mach10": 62, "mach3": 82, "mach5": 57, "mach6": 14}

condABtoInt = {'SA': -3, 'A': -2 , 'WA': -1, 'N': 0, 'WB': 1, 'B': 2, 'SB': 3}
condIntToAB = {condABtoInt[ab]: ab for ab in condABtoInt}
condTempToInt = {'VL': -2, 'L': -1, 'L+RT': -0.5,  'RT': 0, 'H': 1, 'VH': 2}
condIntToTemp = {condTempToInt[temp]: temp for temp in condTempToInt}


class WrongData(Exception):
    pass


class WrongOptions(Exception):
    # when combination of used cmd arguments has no sense
    pass


class NegativeMassException(Exception):
    # when after step calculation one of mass is negative
    pass
