## Information about helper scripts
In this directory there are deposited scipts for:
- calculating Mayr rate for given reaction (mayr.py)
- predicting and rank the most acidic C-H site in molecule (pkarank.py)
- converting json with reaction pathway to pickle which can be used as input to yieldPrediction.py


## Details about scripts:
### pkarank.py
The scipt use precalculated pka values of atoms in given topological enviroment to predict most acidic sites in the molecule. 
The pka values were obtained from our [pKa predictor](https://pubs.acs.org/doi/full/10.1021/jacs.9b05895). The program works as follow:
  - for each carbon atom in input molecule calculates fingerprint centered on the atom. By default series of fingerprints with radius 1, 2 and 3 are calculated.
  - search for matching fingerprint in provided database, starting from fingeprints with the largest radius, and return found pKa

User can provides his/her own database with other radius and values from other source(s). 

Usage: `python ./pkarank.py --smiles SMILES`, for example: `python ./pkarank.py --smiles 'CCC(=O)OCC'`

The command needs to load database (table) with pkas values hence depending on your computer it may takes up to few minues.

### mayr.py
The scipt use Mayr table from dataMayr subdirecotry. The data was collected from [Mayr's group webpage](https://www.cup.lmu.de/oc/mayr/reaktionsdatenbank/).
For details of algoritm see Figure 4 in the paper.

Usage:  `python mayr.py --rxsmiles RXSMILES --rxsmarts RXSMARTS` 

where RXSMILES is SMILES of reaction and RXSMARTS is SMARTS of reaction template,

for example: `python mayr.py --rxsmiles 'O=C=O.[O-]c1ccc(I)cc1>>O=C([O-])Oc1ccc(I)cc1' --rxsmarts '[O:1]=[C:2]=[O:3].[O-:4][c:5]>>[c:5][O+0:4][*:2](=[O:3])[O-:4]'`

See `python mayr.py --help` for optional parameters.

### converJSONtoPickle.py
The script use json with reaction pathway (see format description below) and generates pickle. 

Basic usage `python ./convertJSONtoPickle.py -i PATHWAY.json --adb ../inputdata/reactionInfo.json -o PATHWAY.pickle` 

where `PATHWAY.json` is user-generated file with information about pathway,  `../inputdata/reactionInfo.json` is file with information about reaction and `PATHWAY.pickle`
Input json file need to have following keys (see file [exampleconverter/rx31_.C2_F3_I1_O2_c6.allgens.json](exampleconverter/rx31_.C2_F3_I1_O2_c6.allgens.json) for example):
  - graph - list of edges in reaction pathway. Nodes in graph are compounds (smiles string) and reaction number (int)
  - rxes - dict where keys are reaction numbers (as defined in 'graph') whereas value is dict with following keys:
    * smiles - reaction smiles (without duplicated substrates and sideproducts)
    * rxid - reaction number as defined in reactionInfo.json
    * DE - ΔE of reaction 
    * byproducts - byproducts from reactants (e.g. when butyllithium is used for deprotonation most likely byproducts is butane)
    * mappedRX - full reaction (with all side-prodcuts and duplicated substrates) with mapping and proper stoichiometry
    * fullRxSmiles - full reaction as above but without atom mapping
    * rate - (optional) Mayr rate calculated by mayr.py script
  - rxesL3 - as rxes but for L3 level of analysis
  - graphL3 - as graph but for L3 level of analysis
  - rxesL4 - as rxes but for L4 level of analysis
  - graphL4 - as graph but for L4 level of analysis
  - targetSmiles - smiles of final product of pathway
  - catalysts - list of catalysts
  - isOnepot - true/false defined is reaction is one pot or MCR
  - highlyReactive - list of highly reactive species
  - initC - (optional) if given initial concentration of given substrates will be overwriten with value given here
  - formAB - dict with protonated/deprotonated compounds present in reaction pathway. Key is protonated/deprotonate smiles and value is list `[formType, mainForm]`
    where formType is 'acid' or 'base' and mainForm is smiles after quenching/neutralisation.
