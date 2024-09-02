## Information about helper scripts
In this directory there are deposited scipts for:
- calculating Mayr rate for given reaction (mayr.py)
- predicting and rank the most acidic C-H site in molecule (pkarank.py)
- converting json with reaction pathway to pickle which can be used as input to yieldPrediction.py


## Details about scripts:
- ### pkarank.py
The scipt use precalculated pka values of atoms in given topological enviroment to predict most acidic sites in the molecule. 
The pka values were obtained from our [pKa predictor](https://pubs.acs.org/doi/full/10.1021/jacs.9b05895). The program works as follow:
  - for each carbon atom in input molecule calculates fingerprint centered on the atom. By default series of fingerprints with radius 1, 2 and 3 are calculated.
  - search for matching fingerprint in provided database, starting from fingeprints with the largest radius, and return found pKa

User can provides his/her own database with other radius and values from other source(s). 

Usage: `python ./pkarank.py --smiles SMILES`, for example: `python ./pkarank.py --smiles 'CCC(=O)OCC'`


- ### mayr.py
The scipt use Mayr table from dataMayr subdirecotry. The data was collected from [Mayr's group webpage](https://www.cup.lmu.de/oc/mayr/reaktionsdatenbank/).
For details of algoritm see Figure 4 in the paper.

Usage:  `python mayr.py --rxsmiles RXSMILES --rxsmarts RXSMARTS` 

where RXSMILES is SMILES of reaction and RXSMARTS is SMARTS of reaction template,

for example: `python mayr.py --rxsmiles 'O=C=O.[O-]c1ccc(I)cc1>>O=C([O-])Oc1ccc(I)cc1' --rxsmarts '[O:1]=[C:2]=[O:3].[O-:4][c:5]>>[c:5][O+0:4][*:2](=[O:3])[O-:4]'`

See `python mayr.py --help` for optional parameters.

- ### converJSONtoPickle.py
The script use json with reaction pathway (see format description below) and generates pickle. Basic usage `python ./convertJSONtoPickle.py -i PATHWAY.json --adb ../inputdata/reactionInfo.json -o PATHWAY.pickle` 
where `PATHWAY.json` is user-generated file with information about pathway,  `../inputdata/reactionInfo.json` is file with information about reaction and `PATHWAY.pickle`
