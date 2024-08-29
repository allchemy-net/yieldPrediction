## Introduction
This repository provides code and data for reaction yield prediction. 
The model is described in the paper "Estimation of multicomponent reactions’ yields from networks of mechanistic steps" and is also used and mentioned in "Systematic, computational discovery of multicomponent reactions and one-pot sequences." by Rafał Roszak, Louis Gadina, Maciej Krzeszewski, Agnieszka Wołos, Ahmad Makkawi, Barbara Mikulak-Klucznik, Olena Vakuliuk, Yasemin Bilgi, Karol Molga, Patrycja Gołębiowska, Mariusz Tasior, Oskar Popik, Tomasz Klucznik, Sara Szymkuć, Martyna Moskal, Sebastian Baś, Mohammad B. Teimouri, Rafał Frydrych, Jacek Mlynarski, Daniel T. Gryko & Bartosz A. Grzybowski

## Deposited data
All data are stored in inputdata directory and consists of:
- modelParams.json - optimized model parameters together with program parameters used during training and evaluation
- reactionInfo.json - additional information on reactions used in all deposited reaction paths like: general reaction conditions (strongly basic, basic, mildly basic, neutral, mildly acidic, strongly acidic, and whether the reaction requires Lewis acid), class of solvent (protic-aprotic and polar-nonpolar), temperature range (VL: < -20 °C, L: -20 to 20 °C, RT, H: 40 to 150 °C, VH: > 150 °C), water tolerance (yes, no, water is necessary), typical reaction speeds (very slow, slow, fast, very fast, and unknown if inconsistent data have been reported in the literature) and more.
- reaction paths - there are three types of path:
  - training set (deposited in the subdirectory trainingSet), described and used in "Estimation of multicomponent reactions’ yields from networks of mechanistic steps" paper
  - test set (deposited in the subdirectory testSet) described and used in both papers
  - additional pathways (deposited in subdirectory additalReactions) presented in Figure 6 of the paper "Systematic, computational discovery of multicomponent reactions and one-pot sequences."

For details about deposited paths see [READMEpaths.md](inputdata/READMEpaths.md)


## Code
### Overview
The code for predicting the yield from a given mechanistic path is stored in the main directory. The code requires as input:
- reaction database (database needed for calculating all deposited paths is stored in inputdata/reactionInfo.json file)
- model parameters (trained model is stored in inputdata/modelParams.json)
- information about path(s) to be calculated (example pickle files are strored in inputdata/testSet, inputdata/trainingSet and inputdata/additionalReactions)

To run the code use yieldPrediction.py script (see example usage below). For a full list of options run `python yieldPrediction.py --help`. However an easier way to change options is to edit the namespace file.
### Requirements:
- Python >= 3.9 (tested on Python 3.9, 3.10 and 3.11)
- networkX 2.x (tested on 2.5)
- cachetools (tested on 5.3.1 and 5.5.0)
- numpy (tested on 1.23.1 and 2.0.0)
  
### Example usage:
Calculate yield for all deposited paths with precalculated model parameters and original reaction data base:
 `python ./yieldPrediction.py --namespace inputdata/modelParams.json  --loadrxdb inputdata/reactionInfo.json --loadpickles inputdata/*/*.pickle`
