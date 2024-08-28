## Introdution
This repository provides code and data for reaction yield prediction. 
The model is described in the paper "Estimation of multicomponent reactions’ yields from networks of mechanistic steps" and is also used and mentioned in "Systematic, computational discovery of multicomponent reactions and one-pot sequences." by Rafał Roszak, Louis Gadina, Maciej Krzeszewski, Agnieszka Wołos, Ahmad Makkawi, Barbara Mikulak-Klucznik, Olena Vakuliuk, Yasemin Bilgi, Karol Molga, Patrycja Gołębiowska, Mariusz Tasior, Oskar Popik, Tomasz Klucznik, Sara Szymkuć, Martyna Moskal, Sebastian Baś, Mohammad B. Teimouri, Rafał Frydrych, Jacek Mlynarski, Daniel T. Gryko & Bartosz A. Grzybowski

## Deposited data
All data are stored in inputdata directory and consists of:
- modelParams.json - optimized model paramethers toghether with program paramethers using during training and evaluation
- reactionInfo.json - information about all reactions calculated used in all deposited reaction paths
- reaction paths, there are three types of path:
  - training set (deposited in subdirectory trainingSet), described and used in "Estimation of multicomponent reactions’ yields from networks of mechanistic steps" paper
  - test set (deposited in subdirectory testSet) described and used in both papers
  - additional reactions (deposited in subdirectory additalReactions) presented in Figure 6 in "Systematic, computational discovery of multicomponent reactions and one-pot sequences."

For details about deposited paths see [READMEpaths.md](inputdata/READMEpaths.md)


## Code
### Overview
Code for prediction yield from path are stored in main directory. The code requires as a input:
- reaction databased (databased needed for calculating all deposited paths is stored in inputdata/reactionInfo.json file)
- model parameters (trained model is stored in inputdata/modelParams.json)
- information about path(s) which be calculated (example pickle files are strored in inputdata/testSet, inputdata/trainingSet and inputdata/additionalReactions)

To run the code use yieldPrediction.py script (see example usage below), for full list of options run `python yieldPrediction.py --help` however easier way to change options is to edit namespace file.
### Requirements:
- Python 3.x (tested on Python 3.9 and 3.11)
- networkX 2.x (tested on 2.5)
- cachetools (tested on 5.3.1 and 5.5.0)
- numpy (tested on 1.23.1 and 2.0.0)
  
### Example usage:
Calculate yield for all deposited paths with precalculated model parameters and original reaction data base:
 `python ./ yieldPrediction.py --namespace inputdata/modelParams.json  --loadrxdb inputdata/reactionInfo.json --loadpickles inputdata/*/*.pickle`
