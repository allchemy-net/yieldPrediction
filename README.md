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

In utils/ directory are available helper scripts for generating input files. See [readme in utils/](utils/README.md) for details

### Requirements:
- Python >= 3.9 (tested on Python 3.9, 3.10 and 3.11)
- networkX 2.x (tested on 2.5)
- cachetools (tested on 5.3.1 and 5.5.0)
- rdkit (tested on 2022.03)
- numpy (tested on 1.23.1 and 2.0.0)
- scikit-learn (tested on 0.24.2)
The easiest way to install dependencies is [using conda](https://www.rdkit.org/docs/Install.html)

### Instalation:
No instalation is required to use the software. Please install required library and download/clone this repository. 

### Example usage:
Calculate yield for all deposited paths with precalculated model parameters and original reaction data base:
 `python ./yieldPrediction.py --namespace inputdata/modelParams.json  --loadrxdb inputdata/reactionInfo.json --loadpickles inputdata/*/*.pickle`

Abovementioned command will takes several minutes to complite on moder laptop/desktop.

### Output format:
```
$ python ./yieldPrediction.py --namespace inputdata/modelParams.json  --loadrxdb inputdata/reactionInfo.json --loadpickles inputdata/*/*.pickle
copasi mode not available
calculation running with following configs/options:  Namespace(input=None, longCutFile='/code/modules/yieldpred/mech_shortcuts.tsv', heatFormation='/code/modules/yieldpred/policzoneMopaciem', reactiveSmarts='/code/modules/yieldpred/reactive.smarts', adb=None, MayrFilesExp=['/code/modules/yieldpred/Nucleophiles.tsv', '/code/modules/yieldpred/Electrophiles.tsv'], MayrFilesPred='/code/modules/yieldpred/mayr_params_from_luosz', sbsStoichFile=None, masterconfig='', namespace='', mutexRxFile='/code/modules/yieldpred/mutexrx.tsv', baseUsage='/code/modules/yieldpred/baseinfo.tsv', debug=0, pprint=True, exportGraph=None, saveDataForPlot='', archiveResult=False, paropt='N', scan='', optSteps=27500, openboxSurrogate='gp', disableOpt='', applyHistory='', fitTo='12', transferLearning='', epsilon=1e-05, numprocs=64, vfrate=1000000.0, eqrate=10000.0, disableLongCut=1, includeAdditives=False, temperature=300, multiplyCat=50, additivesScale=0.0001, removeLateSides='rm1', addReverse='YnoTAB', filterReverse='arom', stoichiometry='optimal', removeDuplicates='sbsesProd', overwriteStoich='', propagationMode='all', maxiter=600, includeConditions='minimal', onepotfile='/code/modules/yieldpred/onepot.txt', catDef='fixPd', catAddTime='peak', includepka='binary1', minpkaprob=0.8, rmRxMode='N', scaleByADBspeed=[2.0859169241895286, 0.37654380677490384, 0.010040453467482027], scaleSideInOldCondition=0.9990131459471931, scaleWrongConditionRx=0.0223767532762191, reverseSpeed=0.5596350930308015, reverseSpeedMode='simple', forwardReversableScaler=2.793282597102021, fromTargetYield=0.0007391942260074859, fromTargetDef='sideWithTautoDeProt', sidemechYield=0.4494087444808005, defaultYield=0.899894020761921, mainFormYield=0.0012094162963565945, maxConsumption=0.20037937481947318, maxStepYield=0.95, scaleWaterRx='whenWaterSbs', waterRxScaleDesc=299.94937543796044, stoichScale='c100', c100=50, redox='rmRedOxWater', shortcutspeed=0.9970772628774253, heatFormationBoltzman=False, scaleByDE=0.9986455110073102, mayrMatchMode='core', useMayr='exp', scaleByMayr=0.8684788056458217, constMayrSN=1.0, maxAllowDE=3, maxAllowDEtr=1, useEvansPolanyi='Y', maxAllowMayrNE=30, mayrAndDEMode='expTexp', boostIntramol=1.0011591951741015, mutexRxMode='sidemech', fixTautomerRatio='constEnolate', tautomerProductIs='lessStable', setTautomerRatio=0.5901492904849115, setTautoEnolateRatio=1, fixAcidBaseRatio='condIgnoreWeak', setAcidBaseRatio=1.4877108665539922, yieldDefMode='maxWithTautoAndDeProt', loadpickles=['inputdata/additionalReactions/mol13_.B3_C16_N2_c40_n9_r3.allgens.pickle', 'inputdata/additionalReactions/mol14_.C1_N1_c35_n1.allgens.pickle', 'inputdata/additionalReactions/mol15_.B1_C6_O1_c25_n1_r1.allgens.pickle', 'inputdata/additionalReactions/mol16_.C10_N1_c24_n2.allgens.pickle', 'inputdata/testSet/mach10_.C8_N3_O2.allgens.pickle', 'inputdata/testSet/mach1_.C13_O1_S1_c2_i1_n3.allgens.pickle', 'inputdata/testSet/mach2_.C18_O1_S1_c8_i1_n3.allgens.pickle', 'inputdata/testSet/mach3_.C9_O1_c6.allgens.pickle', 'inputdata/testSet/mach4_.C10_N1_O4_c16.allgens.pickle', 'inputdata/testSet/mach5_.C10_F3_O2_S1_c12.allgens.pickle', 'inputdata/testSet/mach6_.C10_N1_O1_S1_c12.pathid_1.allgens.pickle', 'inputdata/testSet/mach7_.C12_O5_c12.pathid_0.allgens.pickle', 'inputdata/testSet/mach8_.C14_O2.allgens.pickle', 'inputdata/testSet/mach9_.C10_O1_c6.allgens.pickle', 'inputdata/trainingSet/rx12_.C10_N1_O3_c3_n1_s1.allgens.pickle', 'inputdata/trainingSet/rx13_.C5_O4_c11_n4.allgens.pickle', 'inputdata/trainingSet/rx14_.C13_O5.pathid_1.allgens.pickle', 'inputdata/trainingSet/rx15_.C11_N2_O2_c21_n1.pathid_1.allgens.pickle', 'inputdata/trainingSet/rx16_.C10_N2_O2_c16_n2.allgens.pickle', 'inputdata/trainingSet/rx17_.C3_c23_n5.allgens.pickle', 'inputdata/trainingSet/rx18_.C5_H1_N1_O3_P1_c18.allgens.pickle', 'inputdata/trainingSet/rx19_.C7_N2_O2_c10_o1.allgens.pickle', 'inputdata/trainingSet/rx20_.B1_C13_N1_c12_r1.allgens.pickle', 'inputdata/trainingSet/rx21_.C15_N1_O6_c3_n1_s1.pathid_1.allgens.pickle', 'inputdata/trainingSet/rx22_.C5_N3_O1_c6.pathid_0.allgens.pickle', 'inputdata/trainingSet/rx23_.C7_O1_c6_l3.allgens.pickle', 'inputdata/trainingSet/rx24_.C14_N4_O1_c6.allgens.pickle', 'inputdata/trainingSet/rx25_.C5_H1_N1_c14_n1.allgens.pickle', 'inputdata/trainingSet/rx26_.C8_N1_O4.allgens.pickle', 'inputdata/trainingSet/rx27_.C10_N1_O1_c12.allgens.pickle', 'inputdata/trainingSet/rx28_.C8_N2_O3_c12.pathid_0.allgens.pickle', 'inputdata/trainingSet/rx29_.C14_O1.allgens.pickle', 'inputdata/trainingSet/rx30_.C10_F3_O5_c6.allgens.pickle', 'inputdata/trainingSet/rx31_.C2_F3_I1_O2_c6.allgens.pickle'], loadrxdb='inputdata/reactionInfo.json')
64 process(es) will be started
RESULTS
input name	L2-level yield	L3-level yield	L4-level yield
mach1    	 46.3    	 53.8    	 53.0
mach2    	 47.5    	 28.8    	 28.3
mach3    	 74.4    	 68.3    	 67.6
mach4    	 98.3    	 98.9    	 98.9
rx18    	 99.5    	 99.5    	 99.5
mol16    	 3.5    	 2.2    	 4.3
mach7    	 53.3    	 53.8    	 77.9
mach10    	 84.3    	 47.0    	 50.5
rx20    	 12.2    	 95.5    	 95.5
rx19    	 70.7    	 60.3    	 47.6
mach8    	 49.8    	 48.1    	 29.7
rx14    	 49.0    	 49.0    	 49.0
mach6    	 36.6    	 11.6    	 0.7
rx21    	 26.6    	 20.6    	 19.2
rx13    	 49.0    	 36.9    	 35.8
rx25    	 28.4    	 11.8    	 4.8
rx23    	 99.5    	 27.0    	 26.9
rx26    	 49.3    	 41.5    	 38.2
rx15    	 65.5    	 55.9    	 55.7
mach5    	 95.6    	 39.6    	 29.2
rx16    	 73.8    	 76.4    	 76.8
mach9    	 99.5    	 35.1    	 8.1
rx24    	 45.0    	 20.5    	 27.7
mol14    	 26.5    	 25.3    	 22.4
rx12    	 76.7    	 77.2    	 76.6
rx31    	 99.9    	 99.9    	 99.9
rx17    	 41.4    	 71.9    	 75.3
rx27    	 98.1    	 78.7    	 77.3
rx22    	 98.7    	 46.3    	 63.8
rx29    	 97.6    	 99.8    	 99.8
rx28    	 64.4    	 55.0    	 55.0
rx30    	 85.6    	 62.7    	 15.0
mol15    	 18.3    	 2.3    	 4.2
mol13    	 17.2    	 15.4    	 6.6
```
Program print predicted yield (in %) for L2, L3 and L4 level of analysis. Input name is begining of file name up to fiest underscore (`_`) character.
