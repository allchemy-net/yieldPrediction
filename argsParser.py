import argparse
import consts


def parseArgs():
    args, rangeDict = _parseArgs()
    printWarnings(args)
    # check correctness of arguments
    if args.useMayr != 'N':
        if args.scaleByMayr == 0:
            raise consts.WrongOptions
        if args.useMayr == 'exp':
            if args.mayrMatchMode == 'exact':
                raise consts.WrongOptions
        elif args.useMayr == 'luosz':
            if args.mayrMatchMode != 'exact':
                raise consts.WrongOptions
        else:
            raise NotImplementedError
        # calcRxRate(allSmiles, mayrData)
    if args.c100 <= 0:
        raise consts.WrongOptions

    # disableOpt
    if not args.disableOpt:
        pass
    elif args.disableOpt == 'noBoostIntramol':
        del rangeDict['boostIntramol']
        args.boostIntramol = 1.0
    elif args.disableOpt == 'noWater':
        args.waterRxScaleDesc = 1.0
        del rangeDict['waterRxScaleDesc']
    elif args.disableOpt == 'onlyMayr':
        del rangeDict['boostIntramol']
        args.boostIntramol = 1.0
        args.waterRxScaleDesc = 1.0
        del rangeDict['waterRxScaleDesc']
        args.scaleByADBspeed = [1.0, 1.0, 1.0]
        del rangeDict['scaleByADBspeed']
        del rangeDict['sidemechYield']
        del rangeDict['defaultYield']
        args.defaultYield = 1.0
        args.sidemechYield = 1.0
        args.scaleSideInOldCondition = 1.0
        del rangeDict['scaleSideInOldCondition']
        args.useEvansPolanyi = 'N'
        del rangeDict['maxAllowMayrNE']
    elif args.disableOpt == 'noMayr':
        args.useMayr = 'N'
        del rangeDict['scaleByMayr']
        args.scaleByMayr = 1.0
    elif args.disableOpt == 'noSpeedClass':
        args.scaleByADBspeed = [1.0, 1.0, 1.0]
        del rangeDict['scaleByADBspeed']
    elif args.disableOpt == 'noSideMech':
        del rangeDict['sidemechYield']
        del rangeDict['defaultYield']
        args.defaultYield = 1.0
        args.sidemechYield = 1.0
    elif args.disableOpt == 'noSideOldCondScale':
        args.scaleSideInOldCondition = 1.0
        del rangeDict['scaleSideInOldCondition']
    elif args.disableOpt == 'noRevScal':
        del rangeDict['reverseSpeed']
        args.reverseSpeed = 1.0
        del rangeDict['forwardReversableScaler']
        args.forwardReversableScaler = 1.0
    elif args.disableOpt == 'noDE':
        args.useEvansPolanyi = 'N'
        del rangeDict['maxAllowDE'] 
        del rangeDict['maxAllowDEtr']
    elif args.disableOpt == 'noDEnoSideOldCond':
        args.scaleSideInOldCondition = 1.0
        del rangeDict['scaleSideInOldCondition']
        args.useEvansPolanyi = 'N'
        del rangeDict['maxAllowDE']
        del rangeDict['maxAllowDEtr']
    elif args.disableOpt == 'no2c':
        args.stoichScale = 'uniq'
    else:
        print("ARGS", type(args.disableOpt), args.disableOpt, "::")
        raise NotImplementedError
    return args, rangeDict


def _parseArgs():
    rangeDict = dict()
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # #########################################################################
    # options which involved external input files
    parser.add_argument('--input', type=str, nargs='+', help='file format is {PREFIX}.calcedGen_{GEN}.pickle '
                        'where prefix is any string (but same for one calculations) and gen is int starts with 0')
    parser.add_argument('--loadpickles', type=str, nargs='+', help='file name needs to starts with {PREFIX}_')
    parser.add_argument('--longCutFile', type=str, default='mech_shortcuts.tsv')
    parser.add_argument('--heatFormation', default='policzoneMopaciem', type=str, help='file with calculated with mopac PM6 heat of formation')
    parser.add_argument('--reactiveSmarts', type=str, default='reactive.smarts', help='list of smarts of reactive groups (used by --removeLateSide)')
    parser.add_argument('--adb', default='BazaReakcji23.11.24_ADB', type=str)
    parser.add_argument('--loadrxdb', type=str, help='load reaction database from JSON file')
    parser.add_argument('--MayrFilesExp', type=str, nargs=2, default=['Nucleophiles.tsv', 'Electrophiles.tsv'],
                        help="files with Mayr experimental params. Sara's format collected from Mayr webpage")
    parser.add_argument('--MayrFilesPred', type=str, default='mayr_params_from_luosz', help='files with calculated Mayr params')
    parser.add_argument('--sbsStoichFile', type=str, help="file with stoichiometry of substrates in Aga's format; remember to define --stoichiometry")
    parser.add_argument('--masterconfig', type=str, default='', help='when set configuration is readed from the file and rest from cmd line.'
                        '*NOTE* file has higher priority than cmd args; expected format JSON')
    parser.add_argument('--namespace', type=str, default='', help='full options readed from namespace only following option can be readed from cmd: '
                        '--input, --saveDataForPlot, --debug --numprocs --exportGraph')
    parser.add_argument('--mutexRxFile', type=str, default='mutexrx.tsv', help='file with `mutual exclusion` id est when one is in main path '
                        'other is unlikely to proceed. Exact behaviour is set by --mutexRxMode option')
    parser.add_argument('--baseUsage', type=str, default='baseinfo.tsv', help='file with number of base equivalents')

    # #########################################################################
    # options which change output information/format
    parser.add_argument('--debug', type=int, default=0)
    parser.add_argument('--pprint', action='store_true')
    parser.add_argument('--exportGraph', type=str)
    parser.add_argument('--saveDataForPlot', type=str, default='', help='full data in each step will be stored in file under given prefix')
    parser.add_argument('--archiveResult', action='store_true',
                        help='generate output filename based on full args and place there input args and obtained results')

    # ##########################################################################
    # parameter optimization
    parser.add_argument('--paropt', default='N', type=str, choices=['R2', 'R2_0', 'N', 'MAE', 'MSE', 'RMSE', 'pearson'],
                        help='params optimization: N - dont perform opt; R2 - minimize 1-R2 for linear model y=ax+b; '
                        'MAE - minimize MAE (mean absolute error); MSE - minimize MSE (mean squared error)'
                        'R2_0 - minimize 1-R2 for linear model y=ax; pearson - minimize 1 - r_pearson; ')
    parser.add_argument('--scan', default='', type=str, help='perform param scan based on given file;. File format json like in rangeDict')
    parser.add_argument('--optSteps', default=850, type=int, help='number of params optimization steps')
    parser.add_argument('--openboxSurrogate', default='gp', type=str, choices=['gp', 'prf', 'auto'],
                        help='surrogate type for openbox. Available options: '
                        ' gp - gaussian process - default, works well but become slow after few hunders points'
                        ' prf - probabilistic random forest - need additional lib - better scalability but need many points to get nice results'
                        ' auto - gaussian process and after 300 points switch to prf (surrogate model  auto in openbox)')
    parser.add_argument('--disableOpt', default='', type=str, choices=['', 'noMayr', 'onlyMayr', 'noBoostIntramol', 'noSpeedClass',
                          'noSideMech',  'noRevScal', 'noSideOldCondScale', 'noWater', 'noDE', 'noDEnoSideOldCond', 'no2c'],
                        help='Disable some features of model (this is mainly for statistic purpose). Available option: '
                        ' onlyMayr - disable all correction from paper and based on Mayr'
                        ' noMayr - disable Mayr correction'
                        ' noBoostIntramol - disable boost for intramolecular ring formation (correction 6. in paper)'
                        ' noSpeedClass - disable different speed rate for VS, S, F reaction (correction 1. in paper)'
                        ' noSideMech - dont differentiate between normal and sidemech reactions (correction 2. in paper) '
                        ' noWater - dont scale down side-reaction with water when water not required on main path (correction 3. in paper)'
                        ' noSideOldCondScale - dont scale side-reaction which perform in old conditions, during change of condition in '
                        '    main path (correction 4. in paper)'
                        ' noRevScal - dont add special treatment for reversible reactions (correction 5. in paper)'
                        ' noDE - dont include DE'
                        ' no2c - '
                        ' noDEnoSideOldCond - combined noDE and noSideOldCondScale')
    parser.add_argument('--applyHistory', default='', type=str, help='continue parameters optimization using history from given file')
    parser.add_argument('--fitTo', default='12', choices=['0', '1', '2', '12', '012', '01', '02'], help='Define to which level of calculation perform fiting.'
                        ' 0 - fit to L2 (no additional calculation); 1 - fit to L3 (1g of additional calculation); 2 - fit to L4 (2g of '
                        ' additional calculation). Any combination of 0,1,2 are allowed, e.g. 12 (default) mean fit to L3 and L4')
    parser.add_argument('--transferLearning', default='', type=str, nargs='*', help='read data from previous optimization and make transfer learning '
                        'based on such data')

    # #########################################################################
    # minor technical options (no effect on chemistry)
    parser.add_argument('--epsilon', type=float, default=1e-5)
    parser.add_argument('--numprocs', type=int, default=64)
    parser.add_argument('--vfrate', type=float, default=1e6,
                        help='rate of VF reaction should be much higher that normal rx (neeeded only for --propagationMode copasi)')
    parser.add_argument('--eqrate', type=float, default=1e4, help='base rate of fixed equilibrium (tautomer and acid base equilibrium) '
                        ' this make sense only for --propagationMode copasi and should be much higher than normal rx')

    # #########################################################################
    # initial conditions (concentration, connectivity in the graph)
    parser.add_argument('--disableLongCut', type=int, default=1, help='0 = do not disable longcut, 1 = remove longcut reactions;'
                        'longcut is opposite to shortcut; i.e. sequence of steps which is also present in graph as one step')
    parser.add_argument('--includeAdditives', action='store_true')
    parser.add_argument('--temperature', default=300, type=int, help='when --heatFormation is enable use given temperature for Boltzman distribution')
    parser.add_argument('--multiplyCat', default=50, type=int, help='when substrate stoichiometry is read from file '
                        'multiply catalysts stoiciometry by given number')
    parser.add_argument('--additivesScale', type=float, default=0.0001, help='set concentration of additives to: sbs_conentration * this factor')
    parser.add_argument('--removeLateSides', type=str, default='rm1', choices=['N', 'rm1', 'rm2'],
                        help='remove late side reactions with reactive intermediates/substrates, Options: '
                        'N - dont remove late side; rm1 - remove reaction where reactive intermediate is used in generation 1+ compare to main path '
                        ' e.g. in main path compounds is used in 3rd generation so any reaction with the compound in 4th and further generation will removed; '
                        'rm2 - remove reaction used in generation 2+ compare to main path (see example above)')
    parser.add_argument('--addReverse', type=str, default='YnoTAB', choices=['N', 'Y', 'wildToMainOnly', 'DE', 'YnoTAB'],
                        help='add reverse reaction. Allow values: N - dont add; '
                        'Y - add for all reversible reaction (as deined in allchemy DB); '
                        'DE - add for all reversible with speed based on energy diff; '
                        'YnoTAB - add for reversible but not for (de)protonation and tautomerisation make sense when --fixTautomerRatio and '
                        ' and --fixAcidBaseRatio is set'
                        'wildToMainOnly - add reversible only for reaction which produce wild form, speed of reverse is define by --reverseSpeed')
    parser.add_argument('--filterReverse', type=str, default='arom', choices=['N', 'arom'], help='Remove selected type of reversible connections. Options:'
                        ' N - dont remove any reverse connections (so reverse connection is fully controlled by --addReverse option)'
                        ' arom - remove reversible connection when forward reaction involve aromatization')
    parser.add_argument('--stoichiometry', default='optimal', type=str, choices=['equal', 'optimal', 'fromFile'], help='initial concentration of substrates: '
                        'equal - all have the same values;  optimal - initial concentration proportional to stoichiometry/usage in main path; '
                        'fromFile - use values defined in --sbsStoichFile')
    # rangeDict['addReverse'] = {'values': ('N', 'Y', 'DE'), 'default': 'Y'}  # dont use wildToMainOnly as incompatible with acidBase equlibrium
    parser.add_argument('--removeDuplicates', default='sbsesProd', type=str, choices=['N', 'sbsesProd'],
                        help='control if/when add different reaction from/to the same cmds. N - no ignoring add all rxes; '
                        'sbsesProd - add only one rx with the same all substrates and product')
    parser.add_argument('--overwriteStoich', default='', type=str, help='overwrite stoichiometry for given graph and compound(s). Accepted values:'
                        ' empty string or filename in such format: graphname<tab>canonical_smiles<tab>concntration')

    # #########################################################################
    # main mode of network propagation
    parser.add_argument('--propagationMode', choices=['front', 'all', 'copasi'], default='all',
                        help='front - means propagate network only to new nodes and stop when all reached. Final concentrations are reported '
                        'all - propagates all possible reaction in each step, list of values at each step is reported '
                        'copasi - use basico library (based on copasi) to solv equetion (NEW) and much faster')
    parser.add_argument('--maxiter', type=int, default=500, help="maximal iteration, when -includeConditions != no this is # of iteration per conditions")
    rangeDict['maxiter'] = {'min': 10, 'max': 600, 'step': 10, 'default': 390}
    parser.add_argument('--includeConditions', type=str, choices=['no', 'minimal', 'minimalSide', 'minimalMild'], default='minimal',
                        help='Define if reaction conditions should be included in network propagation. Available options: '
                        'no - all reaction ca proceed any time regardless its conditions'
                        'minimal - detect minimal conditions changes range and all not matching according to --scaleWrongConditionRx options.'
                        '  note in this mode range of conditions can be allowed'
                        'minimalMild - as minimal but when set of conditions are allowed select this one with is closer to next/prev conditions '
                        '   and in second step (if needed) those with is closer to neural and root temperature'
                        'minimaSide - as minimal but when more than one conditions is allowed at particular step '
                        '   select from matched conditions in such a way to minimize possible side reactions. '
                        '   This also can product more than one allowed conditions per step')
    parser.add_argument('--onepotfile', type=str, default='', help='By default all reactions networks are treated as multicomponent reaction id est '
                        ' all substrates are available from the begining (except catalysts which is treated as defined by --catDef and --catAddTime), '
                        ' this option allow to give file with name of result which should be treated as onepot, id est substrates are added just in time')
    parser.add_argument('--catDef', default='Pd', type=str, choices=['Pd', 'fixPd'], help='definition of catalyst for --catAddTime option. Allowed options: '
                        'Pd - only smiles [Pd] is treat as catalyst')
    parser.add_argument('--catAddTime', default='init', type=str, choices=['init', 'jitEarly', 'jitLate', 'peak'],
                        help='define when catalyst is added to reaction mixture. Aloved values: '
                        ' init - available from begining as all other substrates '
                        ' jitEarly - available at stage (minimal generation) when needed;  make sense only with --propagationMode=front, '
                        '   e.g. when --maxiter is 100 and catalyst is used in 3rd generation reaction with catalysts will be allowed only '
                        '  from step 300 onward, when diverse nonbalanced branches (e.g. one branch has 2 rx other 3 rx) reaction is available '
                        '   in earlier generation and onward. In other words: minimal generation counting from substrate when catalysts is used'
                        ' jitLate - as above but in nonbalanced branches available at later generation. I.e. minimal distance from target/final product. '
                        '   This value of the option also affect condition order when --includeConditions=minimal in such a way that conditions for '
                        '   branch with catalyst will be later (i.e. first conditions for branch without catalysts and then for this with catalyst) '
                        ' peak - add catalyst just after (1st) peak of concentration of needed intermediate')
    parser.add_argument('--includepka', type=str, default='binary1', choices=['N', 'binary', 'binary1', 'binary0'],
                        help='how to treat deprotonation; N - do nothing, '
                        'binary - disable reactions with other deprotonation number, e.g. when 2 eq of base is used dianion is formed and '
                        ' only 2nd most acidic position reacts other reaction involing deprotonation will be removed from consideration, '
                        ' when no info about number of eq given dont remove reactions at all; '
                        'binary1 - like above but when no information given allow only reaction on most acidic side; '
                        'binary0 - like above but when no information given dont allow deprotonation at all (usefull only for tests);')
    parser.add_argument('--minpkaprob', type=float, default=0.8, help='minimal pka probability for 1st position to ban reaction on other possision '
                        'or to ban this one if 2 eq of base used, make sense only when --includepka != N')
    parser.add_argument('--rmRxMode', type=str, default='N', choices=['N', 'mainAwild', 'mainAwildALA'], help='Remove less probable reactions: '
                        ' N - dont remove (disable); mainAwild - remove reaction between main form and its anionic wild form '
                        ' mainAwildALA - as mainAwild plus additional removal sidereaction with LA when A or SA on main path and vice versa '
                        '  (i.e. no SA/A in side when LA on main path)')

    # ########################################################################
    # rate/scale based on allchemy db
    parser.add_argument('--scaleByADBspeed', type=float, nargs=3, default=[1.0, 1.0, 1.0], help='scale speed for F, S and VS reactions')
    # 1 correction: r_rate_class
    rangeDict['scaleByADBspeed'] = {'min': (1.25, 0.1, 0.01), 'max': (2.5, 0.5, 0.1), 'step': (0.01, 0.01, 0.01), 'default': (1.3, 0.4, 0.1)}
    parser.add_argument('--scaleSideInOldCondition', type=float, default=1.0, help='scale side rx which are in "old" condition. '
                        'E.g. rx1->rx2, rx1->side1 when rx1 and side1 are in cond1 and rx2 in cond2. Without scaling mass easily go to side1 and when '
                        'conditions change to cond2 no substrate left for rx2. Scaling slow down side1 which faciitate rx2')
    rangeDict['scaleSideInOldCondition'] = {'min': 0.00001, 'max': 1, 'default': 0.0001}
    parser.add_argument('--scaleWrongConditionRx', type=float, default=1e-5,
                        help='when --includeConditions != no multiply wrong condition by given factor, should be in range 0 (no rx at all) up to 1')
    # 4 correction: r_cond
    rangeDict['scaleWrongConditionRx'] = {'min': 0.000001, 'max': 1, 'default': 0.001}
    parser.add_argument('--reverseSpeed', type=float, default=1.5, help='speed of reversible reaction, it works only when --addReverse is properly set')
    # 5 correction: r_rev (1st part  of 2)
    rangeDict['reverseSpeed'] = {'min': 0.001, 'max': 2.0, 'default': 0.3}
    parser.add_argument('--reverseSpeedMode', type=str, default='simple', choices=['simple', 'full'],
                        help='simple means no Mayr scalling for reverse rx')
    parser.add_argument('--forwardReversableScaler', type=float, default=1.0, help='scale forward reactions which are reversible')
    # 5. correction: r_rev (2nd part of 2)
    rangeDict['forwardReversableScaler'] = {'min': 0.001, 'max': 2.8, 'default': 1.0}  # was 5.0
    parser.add_argument('--fromTargetYield', type=float, default=0.01, help='yield should positive but below 1.0')
    rangeDict['fromTargetYield'] = {'min': 0.000001, 'max': 0.75, 'default': 0.0001}
    parser.add_argument('--fromTargetDef', type=str, default='sideWithTautoDeProt', choices=['strict', 'sideWithTautoDeProt'],
                        help='Define which reactions are scaled by --fromTargetYield. Available options: '
                        ' "strict" - only one exact smiles the target is treated as target'
                        ' "sideWithTautoDeProt" - scale all side reaction(s) (i.e. not those from main path) which involve target or '
                        '   any of (de)protonated or tautomeric form of target')
    parser.add_argument('--sidemechYield', type=float, default=0.1, help='yield should positive but below 1.0')
    # 2 correction: : r_class
    rangeDict['sidemechYield'] = {'min': 0.00001, 'max': 1, 'default': 0.01}
    parser.add_argument('--defaultYield', type=float, default=0.75, help='yield should positive but below 1.0')
    rangeDict['defaultYield'] = {'min': 0.003, 'max': 0.9, 'default': 0.1}
    parser.add_argument('--mainFormYield', type=float, default=1.0, help='yield should positive but below 1.0')  # 0.1
    rangeDict['mainFormYield'] = {'min': 0.001, 'max': 0.1, 'default': 0.01}
    parser.add_argument('--maxConsumption', type=float, default=0.9, help='yield should positive but below 1.0. '   # TODO
                        'Max consumption of substrate, it is used for initial scalling when compound is involved in 2+ reactions')
    rangeDict['maxConsumption'] = {'min': 0.2, 'max': 0.95, 'default': 0.9}

    parser.add_argument('--maxStepYield', type=float, default=0.95, help='maximal yield of step allowed value range (0.0, 1.0). '
                        'This is apply at the end after all other scalling')
    parser.add_argument('--scaleWaterRx', type=str, default='whenWaterSbs', choices=['no', 'whenWaterSbs', 'whenSbsNoRev'],
                        help='Scale reaction based on water needed. Allowed Options: "no" - dont scale at all; '
                        '"whenWaterSbs" - scale reactions which use water as substrate when water is not needed in main path; '
                        '"whenSbsNoRev" - scale competing reactions when water not involved in mainPath and not in rev. rx')
    parser.add_argument('--waterRxScaleDesc', type=float, default=1.0, help='divide yield/rate by this const')
    # 3. correction: r_water
    rangeDict['waterRxScaleDesc'] = {'min': 100.0, 'max': 300.0, 'default': 150.0}
    # rangeDict['scaleWaterRx'] = {'values': }
    parser.add_argument('--stoichScale', type=str, default='c100', choices=['N', 'simple2', 'simple', 'uniq', 'c100'],
                        help='Scale reaction yield/rate based on stoichiometry: N - no scalling;'
                        ' simple - rate(yield) / number_of_components; simple2 - rate(yield) / (number_of_compounds - 1)'
                        ' uniq - rate(yield) / number_of_uniq_components'
                        ' c100 - rate(yield) * conc_sbs1/100 * conc_sbs2/100 * ... multiply by concentration nonlimiting substrates')
    parser.add_argument('--c100', type=float, default=50, help='for --stoichScale c100 scale 2c rx when second (more abundant) reactants is below this number')
    parser.add_argument('--redox', type=str, default='rmRedOxWater', choices=['N', 'NAVS', 'NAsidemech', 'rm', 'rmRedOxWater'],
                        # gdy main: 1) NA -> side/competing ox -> sidemech  2) NA -> side/competing red -> sidemech
                        help='Use red/ox information to slow down or remove unlikely reactions, ie. when on main path'
                        ' there is no oxidative (reductive) conditions (i.e. NA) then slow down (remove) side/competing reactions'
                        ' which requires reductive (oxidative) conditions. Note when Ox/Red is on main path then Red/Ox are removed'
                        ' from alowed conditions in allchemy search results. When one use data from other sources preprocessing'
                        ' which remove wrong Ox/Red competing reactions is needed. Available options: '
                        ' N - dont use this correction, i.e. dont remove/slow down any Ox/Red reactions '
                        'NAVS - make noncompatible side/competing reactions very slow '
                        'NAsidemech - make noncompatible reaction sidemech'
                        'rm - completely remove from graph all ox and red reaction when NA on path'
                        'rmRedOxWater - remove nicompatible Red reaction whereas incompatible Ox reaction scale down'
                        ' using value of --waterRxScaleDesc')
    parser.add_argument('--shortcutspeed', type=float, default=1, help='scale reaction defined as shortcut by this number')
    rangeDict['shortcutspeed'] = {'min': 0.1, 'max': 1.0, 'default': 0.9, 'step': 0.1}

    # ########################################################################
    # rate/scale based on additionl/extended information not from allchemy db
    parser.add_argument('--heatFormationBoltzman', action='store_true', help='calc rate of competing reaction using Boltzman equation '
                        'This is suitable for thermodynamic equlibrion here it rather lowering correlation so dont switch it on')
    parser.add_argument('--scaleByDE', type=float, default=0.1, help='scale constant for reaction with DE lower than --maxAllowDE '
                        'for DE > 0 (endothermic): rate *= 1 - (DE * scale) and for DE < 0: rate *= 1 + (DE * scale) '
                        'scale should be possitive (reasonable range: 0.001-1); when 0, scalling is not used')
    rangeDict['scaleByDE'] = {'min': 0.000001, 'max': 1, 'default': 0.001}
    parser.add_argument('--mayrMatchMode', type=str, choices=['exact', 'core', 'coreAnion'], default='core', help='matching mode: '
                        'exact - canon smiles (string) comparison;  core - find most similar among those with the same reacting core'
                        'coreAnion - when no matching for nucleophile and rxid is in mayrNeutralToAnion perform deprotonation in rx core and matching')
    parser.add_argument('--useMayr', type=str, choices=['N', 'exp', 'luosz'], default='exp', help='source of Mayr paramethers: '
                        'N - dont include Mayr; exp - use experimental data (low coverage) luosz - use values predicted by Luo 10.1002/cphc.202300162')
    parser.add_argument('--scaleByMayr', type=float, default=2.8, help='scale reaction based on Mayr params N + E '
                        ' scaling for N+E > 0:  scaling=(N+E)*scaleByMayr then yield *= scaling'
                        ' for N+E < 0: scaling= 1/ ((N+E)*scaling) then yield *= scaling; when scaling or N+E is 0.0 no scalling is apply'
                        ' N+E is range ca. (-40, 40) but usually in (-10, 10) so resonable value is range ca. (0.01, 1). Lower values is ')
    rangeDict['scaleByMayr'] = {'min': 0.0001, 'max': 2, 'default': 0.87}
    parser.add_argument('--constMayrSN', type=float, default=1.0, help='when != 0.0 use const value instead of real Sn paramether')
    parser.add_argument('--maxAllowDE', default=10, type=int, help='calculate relative yield between competing reaction based on energy '
                        'difference (DE) only when max(DE) - min(DE) is below given value (kJ/mol). This option make sense only '
                        'when --heatFormation is enable')
    rangeDict['maxAllowDE'] = {'min': 3, 'max': 100, 'default': 20}
    parser.add_argument('--maxAllowDEtr', default=8, type=int, help='')
    rangeDict['maxAllowDEtr'] = {'min': 1, 'max': 50, 'default': 8}
    # 7. correction r_polanyi
    parser.add_argument('--useEvansPolanyi', type=str, default='Y', choices=['Y', 'N'], help='use Evans-Polanyi: Y - use; N - dont use')
    parser.add_argument('--maxAllowMayrNE', default=10, type=int, help='TODO')
    rangeDict['maxAllowMayrNE'] = {'min': 2, 'max': 30, 'default': 20}
    parser.add_argument('--mayrAndDEMode', default='expTexp', help='',
                        choices=['linear', 'exp', 'expabs', 'expsigmoid', 'expKsigmoid', 'expLELU', 'expLOLI', 'expTexp'])
    # 6. correction - r_ring
    parser.add_argument('--boostIntramol', type=float, default=2.0, help='boost reaction with are intramolecular and close '
                        'ring of size 3, 5, 6 or 7 members')
    rangeDict['boostIntramol'] = {'default': 1.5, 'min': 1.0, 'max': 2.8}  # max was 5.0 then 2.8
    parser.add_argument('--mutexRxMode', type=str, default='sidemech', choices=['N', 'VS', 'sidemech'],
                        help='how mutexRx should be treated: N - mutexRx not used; '
                        'sidemech - mutex reaction which are not on main path will be treated as sidemech '
                        'VS - mutex reactions not from main path treated as VS (Very Slow)')

    # ########################################################################
    # fixed ratio fo tautomerisation and acid-base-neutral equlibrium
    parser.add_argument('--fixTautomerRatio', type=str, default='constEnolate', choices=['no', 'constEnolate', 'const', 'DE'],
                        help='Fix ratio of tautomer: no - dont fix, DE - as defined by Boltzman const - as defined in setTautomerRatio, '
                        'constEnolate - as defined in setTautomerRatio and in setTautoEnolateRatio - for enolate tautomerisations')
    parser.add_argument('--tautomerProductIs', type=str, default='lessStable', choices=['lessStable', 'rightHand', 'moreStable'],
                        help='tautomerisation product (as used in --setTautomerRatio) is "lessStable" or "rightHand" side of rx compound')
    parser.add_argument('--setTautomerRatio', default=0.1, type=float,
                        help='when --fixTautomerRatio is set to "const" product/substrate ratio is defined by this option')
    rangeDict['setTautomerRatio'] = {'min': 0.001, 'max': 1.0, 'default': 0.1}
    parser.add_argument('--setTautoEnolateRatio', default=1, type=float, help='enolate/ketone (C anion) ratio')
    parser.add_argument('--fixAcidBaseRatio', type=str, default='condIgnoreWeak', choices=['N', 'Y', 'cond', 'condIgnoreWeak'],
                        help='keep fixed ratio between acid/base and main form')
    parser.add_argument('--setAcidBaseRatio', type=float, default=1.0, help='ratio: acid_or_base / main_form')
    rangeDict['setAcidBaseRatio'] = {'min': 0.01, 'max': 1.5, 'default': 1.0}
    parser.add_argument('--yieldDefMode', type=str, default='maxWithTautoAndDeProt', #AndTautoDeProt',
                        choices=['max', 'maxWithTauto', 'last', 'lastWithTauto', 'lastWithTautoAndDeProt', #'maxWithTautoAndDeProt',
                                 'maxWithTautoAndDeProtAndTautoDeProt'],
                        help='How yield of final product is defined. Allowed values: max - max observed value during simulation; '
                        'maxWithTauto - like max but combine all tautomeric forms of final product; '
                        'last - last observed value (at the end of simmulation) '
                        'lastWithTauto - like "last" but sum up all tautomeric forms of final product; '
                        'lastWithTautoAndDeProt - like "last" but sum up all tautomeric forms and all deprotonated and protonated forms '
                        'maxWithTautoAndDeProt - like "max" but sum up all tautomeric forms and all deprotonated and protonated forms ')
                        # above includes also deprotonated/protonated form of tautomers of target
                        #'maxWithTautoAndDeProtAndTautoDeProt - like "maxWithTautoAndDeProt" but includes also (de)protonated forms of tautomers')

    # ########################################################################
    # finally all done and parse it
    args = parser.parse_args()
    return args, rangeDict


def printWarnings(args):
    if args.propagationMode == 'front':
        print("front propagationMode gives low quality results, use 'all' instead")
    if args.propagationMode == 'copasi':
        print("copasi propagationMode is not fully tested use 'all' mode instead")
    if args.c100 > 100:
        print(f"high value of {args.c100} probably 2c reaction will be underrated")
