import json, pickle
from networkx.classes.coreviews import AtlasView
from networkx.readwrite import json_graph

class SetEncoder(json.JSONEncoder):
    def default(self, obj):
        #print("BOO", obj, type(obj))
        if isinstance(obj, (set, frozenset)):
            retobj = list(obj)
            #print("RET", type(retobj))
            return retobj
        if isinstance(obj, dict):
            hasT = False
            for key in obj:
                if isinstance(key, tuple):
                    print("==", key)
                    hasT = True
            # if hasT:
            return {str(key): str(obj[key]) for key in obj}
        if isinstance(obj, AtlasView):
            obj = {k: obj[k] for k in obj}
            print("BNNN", type(obj), obj)
            return obj
        return json.JSONEncoder.default(self, obj)


def exportRXDBtoJSON(data):
    torm = ('inco', 'incoLiberal', 'rxSmarts', 'reactsSmiles', 'reactive', "allIzolatedFragmentMolInSubstrates", 'conflictSubst', 'smiles',
            'reactingAtomsIdx', 'rx_reactive', 'multistep', 'preFilters', 'aromaticRx', 'absStereo', 'alterByproducts', 'alterCond',
            'tandem', 'fixedSb', 'typicalSubst', 'avgYield', 'yield', 'rxspec', 'multicompUsed', 'solvent', 'aromaticType', 'ringSize', 
            'pmi', 'pmiBest', 'pmiWorst', 'pmiYieldWorst', 'pmiYield', 'pmiYieldBest', 'additionalReagents', 'additionalReagentsForAlter',
            'whenUseFixed', 'mainProductRingAtoms', 'intermol', 'needAlterSolvCond', 'same', 'rxname', 'ref', 'module')
    todelete = []
    for rx in data['rxes']:
        if 'mech' not in data['rxes'][rx]['module']:
            todelete.append(rx)
            continue
        for elem in torm:
            del data['rxes'][rx][elem]
    for rx in todelete:
        del data['rxes'][rx]
    fh = open(f'rxdb,json', 'w')
    json.dump(data, fh, cls=SetEncoder, indent=2)
    fh.close()
    return data


def exportToJSON(data):
    gens = dict()
    for cmd in data:
        #print("GRAPH", type(data[cmd][2]['graph']))
        data[cmd] = formatDict(data[cmd])
        for num in data[cmd]:
            graph = json_graph.node_link_data(data[cmd][num]['graph'])
            print("GR", type(graph))
            data[cmd][num]['graph'] = graph
        graph = json_graph.node_link_data(data[cmd][0]['common']['mainPathInfo']['graph'])
        data[cmd][0]['common']['mainPathInfo']['graph'] = graph
        del data[cmd][0]['data']['competing']['smiMolMap']
        print(cmd, data[cmd].keys())
        fh = open(f'{cmd},json', 'w')
        json.dump(data[cmd], fh, cls=SetEncoder, indent=2)
        fh.close()


def formatDict(dct):
    newDict = dict()
    for k in dct:
        if isinstance(k, (set, tuple, list)):
            newK = str(k)
        else:
            newK = k
        if isinstance(dct[k], dict):
            newV = formatDict(dct[k])
        else:
            newV = dct[k]
        newDict[newK] = newV
    return newDict


def exportToPickle(data, initStoich):
    print("INIT STO", initStoich)
    for cmd in data:
        #print("GRAPH", type(data[cmd][2]['graph']))
        del data[cmd][0]['data']['competing']['smiMolMap']
        pkl = {'data': data[cmd], 'initC': initStoich.get(cmd, None)}
        fh = open(f'{cmd},all.pickle', 'wb')
        pickle.dump(pkl, fh, protocol=0)
        fh.close()


def loadallPickle(fname):
    fh = open(fname, 'rb')
    pkl = pickle.load(fh)
    fh.close()
    # print(fname, "===PKL==", type(pkl), pkl.keys())
    return pkl['data'], pkl['initC']


def loadrxdb(fname):
    fn = open(fname, 'r')
    data = json.load(fn)
    return data
