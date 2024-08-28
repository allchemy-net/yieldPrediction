import sys
from argparse import Namespace


def readOpenboxHistory(files, args):
    for fname in files:
        _parseSingleFile(fname)


def _parseSingleFile(fn):
    fh = open(fn)
    for line in fh:
        if line.startswith(' {'):
            results, space = line.split(' Namespace')
            results = eval(results)
            space = eval(f'Namespace{space}')
            print(results.keys(), type(space))
    fh.close()


if __name__ == "__main__":
    readOpenboxHistory(sys.argv[1:], None)
