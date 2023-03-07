#!/usr/bin/env python

from rdkit import Chem
from espsim import EmbedAlignScore
import argparse,sys

def parse_command_line(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--target', help='target molecule', required=True)
    parser.add_argument('--test', help='query molecule', required=True)
    return parser.parse_args()

args = parse_command_line(sys.argv)

target_mol = args.target
test_mol = args.test

data = [m for m in Chem.SDMolSupplier(test_mol, removeHs=False)]
target = [n for n in Chem.SDMolSupplier(target_mol, removeHs=False)]

shapesim, espsim = EmbedAlignScore(target[0], data, renormalize=True,
                           metric='tanimoto')
out=(shapesim, espsim)

with open('output.txt','w') as f:
	f.write(str(out))
	f.close()
