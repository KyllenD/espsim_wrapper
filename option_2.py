#!/usr/bin/env python

from rdkit import Chem
from rdkit.Chem import rdMolAlign
from rdkit.Chem import rdMolDescriptors
from espsim import GetEspSim, GetShapeSim
import argparse,sys

def parse_command_line(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--target', help='target molecule', required=True)
    parser.add_argument('--test', help='query molecule', required=True)
    return parser.parse_args()

def align(prbMol, refMol, prbCrippen=None, refCrippen=None, i=-1, j=-1):
    if prbCrippen is None:
        prbCrippen = rdMolDescriptors._CalcCrippenContribs(prbMol)
    if refCrippen is None:
        refCrippen = rdMolDescriptors._CalcCrippenContribs(refMol)
    alignment = rdMolAlign.GetCrippenO3A(prbMol, refMol, prbCrippen, refCrippen, i, j)
    alignment.Align()
    
args = parse_command_line(sys.argv)

target_mol = args.target
test_mol = args.test

data = [m for m in Chem.SDMolSupplier(test_mol, removeHs=False)]
target = [n for n in Chem.SDMolSupplier(target_mol, removeHs=False)]

with open('output.txt','w') as f:
  for mol in data:
    align(mol, target[0])
    espsim = GetEspSim(mol, target[0], renormalize=True, metric='tanimoto')
    shapesim = GetShapeSim(mol, target[0])
    out=(shapesim, espsim)
    f.write(str(out))
  f.close()
