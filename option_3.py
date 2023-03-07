#!/usr/bin/env python

import sys
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import espsim
from espsim import EmbedAlignConstrainedScore, readSdfFile
import pandas as pd
import argparse,sys

def parse_command_line(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--target', help='target molecule', required=True)
    parser.add_argument('--test', help='query molecule', required=True)
    return parser.parse_args()

def common_core(mols):
    #Here you can customize the rdFMCS parameters (such as setting bond order comparisons to exact), feel free to change this:
    atomComp = rdFMCS.AtomCompare.CompareElements

    atomParam = rdFMCS.MCSAtomCompareParameters()
    atomParam.MatchChiralTag = False
    atomParam.MatchFormalCharge = True
    atomParam.MatchValences = False
    atomParam.RingMatchesRingOnly = False

    bondComp = rdFMCS.BondCompare.CompareOrderExact

    bondParam = rdFMCS.MCSBondCompareParameters()
    bondParam.CompleteRingsOnly = True
    bondParam.RingMatchesRingOnly = False

    opt = rdFMCS.MCSParameters()
    opt.MaximizeBonds = True
    opt.Timeout = 3600
    opt.AtomCompareParameters = atomParam
    opt.BondCompareParameters = bondParam
    opt.SetAtomTyper(atomComp)
    opt.SetBondTyper(bondComp)
    opt.Threshold = 1.0
    opt.Verbose = False
    opt.seedSmarts = ""

    mcs = rdFMCS.FindMCS(mols, opt)
    return mcs

args = parse_command_line(sys.argv)

target_mol = args.target
test_mol = args.test

data = [m for m in Chem.SDMolSupplier(test_mol, removeHs=False)]
target = [n for n in Chem.SDMolSupplier(target_mol, removeHs=False)]

#Only embed target once:
AllChem.EmbedMolecule(target[0], AllChem.ETKDGv2()) #Acutally, there are coordinates in the SDF, if you want to use those, comment out this line and the next!
AllChem.UFFOptimizeMolecule(target[0]) 

with open('output.txt','w') as f:
  for mol in data:
    mols = target + [mol]
    mcs = common_core(mols)
    mcsmol = Chem.MolFromSmarts(mcs.smartsString)
    core = AllChem.DeleteSubstructs(AllChem.ReplaceSidechains(mols[0],mcsmol),Chem.MolFromSmiles('*'))
    shapesim, espsim = EmbedAlignConstrainedScore(mols[0], [mols[1]], core, renormalize=True,metric='tanimoto') # I again omitted setting NumConfs to 1, because I think you don't want that
    out=(shapesim, espsim)
    f.write(str(out))
  f.close()

