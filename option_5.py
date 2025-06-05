#/usr/bin/env python

from rdkit import Chem
from rdkit.Chem import rdMolAlign
from rdkit.Chem import rdMolDescriptors,rdmolops
from rdkit.Chem import AllChem
from espsim import GetEspSim, GetShapeSim, EmbedAlignScore
from espsim.electrostatics import GetIntegralsViaGaussians,GetIntegralsViaMC,Renormalize
from espsim.electrostatics import GetMolProps
import mmap
from openbabel import pybel, openbabel as ob 
import argparse,sys
import pandas as pd 
import os
import numpy as np
import psi4

def mmap_read_sdf(filename):
    with open(filename, 'r+b') as f:
        mm = mmap.mmap(f.fileno(), 0)
        supplier = Chem.ForwardSDMolSupplier(mm,sanitize=False,removeHs=False)
        return [mol for mol in supplier if mol]

def parse_command_line(argv):
	parser = argparse.ArgumentParser()
	parser.add_argument('--target', help='target molecule', required=True)
	parser.add_argument('--metric', help='similarity metric (tanimoto/carbo)', required=True)
	parser.add_argument('--chargemethod', help='Partial charge distribution to use: Gasteiger (default), mmff, ml, RESP', required=False)
	parser.add_argument('--sort', help='sort output by (ESP/Shape/ESP*Shape similarity)', required=True)
	parser.add_argument('--sanitize', help='sanitize input or not', required=False)
	parser.add_argument('--outfile', help='galaxy_ouput', required=False)
	parser.add_argument('--method', help='QM method', required=False)
	parser.add_argument('--basis', help='QM basis set', required=False)
	return parser.parse_args()

def crip_align(prbMol, refMol, prbCrippen=None, refCrippen=None, i=-1, j=-1):
	if prbCrippen is None:
		prbCrippen = rdMolDescriptors._CalcCrippenContribs(prbMol)
	if refCrippen is None:
		refCrippen = rdMolDescriptors._CalcCrippenContribs(refMol)
	alignment = rdMolAlign.GetCrippenO3A(prbMol, refMol, prbCrippen, refCrippen, i, j)
	alignment.Align()
    
args = parse_command_line(sys.argv)

target_mol = args.target
test_mol=args.target
met= args.metric
charge=args.chargemethod
sortby=str(args.sort)
san_opt=bool(args.sanitize)
qm_meth=args.method
basis=args.basis
outfile=args.outfile


if sortby=="ESP":
	sortk=2
elif sortby=="Shape":
	sortk=3
else:
	sortk=4

data = mmap_read_sdf(target_mol)

if charge=="psi4":
    psi4.set_memory('240 GB')
    psi4.set_num_threads(16)
    psi4.set_options({'geom_maxiter':5000,'maxiter':5000})
    refCoor,refCharge,refmol2=GetMolProps(data[0],0,[],'gasteiger')
    mol=data[0]
    mol_block=Chem.MolToMolBlock(mol)
    pybel_mol = pybel.readstring("mol", mol_block)
    #ref_netcharge=round(sum(pybel_mol.calccharges()))
    xyz=Chem.rdmolfiles.MolToXYZBlock(mol,confId=0)
    refCharge=psi4Charges(xyz,pybel_mol.charge,pybel_mol.spin,basis,qm_meth,6)
    refVdw = np.array([Chem.GetPeriodicTable().GetRvdw(a.GetAtomicNum()) for a in renum_ref.GetAtoms()]).reshape(-1,1)

if charge=="gasteiger":
    refCoor,refCharge=GetMolProps(data[0],0,[],'gasteiger')
    molblock = Chem.MolToMolBlock(data[0])
    obmol = pybel.ob.OBMol()
    obconv = pybel.ob.OBConversion()
    obconv.SetInFormat("mol")
    obconv.ReadString(obmol, molblock)
    pybel_mol = pybel.Molecule(obmol)
    charge_model = ob.OBChargeModel.FindType("Gasteiger")
    charge_model.ComputeCharges(pybel_mol.OBMol)
    refCharge=np.array([atom.partialcharge for atom in pybel_mol])

for mol in data:
    crip_align(mol,data[0])
    name=mol.GetProp("_Name")
    shapesim = GetShapeSim(mol, data[0])
    if charge=="psi4":
        prbCoor,prbCharge=GetMolProps(mol,0,[],'gasteiger')
        mol=renum_mol
        mol_block=Chem.MolToMolBlock(mol)
        pybel_mol = pybel.readstring("mol", mol_block)
        xyz=Chem.rdmolfiles.MolToXYZBlock(mol,confId=0)
        #prb_netcharge=round(sum(pybel_mol.calccharges()))
        prbCharge=psi4Charges(xyz,pybel_mol.charge,pybel_mol.spin,basis,qm_meth,6)
        for i, atom in enumerate(pybel_mol.atoms):
            atom.OBAtom.SetPartialCharge(prbCharge[i])
        pybel_mol.title=name
        mol2_content = pybel_mol.write(format="mol2")
        with open(name+'_resp_charged.mol2', "w") as f:
            f.write(mol2_content)
        prbVdw = np.array([Chem.GetPeriodicTable().GetRvdw(a.GetAtomicNum()) for a in mol.GetAtoms()]).reshape(-1,1)
        #similarity=GetIntegralsViaMC(prbCoor,refCoor,prbCharge,refCharge,prbVdw,refVdw,met,10,1)
        similarity=GetIntegralsViaGaussians(prbCoor,refCoor,prbCharge,refCharge,met)
        espsim=Renormalize(similarity,met,None)
    else:
        prbCoor,prbCharge=GetMolProps(mol,0,[],'gasteiger')
        molblock = Chem.MolToMolBlock(mol)
        obmol = pybel.ob.OBMol()
        obconv = pybel.ob.OBConversion()
        obconv.SetInFormat("mol")
        obconv.ReadString(obmol, molblock)
        pybel_mol = pybel.Molecule(obmol)
        charge_model = ob.OBChargeModel.FindType("Gasteiger")
        charge_model.ComputeCharges(pybel_mol.OBMol)
        prbCharge=np.array([atom.partialcharge for atom in pybel_mol])
        similarity=GetIntegralsViaGaussians(prbCoor,refCoor,prbCharge,refCharge,met)
        espsim=Renormalize(similarity,met,None)
		#espsim = GetEspSim(mol, data[0], renormalize=True, metric=met, partialCharges=charge)
    product = espsim * shapesim
    os.system("echo {0},{1},{2},{3} >> {4}".format(name,espsim,shapesim,product,outfile))

os.system("cat {0} | sed '/nan/d' | sort -g -k {1} -t , -r > sorted.csv".format(outfile,sortk))
os.system("sed -i '1s/^/Name,ESP Similarity,Shape Similarity,ESP_Shape Similarity\\n/' sorted.csv")

