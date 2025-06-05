# espsim_wrapper
Scipts required to run esp-shape screening

Requires a conda env with:
- espsim (https://github.com/hesther/espsim)
- openbabel

Copy the files in this repo and run 

`bash main_run.sh target.sdf database.sdf flexible/rigid/rdkit metric(tanimoto/carbo) chargemethod(gasteiger/mmff/espaloma) sortmethod(ESP*Shape/ESP/Shape) extract(yes/no) numbertoextract`

Mol2 files can also be used as inputs for target.sdf and database.sdf
Flexible option performs flexible alignemnt with LSalign, rigid performs rigid alignment with LSalign and rdkit performs rigid alignment using crippen alignemnt in rdkit
Currently only gasteiger charges are supported
Extarction options take the top number of molecules based on the sortmethod selected
