The files aims to transform the files containing the SMILES (Simplified molecular-input line-entry system) of molecules into the fragments counts for isodesmic 
reaction establization and energy prediction by pretrained group additivity model.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Directly run 'file_transformation.py' to transform the SMILE file.

In the python file, the line 'file_transform('test_smiles.txt',atoms, mode = 'regular')' allows for different configuration.

1) change input files by modifying 'test_smiles.txt'.

2) change fragments by altering 'mode = 'regular'' to 'hybridization' or 'pathway fps'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prerequisites:

rdkit module is required for the fragments counting. 
https://www.rdkit.org/docs/Install.html

To install rdkit
type 'conda install -c https://conda.anaconda.org/rdkit rdkit' in the command prompt.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Two database are included in the 'database' folder for benchmarking ('NIST_ene.txt') or isodesmic reaction extension ('QM7.txt').

1) 'NIST_ene.txt' file contains different molecules and their experimental heat of formation energies.

2) 'QM7.txt' file contains molecules consisted of C, O, H atoms from GDB 7 database.

Training model is built with group additiviy and LASSO based on molecules from QM9 database. Especially, molecules with multiple rings or with heavy atoms more than 9 are excluded, which leads to 22k molecules as training set.