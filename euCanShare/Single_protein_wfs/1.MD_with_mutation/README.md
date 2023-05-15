# MD with mutation 

## Quick installation and run

Go to workflow folder and install conda environment:

```bash
export KEY_MODELLER="HERE YOUR MODELLER KEY"
conda env create -f environment.yml
conda activate single_protein_wf1
```

To install it in an HPC environment, do the same after loading the corresponding Conda or Miniconda module.

See options for worklow:

```bash
vi input_HPC.yml
python biobb_md_setup_mutation.py --help
```

See [biobb documentation](https://mmb.irbbarcelona.org/biobb/documentation/source) for additional properties not included in the YAML configuration file.

To run in an HPC environment adapt the run_HPC.sl and input_HPC.yml scripts and send a job to the slurm queue:

```bash
sbatch run_HPC.sl
```

To run locally, modify run_local.sh and input_local.yml if needed:

```bash
./run_local.sh
```

The output will be generated in the "working_dir_path" folder selected in the corresponding YAML input. The global log files will be in "/working_dir_path/log.out" and "/working_dir_path/log.err". Each step will have its own log files and output in a separate folder inside "/working_dir_path".

## Description

This workflow has several steps. The input for the workflow can be (1) a pdb file to be fixed and prepared in steps 1-3 (given through the YAML configuration file or the command line arguments). Or (2) an already prepared gromacs structure file and .zip topology files ready to be minimized (given through command line arguments). The YAML configuration file will contain the default settings and paths of the workflow. The command line arguments can be used to provide some inputs and settings that will be prioritized over those in the YAML configuration file.

- **Step 1**: extraction of structure from PDB. Provide the input pdb file and chain to be extracted through the command line arguments of the workflow or through the paths and properties of step 1 section in the YAML configuration file. The workflow will always prioritize the inputs from command line arguments.

**Steps 2 (A-I)**: steps to fix different possible defects in the input pdb structure. See below.

- **Step 2A**: fix alternative locations. Choose one of the possible alternative positions of residues presenting alternative locations. To avoid this step, comment the 'altlocs' property.

- **Step 2B**: mutations of initial pdb structure. Mutations can be requested through the mutation_list property of the YAML configuration file as a single string of mutations separated by commas (no spaces). Where each mutation is defined by string with the following format: "Chain:Wild_type_residue_name Residue_number Mutated_type_residue_name". The residue name should be a 3 letter code starting with capital letters, e.g. "A:Arg220Ala". Alternatively, they can be given through the mutation_list command line argument. If no mutation is desired leave an empty string ("") or comment the mutation_list property.

- **Step 2C**: download a canonical FASTA file from the Protein Data Bank. This FASTA sequence file is then used to model missing backbone atoms in step 2D. Internet connection and a PDB code are required for this step. This step is executed only if "--fix_backbone" is used in the command line arguments. 

- **Step 2D**: add missing backbone heavy atoms using biobb_structure_checking and Modeller suite. A modeller license key and the previous FASTA file are required for this step. This step is executed only if "--fix_backbone" is used in the command line arguments.  

- **Step 2E**: add missing side chain atoms using biobb_structure_checking (and Modeller suite if a license key is provided).

- **Step 2F**: Add missing disulfide bonds. Use carefully. This step is executed only if "--fix_ss" is used in the command line arguments. 

- **Step 2G**: flip the clashing amide groups to avoid clashes.

- **Step 2H**: fix stereochemical errors in residue side-chains changing their chirality.

- **Step 2I**: renumber atomic indices.

**Steps 3 - 7**: preparation of the system for minimization and MD with GROMACS.

- **Step 3**: uses pdb2gmx to obtain a gromacs structure file (.gro) and topology file from the prepared structure file in step 2. Hydrogen atoms will be added in this step. A force field and water model are chosen with the force_field and water_type properties. 

- **Step 4**: generate a simulation box around the structure.

- **Step 5**: generate a box of solvent around the structure.

- **Step 6**: prepare step 7

- **Step 7**: randomly replace solvent molecules with monoatomic ions. To prepare the system externally, use the "--input_gro" and "--input_top" command line arguments.

**Steps 8 - 10**: energy minimization

**Steps 10 - 13**: NVT equilibration

**Steps 13 - 16**: NPT equilibration

**Steps 17 - 26**: launch several production trajectories from equilibrated structure (see --n_trajs command line argument). Concatenate the imaged trajectories afterwards.

Note that re-launching the workflow will skip the previously successful steps if restart is True and the output folder is the same. 

