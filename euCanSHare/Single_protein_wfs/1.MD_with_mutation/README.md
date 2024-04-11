# MD with mutation 

![alt text](../../img/MD_setup.png?raw=true)

## Quick installation and run

Go to workflow folder and install the conda environment (running in Nostrum's cluster use the already installed environments located in */shared/work/BiobbWorkflows/envs*)

```bash
export KEY_MODELLER="HERE YOUR MODELLER KEY"
conda env create -f environment.yml
conda activate biobb_sp_md
```

See [biobb documentation](https://mmb.irbbarcelona.org/biobb/documentation/source) for additional properties not included in the YAML configuration file.

To run in an HPC environment adapt the run_HPC.sl and input_HPC.yml scripts and send a job to the slurm queue:

```bash
sbatch run_HPC.sl
```

By default, the output will be generated in the "working_dir_path" folder selected in the YAML configuration file. However, the "--output" command line option will overwrite "working_dir_path". The global log files will be in "output/log.out" and "output/log.err". Each step will have its own log files and output in a separate folder inside the output folder.

## Inputs

### Configuration file

Take a look at the YAML configuration file to see the different properties that can be set.

```bash
vi input_HPC.yml
```

Specially important are: the binary path of GROMACS, the input files and the MODELLER key. Make sure the binary path specified and the module loaded in the run file (HPC only) agree between them.

### Command line arguments

The command line arguments can be used to provide some inputs and settings that will be prioritized over those in the YAML configuration file.

```bash
python biobb_md_setup_mutation.py --help
```

Specially important are: the configuration file path, the input pdb file or the input topology (.zip) and coordinates file (.gro), the number of trajectories to launch and the output path.

## Description

This workflow has several steps. The input for the workflow can be (1) a pdb file to be fixed and prepared. Or (2) an already prepared gromacs structure file and .zip topology files ready to be minimized.

- **Step 1**: extraction of structure from PDB. Provide the input pdb file and chain to be extracted through the command line arguments of the workflow or through the paths and properties of step 1 section in the YAML configuration file. The workflow will always prioritize the inputs from command line arguments.

**Steps 2 (A-I)**: steps to fix different possible defects in the input pdb structure. See below.

- **Step 2A**: fix alternative locations. Provide a list with the choices of alternative locations to keep in the final structure.

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

To make sure the system has been correctly prepared before minimizing or running MD, launch the workflow adding the '--setup_only' command line option. This will stop the workflow before the energy minimization. 

**Steps 8 - 10**: energy minimization (including position restraints on the proteins heavy atoms)

**Steps 10 - 13**: NVT equilibration (including position restraints on the proteins heavy atoms)

**Steps 14 - 16**: NPT equilibration (including position restraints on the proteins heavy atoms)

**Steps 17 - 26**: launch several production trajectories from equilibrated structure (see --n_trajs command line argument). Concatenate the imaged trajectories afterwards. Computation of RMSD (with fitting) with respect to experimental structure and with respect to equilibrated structure (protein backbone atoms). Computation of Radius of gyration (protein backbone atoms) and RMSF (protein heavy atoms).

Note that re-launching the workflow will skip the previously successful steps if restart is True and the output folder is the same. 

