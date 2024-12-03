# MD with GROMACS

This workflow uses BioBBs to fix PDB defects, prepare the MD simulations, equilibrate and execute production runs, do basic analysis of the trajectories and post-process the trajectories.

![alt text](../../img/MD_setup.png?raw=true)

## Quick installation and run

---

Go to the workflow folder and install the conda environment (running in Nostrum's cluster use the already installed environments located in */shared/work/BiobbWorkflows/envs*)

```bash
export KEY_MODELLER="HERE YOUR MODELLER KEY"
conda env create -f environment.yml
conda activate biobb_md
```

See [biobb documentation](https://mmb.irbbarcelona.org/biobb/documentation/source) for additional properties not included in the YAML configuration file.

To run a single call to the workflow in an HPC environment use:

```bash
sbatch run_HPC.sl
```

To run a long MD simulation while respecting the time limit of the HPC jobs use:

```bash
./launch_long_MD.sh
```

## Inputs

---

### Configuration file

Take a look at the YAML configuration file to see the different properties that can be set.

```bash
vi input_HPC.yml
```

Specially important are: the binary path of GROMACS and the MODELLER key. Make sure the binary path specified and the module loaded in the run file (HPC only) agree between them.

### Command line arguments

The command line arguments can be used to provide some inputs and settings that will be prioritized over those in the YAML configuration file.

```bash
python biobb_md_setup_mutation.py --help
```

```
usage: MD Simulation with GROMACS [-h] --config CONFIG_PATH [--setup_only] [--num_parts NUM_PARTS] [--num_replicas NUM_REPLICAS] [--output OUTPUT_PATH] [--input_pdb INPUT_PDB_PATH] [--his HIS]
                                  [--pdb_chains PDB_CHAINS [PDB_CHAINS ...]] [--mutation_list MUTATION_LIST [MUTATION_LIST ...]] [--input_gro INPUT_GRO_PATH] [--input_top INPUT_TOP_PATH]
                                  [--skip_fix_backbone] [--fix_ss] [--fix_amides] [--nsteps NSTEPS] [--final_analysis] [--ligand_parameters LIGAND_PARAMETERS]

options:
  -h, --help            show this help message and exit
  --config CONFIG_PATH  Configuration file (YAML)
  --setup_only          Only setup the system (default: False)
  --num_parts NUM_PARTS
                        Number of parts of the trajectorie (default: 1)
  --num_replicas NUM_REPLICAS
                        Number of replicas (default: 1)
  --output OUTPUT_PATH  Output path (default: working_dir_path in YAML config file)
  --input_pdb INPUT_PDB_PATH
                        Input PDB file (default: input_structure_path in step 1 of configuration file)
  --his HIS             Histidine protonation states list.
  --pdb_chains PDB_CHAINS [PDB_CHAINS ...]
                        PDB chains to be extracted from PDB file (default: chains in properties of step 1)
  --mutation_list MUTATION_LIST [MUTATION_LIST ...]
                        List of mutations to be introduced in the protein (default: None, ex: A:Arg220Ala)
  --input_gro INPUT_GRO_PATH
                        Input structure file ready to minimize (.gro). To provide an externally prepared system, use together with --input_top (default: None)
  --input_top INPUT_TOP_PATH
                        Input compressed topology file ready to minimize (.zip). To provide an externally prepared system, use together with --input_gro (default: None)
  --skip_fix_backbone   Skip the backbone modeling of missing atoms (default: False)
  --fix_ss              Add disulfide bonds to the protein. Use carefully! (default: False)
  --fix_amides          Flip clashing amides to relieve the clashes (default: False)
  --nsteps NSTEPS       Number of steps of the simulation
  --final_analysis      Run the final analysis of the trajectory/ies. Concatenation of the analysis and trajectory, trajectory drying, imaging and fitting (default: False)
  --ligand_parameters LIGAND_PARAMETERS
                        Folder with .itp files for the ligand topology and heavy-atom constraints (default: None). If there is a ligand in the selected chain with a topology in the ligand_parameters
                        folder, it will be included in the simulation.
```

## Description

This workflow has several steps. The input for the workflow can be (1) a pdb file to be fixed and prepared. Or (2) an already prepared gromacs structure file and .zip topology files ready to be minimized.

1. **Extraction of structure from PDB**
    Provide the input pdb file and chain to be extracted through the command line arguments or through the paths and properties of section 1 in the YAML configuration file. The workflow will always prioritize the inputs from command line arguments. Any cofactor present in the selected chain will also be extracted.

2. **Fix PDB defects (A-I)**
    Steps to fix different possible defects in the input pdb structure. See below.

    A. **Fix alternative locations** 
    Provide a list with the choices of alternative locations to keep in the final structure. If no list is given (_null_ value) it will select the alternative location with the highest occupancy (the workflow will use Biopython to do so). 

    B. **Mutate initial pdb structure** 
    Mutations can be requested through the mutation_list property of the YAML configuration file as a single string of mutations separated by commas (no spaces). Where each mutation is defined by string with the following format: "Chain:Wild_type_residue_name Residue_number Mutated_type_residue_name". The residue name should be a 3 letter code starting with capital letters, e.g. "A:Arg220Ala". Alternatively, they can be given through the mutation_list command line argument. If no mutation is desired leave an empty string ("") or comment the mutation_list property.

    C. **Obtain the Sequence in FASTA format** 
    The sequence is then used to model missing backbone atoms in the next step. The workflow first tries to download the canonical FASTA (including all residues for that protein) from the Protein Data Bank. If there is no internet connection, it will try to obtain the sequence from the _SEQRES_ records in the PDB. If there are no _SEQRES_, then only the residues that contain at least one atom in the structure will be included. This step can be skipped including the ```--skip_fix_backbone``` option.  

    D. **Model missing backbone atoms**
    Add missing backbone heavy atoms using _biobb_structure_checking_ and Modeller suite. A modeller license key and the previous FASTA file are required for this step. This step can be skipped including the ```--skip_fix_backbone``` option.  

    E. **Model missing side chain atoms**
    Add missing side chain atoms using _biobb_structure_checking_ (and Modeller suite if a license key is provided).

    F. **Add missing disulfide bonds**
    It changes CYS for CYX to mark cysteines residues pertaining to a [di-sulfide bond](https://en.wikipedia.org/wiki/Disulfide). It uses a distance criteria to determine if nearby cysteines are part of a di-sulfide bridge (_check_structure getss_). Use carefully. This step is executed only if the "--fix_ss" option is used. 

    G. **Relieve clashes flipping amide groups**
    It flips the clashing amide groups to relieve clashes.

    H. **Fix chirality of residues**
    Creates a new PDB file fixing stereochemical errors in residue side-chains changing it's chirality when needed.

    I. **Renumber atomic indices**
    So they start at 1.

3. **Preparation of topology and coordinates for MD**

    A. **pdb2gmx** 
    Uses pdb2gmx to obtain a gromacs structure file (.gro) and topology file from the fixed PDB. Hydrogen atoms will be added in this step, one can choose to ignore the hydrogens in the original structure or not (```ignh``` property). The protonation state of histidines can be provided (```his``` property) in the form of a list of numbers see below. A force field and water model are chosen here.
    
    For the ```his``` property include a string with the protonation states '0 0 1 1 0 0 0', where:

        - 0 : H on ND1 only (HID)
        - 1 : H on NE2 only (HIE)
        - 2 : H on ND1 and NE2 (HIP)
        - 3 : Coupled to Heme (HIS1)
    
    NOTE: this will be automatized within the workflow in the future, right now pdb4amber can be used externally to obtain a guess

    B. **Generate cofactors topology**: if there are any cofactors with parameters in the `--ligand_parameters` folder, use _tleap_ to build the corresponding AMBER topology and coordinates file.

    C. **Convert AMBER topology and coordinates to GMX**: convert AMBER topology and coordinate files to GROMACS format.

    D-F: **Merge cofactors and structure**: if any parameterized cofactors are present, convert the structures to PDB format and concatenate them.

    G-H. **Generate restraints for cofactor heavy atoms**.

    I. **Append cofactor topology to structure topology**.

    J. **Generate simulation box** 
    Generate a simulation box around the structure. Some box types are more efficient than others (octahedron for globular proteins)

    K. **Solvation** 
    Generate a box of solvent around the structure.

    L-M. **Add ions** 
    Randomly replace solvent molecules with monoatomic ions. 

    N. **Convert final topology to PDB**
    This will be used in subsequent post-processing steps.

To prepare the system externally, use the ```--input_gro``` and ```--input_top``` command line arguments.

To make sure the system has been correctly prepared before minimizing or running MD, launch the workflow adding the ```--setup_only``` command line option. This will stop the workflow before the energy minimization. 

4. **Minimize and equilibrate the initial configuration**

    A-D. **Energy minimization** (including position restraints on the protein heavy atoms)

    E-G. **NVT equilibration** (including position restraints on the protein heavy atoms)

    H-J. **NPT equilibration** (including position restraints on the protein heavy atoms)

5. **Production run**

    Launch several production trajectories from equilibrated structure (see `--num_parts` or `--num_replica` command line argument).

6. **Basic analysis**

    Computation of RMSD (with fitting) with respect to experimental structure and with respect to equilibrated structure (protein backbone atoms). Computation of Radius of gyration (protein backbone atoms) and RMSF (protein heavy atoms).

7. **Post-processing**

    Concatenate if the trajectories are parts. Image, dry and fitting.

Note that re-launching the workflow will skip the previously successful steps if restart is True and the output folder is the same. 

