# MD with mutation 

![alt text](../../img/MD_setup.png?raw=true)

## Quick installation and run

Go to workflow folder and install the conda environment (running in Nostrum's cluster use the already installed environments located in */shared/work/BiobbWorkflows/envs*)

```bash
export KEY_MODELLER="HERE YOUR MODELLER KEY"
conda env create -f environment.yml
conda activate biobb_md
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

1. **Extraction of structure from PDB**
    Provide the input pdb file and chain to be extracted through the command line arguments of the workflow or through the paths and properties of step 1 section in the YAML configuration file. The workflow will always prioritize the inputs from command line arguments.

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

**Steps 3 - 7**: preparation of the system for minimization and MD with GROMACS.

3. **pdb2gmx** 
Uses pdb2gmx to obtain a gromacs structure file (.gro) and topology file from the fixed PDB. Hydrogen atoms will be added in this step, one can choose to ignore the hydrogens in the original structure or not (```ignh``` property). The protonation state of histidines can be provided (```his``` property) in the form of a list of numbers see below. A force field and water model are chosen here.
    
    For the ```his``` property include a string with the protonation states '0 0 1 1 0 0 0', where:

        - 0 : H on ND1 only (HID)
        - 1 : H on NE2 only (HIE)
        - 2 : H on ND1 and NE2 (HIP)
        - 3 : Coupled to Heme (HIS1)
    
    NOTE: this will be automatized within the workflow in the future, right now pdb4amber can be used externally to obtain a guess

4. **Generate simulation box** 
Generate a simulation box around the structure. Some box types are more efficient than others (octahedron for globular proteins)

5. **Solvation** 
Generate a box of solvent around the structure.

6. **Prepare ion addition**
Prepares the next step

7. **Add ions** 
Randomly replace solvent molecules with monoatomic ions. 

To prepare the system externally, use the ```--input_gro``` and ```--input_top``` command line arguments.

To make sure the system has been correctly prepared before minimizing or running MD, launch the workflow adding the ```--setup_only``` command line option. This will stop the workflow before the energy minimization. 

**Steps 8 - 10**: energy minimization (including position restraints on the proteins heavy atoms)

**Steps 10 - 13**: NVT equilibration (including position restraints on the proteins heavy atoms)

**Steps 14 - 16**: NPT equilibration (including position restraints on the proteins heavy atoms)

**Steps 17 - 26**: launch several production trajectories from equilibrated structure (see --n_trajs command line argument). Concatenate the imaged trajectories afterwards. Computation of RMSD (with fitting) with respect to experimental structure and with respect to equilibrated structure (protein backbone atoms). Computation of Radius of gyration (protein backbone atoms) and RMSF (protein heavy atoms).

Note that re-launching the workflow will skip the previously successful steps if restart is True and the output folder is the same. 

