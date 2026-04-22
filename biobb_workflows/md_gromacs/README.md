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
vi input.yml
```

Specially important are: the binary path of GROMACS and the MODELLER key. Make sure the binary path specified and the module loaded in the run file (HPC only) agree between them.

### Command line arguments

The command line arguments can be used to provide some inputs and settings that will be prioritized over those in the YAML configuration file.

```bash
python workflow.py --help
```

```

```

## Description

This workflow has several steps. The input for the workflow can be (1) a prepared PDB file or (2) a gromacs structure file and .zip topology files ready to be minimized.

3. **Preparation of topology and coordinates for MD**

    A. **pdb2gmx** 
    Uses pdb2gmx to obtain a gromacs structure file (.gro) and topology file from the fixed PDB. Hydrogen atoms will be added in this step, one can choose to ignore the hydrogens in the original structure or not (```ignh``` property). The protonation state of histidines can be provided (```his``` property) in the form of a list of numbers see below. A force field and water model are chosen here.
    
    For the ```his``` property include a string with the protonation states '0 0 1 1 0 0 0', where:

        - 0 : H on ND1 only (HID)
        - 1 : H on NE2 only (HIE)
        - 2 : H on ND1 and NE2 (HIP)
        - 3 : Coupled to Heme (HIS1)
    
    NOTE: default behavior is to add charged termini - if one wants ACE and NME it should be provided already with the correct atom names - look at the force field being used: /eb/x86_64/software/GROMACS/2023.3-foss-2022a-CUDA-11.7.0-PLUMED-2.9.0/share/gromacs/top/amber99sb-ildn.ff/aminoacids.rtp
    
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

