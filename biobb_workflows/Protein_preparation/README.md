# Protein preparation

This workflow can be used to fix PDB defects, choose protonation states for tritatable residues and prepare the system for simulation. 

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
sbatch run_wf.sl
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
usage: MD Simulation with GROMACS [-h] --config CONFIG_PATH [--input_pdb INPUT_PDB_PATH] [--pdb_code PDB_CODE] [--pdb_chains PDB_CHAINS [PDB_CHAINS ...]]
                                  [--mutation_list MUTATION_LIST [MUTATION_LIST ...]] [--skip_bc_fix] [--modeller_key MODELLER_KEY] [--cap_ter]
                                  [--skip_sc_fix] [--skip_ss_bonds] [--skip_amides_flip] [--ph PH] [--his HIS [HIS ...]] [--keep_hs]
                                  [--pdb_format PDB_FORMAT] [--output OUTPUT_PATH]

options:
  -h, --help            show this help message and exit
  --config CONFIG_PATH  Configuration file (YAML)
  --input_pdb INPUT_PDB_PATH
                        Input PDB file. If not given the workflow will look for it in the YAML config file. Default: None
  --pdb_code PDB_CODE   PDB code to get the canonical FASTA sequence of the input PDB file. If not given the workflow will look for it in the HEADER of the
                        PDB. Default: None
  --pdb_chains PDB_CHAINS [PDB_CHAINS ...]
                        Protein PDB chains to be extracted from PDB file and fixed. Default: A. Example: A B C
  --mutation_list MUTATION_LIST [MUTATION_LIST ...]
                        List of mutations to be introduced in the protein. Default: None. Example: A:Arg220Ala B:Ser221Gly
  --skip_bc_fix         Skip the backbone modeling of missing atoms. Otherwise the missing atoms in the backbone of the PDB structure will be modeled using
                        'biobb_structure_checking' and the 'Modeller suite' (if the Modeller key is given). Note that modeling of missing loops is only
                        possible if the Modeller key is provided. To obtain one register at: https://salilab.org/modeller/registration.html. Default: False
  --modeller_key MODELLER_KEY
                        Modeller key to be used for the backbone modeling of missing atoms. Note that modeling of missing loops is only possible if the
                        Modeller key is provided. To obtain one register at: https://salilab.org/modeller/registration.html
  --cap_ter             Add terminal residues ACE and NME as necessary, preserving existing atoms. Default: False
  --skip_sc_fix         Skip the side chain modeling of missing atoms. Otherwise the missing atoms in the side chains of the PDB structure will be modeled
                        using 'biobb_structure_checking' and the 'Modeller suite' (if the Modeller key is given). Default: False
  --skip_ss_bonds       Skip the addition of disulfide bonds to the protein according to a distance criteria. Otherwise the missing atoms in the side chains
                        of the PDB structure will be modeled using 'biobb_structure_checking' and the 'Modeller suite' (if the Modeller key is given).
                        Default: False
  --skip_amides_flip    Skip the fliping of clashing amide groups in ASP or GLU residues. Otherwise the amide orientations will be changed if needed to
                        relieve clashes using 'biobb_structure_checking'. Note that amide group orientations coming from PDB structures is not reliable in
                        general due to symmetries in the electron density. Default: False
  --ph PH               pH of the system. Used together with a pKa estimation (with propka) to determine the protonation state of titratable residues.
                        Default: 7.0
  --his HIS [HIS ...]   Manual selection of histidine protonation states (delta: 0, epsilon: 1, fully protonated: 2, bound to heme: 3). If given, the pKa
                        estimation and the pH won't be used to protonate histidine residues. Default: None. Example: 0 1 1
  --keep_hs             Keep hydrogen atoms in the input PDB file. Otherwise they will be discarded. Default: False
  --pdb_format PDB_FORMAT
                        PDB format to be used. Options: amber, gromacs. Default: 'amber'
  --output OUTPUT_PATH  Output path. Default: 'output' in the current working directory
```

## Description

1. **Extraction of structure from PDB**

2. **Fix PDB defects (A-I)**
    Steps to fix different possible defects in the input pdb structure. See below.

    1. **Fix alternative locations** 
    Provide a list with the choices of alternative locations to keep in the final structure. If no list is given (_null_ value) it will select the alternative location with the highest occupancy (the workflow will use Biopython to do so). 

    2. **Mutate initial pdb structure** 
    Mutations can be requested through the mutation_list command line argument. Where each mutation is defined by string with the following format: "Chain:Wild_type_residue_name Residue_number Mutated_type_residue_name". The residue name should be a 3 letter code starting with capital letters, e.g. "A:Arg220Ala".

    3. **Obtain the Sequence in FASTA format** 
    The sequence is then used to model missing backbone atoms in the next step. The workflow first tries to download the canonical FASTA (including all residues for that protein) from the Protein Data Bank. If there is no internet connection, it will try to obtain the sequence from the _SEQRES_ records in the PDB. If there are no _SEQRES_, then only the residues that contain at least one atom in the structure will be included. This step can be skipped including the ```--skip_bc_fix``` option.  

    4. **Model missing backbone atoms**
    Add missing backbone heavy atoms using _biobb_structure_checking_ and Modeller suite. A modeller license key and the previous FASTA file are required for this step. This step can be skipped including the ```--skip_bc_fix``` option.  

    5. **Model missing side chain atoms**
    Add missing side chain atoms using _biobb_structure_checking_ (and Modeller suite if a license key is provided).

    6. **Renumber atomic indices**
    So they start at 1.

    7. **Relieve clashes flipping amide groups**
    It flips the clashing amide groups to relieve clashes.

    8. **Fix chirality of residues**
    Creates a new PDB file fixing stereochemical errors in residue side-chains changing it's chirality when needed.

    9. **Add missing disulfide bonds**
    It changes CYS for CYX to mark cysteines residues pertaining to a [di-sulfide bond](https://en.wikipedia.org/wiki/Disulfide). It uses a distance criteria to determine if nearby cysteines are part of a di-sulfide bridge (_check_structure getss_). Use carefully, this step can be skipped using ```--skip_ss_bonds```

    10. **Remove all hydrogens**

    11. **Estimate the pKa of titratable residues with propka**
    Use an empirical method to estimate the pKa of residues considering the local environment. 

    12. **Optimize the hydrogen placement in Histidines with reduce**
    Use reduce from AmberTools to optimize the hydrogen bonds of histidines.

    13. **Select protonation state for titratable residues** 
    Taking into account previous two steps and pH

    14. **Add hydrogens back**
    If output format is amber.