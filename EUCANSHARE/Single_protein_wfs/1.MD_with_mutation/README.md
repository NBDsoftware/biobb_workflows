# MD with mutation 

## Quick installation and run

Go to workflow folder and install conda environment:

```bash
export KEY_MODELLER="HERE YOUR MODELLER KEY"
conda env create -f environment.yml
conda activate single_protein_wf1
```

To install it in an HPC environment, do the same after loading the corresponding module with Conda or Miniconda installed.

See options for worklow:

```bash
vi input.yml
python biobb_md_setup_mutation.py --help
```

See [biobb documentation](https://mmb.irbbarcelona.org/biobb/documentation/source) for additional properties not included in input.yml.

To run in an HPC environment adapt the run_HPC.sl script and send job to slurm queue:

```bash
sbatch run_HPC.sl
```

To run locally, modify run_local.sh if needed and launch:

```bash
./run_local.sh
```

The output will be generated in the "/output" folder by default and the global log files will be in "/output/log.out" and "/output/log.err". Each successful step will have its log files and output in a separate folder inside "/output".

## Description

This workflow has several parts:

- **Step 1**: extraction of structure from PDB and mutation of residue if requested.

- **Steps 2 (A-I)**: try to fix remaining PDB defects after extraction of the relevant structure and add mutations if needed. List of defects that the workflow attempts to fix: alternative locations, missing side chain atoms, SS bonds, clashing amides and errors regarding side chain chirality of residues. Note that fixing missing backbone atoms or missing residues is deactivated as this requires a PDB code and internet connection. Uncomment the relevant parts to use it.

- **Steps 3 - 7**: preparation of the system for minimization and MD with GROMACS.

- **Steps 8 - 10**: minimization

- **Steps 10 - 13**: NVT equilibration

- **Steps 13 - 16**: NPT equilibration

- **Steps 17 - 23**: launch several trajectories from equilibrated structure (see --n_trajs command line option). Concatenate those trajectories afterwards.

Note that re-launching the workflow will skip the previously successful steps if restart is True. 

