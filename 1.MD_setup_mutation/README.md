# MD setup and run

## Quick installation and run

Go to workflow folder and install conda environment:

```bash
export KEY_MODELLER=MODELIRANJE
conda env create -f environment.yml
conda activate eucanshare_wf1
```

See options for worklow:

```bash
vi input.yml
python biobb_md_setup_mutation.py --help
```

See [biobb documentation](https://mmb.irbbarcelona.org/biobb/documentation/source) for additional properties not included in input.yml.

Modify run_local.sh if needed and launch:

```bash
./run_local.sh
```

The output will be generated in the "/output" folder by default and the global log files will be in "/output/log.out" and "/output/log.err". Each successful step will have its log files and output in a separate folder inside "/output".

## Description

This workflow has several main sections, the workflow can be run until the end or until one of the sections (see --until command line option):

- **Section 1**: extraction of structure from PDB and structure check. 

- **Section 2**: try to fix remaining PDB defects after extraction of the relevant structure and add mutations if needed. List of defects that the section attempts to fix: alternative locations, missing backbone atoms or missing residues (provide PDB code), missing side chain atoms, SS bonds, clashing amides and errors regarding side chain chirality of residues. 

- **Section 3**: preparation of the system for minimization and MD with GROMACS and minimization.

- **Section 4**: NVT equilibration

- **Section 5**: NPT equilibration

- **Section 6**: launch several trajectories from equilibrated structure (see --n_trajs command line option). Concatenate those trajectories afterwards.

