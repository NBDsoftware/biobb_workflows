# MD setup and run

## Quick installation and run

Go to workflow folder and install conda environment:

```bash
export KEY_MODELLER= "HERE YOUR MODELLER KEY"
conda env create -f environment.yml
conda activate eucanshare_wf1
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

This workflow has several main sections, the workflow can be run until the end or until one of the sections (see --until command line option):

- **Section 1**: extraction of structure from PDB and structure check. The input is either a pdb code or a pdb file.

- **Section 2**: try to fix remaining PDB defects after extraction of the relevant structure and add mutations if needed. List of defects that the section attempts to fix: alternative locations, missing backbone atoms or missing residues (provide PDB code), missing side chain atoms, SS bonds, clashing amides and errors regarding side chain chirality of residues. 

- **Section 3**: preparation of the system for minimization and MD with GROMACS and minimization.

- **Section 4**: NVT equilibration

- **Section 5**: NPT equilibration

- **Section 6**: launch several trajectories from equilibrated structure (see --n_trajs command line option). Concatenate those trajectories afterwards.

Note that re-launching the workflow will skip the previously successful steps if restart is True. 

## Sequential run

Activate environment:

```bash
conda activate eucanshare_wf1
```

It's a good idea to run the workflow sequentially to check the output of the different steps for a given PDB file or PDB code. Choose lower simulation times to debug. To retrieve the PDB and try to fix all defects:

```bash
python biobb_md_setup_mutation.py -i pdb:4LOO --config input.yml --until fix
```

Check that all possible PDB defects have been taken into account. A good place to start is the log file printed in step2I. Visualize fixed structure saved in last fix step. Check log files and output for different steps in this section (step2 A-I). Then launch minimization:

```bash
python biobb_md_setup_mutation.py -i pdb:4LOO --config input.yml --until min
```

If restart is True, steps 2A-2I will be skipped. The rest of the steps can also be launched sequentially. The workflow will automatically skip any successful step from previous calls. After the minimization we can launch the nvt equilibration:

```bash
python biobb_md_setup_mutation.py -i pdb:4LOO --config input.yml --until nvt
```

Then npt equilibration:

```bash
python biobb_md_setup_mutation.py -i pdb:4LOO --config input.yml --until npt
```

And the production run:

```bash
python biobb_md_setup_mutation.py -i pdb:4LOO --config input.yml --until all
```

