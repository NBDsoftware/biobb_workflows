# Docking (HTVS)

## Quick installation and run

Go to workflow folder and install conda environment:

```bash
conda env create -f environment.yml
conda activate eucanshare_wf3
```

See options for worklow:

```bash
vi input.yml
python biobb_docking_htvs.py --help
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

- **Section 1**: pocket selection and receptor preparation for docking.

- **Section 2**: docking of all the ligands provided through --lig-lib command line argument.

Note that re-launching the workflow will skip the previously successful steps if restart is True. 

## Sequential run

```bash
conda activate eucanshare_wf3
```

A summary of the models and pockets available from the input folder given in 'input_pockets_zip' is printed if an unexisting pocket is chosen.

