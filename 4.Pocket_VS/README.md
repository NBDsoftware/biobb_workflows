# Pocket VS

## Quick installation and run

Go to workflow folder and install conda environment:

```bash
conda env create -f environment.yml
conda activate pocket_vs
```

See options for worklow:

```bash
vi input.yml
python biobb_pocket_vs.py --help
```

See [biobb documentation](https://mmb.irbbarcelona.org/biobb/documentation/source) for additional properties not included in input.yml.

Modify run_local.sh if needed and launch:

```bash
./run_local.sh
```

The output will be generated in the "/output" folder by default and the global log files will be in "/output/log.out" and "/output/log.err". Each successful step will have its log files and output in a separate folder inside "/output".

## Description

This workflow has several main sections, the workflow can be run until the end or until one of the sections (see --until command line option):



Note that re-launching the workflow will skip the previously successful steps if restart is True. 

## Sequential run


