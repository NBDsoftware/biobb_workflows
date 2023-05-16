# Cavity analysis

![alt text](../../img/cavity_analysis.png?raw=true)

## Quick installation and run

Go to workflow folder and install conda environment:

```bash
conda env create -f environment.yml
conda activate protein_protein_wf3
```

To install it in an HPC environment, do the same after loading the corresponding Conda or Miniconda module.

See options for worklow:

```bash
vi input_HPC.yml
python biobb_cavity_analysis.py --help
```

See [biobb documentation](https://mmb.irbbarcelona.org/biobb/documentation/source) for additional properties not included in the YAML configuration file.

To run in an HPC environment adapt the run_HPC.sl and input_HPC.yml scripts and send a job to the slurm queue:

```bash
sbatch run_HPC.sl
```

The output will be generated in the "working_dir_path" folder selected in the corresponding YAML input. The global log files will be in "/working_dir_path/log.out" and "/working_dir_path/log.err". Each step will have its own log files and output in a separate folder inside "/working_dir_path".

## Description

This workflow has several steps. The input for the workflow is a zip file containing pdb structures. Its path can be provided through the YAML configuration file or through the command line arguments. The YAML configuration file will contain the default settings and paths of the workflow. The command line arguments can be used to provide some inputs and settings that will be prioritized over those in the YAML configuration file.

- **Step 0**: extraction of all the poses from the input zip file. The next steps are done for each of the structures.

- **Step 1**: cavity analysis of the structures using Fpocket.

- **Step 2**: filtering of the cavities found according to the criteria defined in the properties of this step.




