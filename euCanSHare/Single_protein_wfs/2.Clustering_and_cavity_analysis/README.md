# Clustering and cavity analysis

![alt text](../../img/clust_cavity_analysis.png?raw=true)

## Quick installation and run

Go to workflow folder and install the conda environment (running in Nostrum's cluster use the already installed environments located in */shared/work/BiobbWorkflows/envs*)

```bash
conda env create -f environment.yml
conda activate biobb_sp_cavity_analysis
```

To install it in an HPC environment, do the same after loading the corresponding Conda or Miniconda module.

See different options for the worklow and modify any if needed:

```bash
vi input.yml
python biobb_clustering_cavity_analysis.py --help
```

Specially important are: the input files, the definition of the rmsd group and the filtering settings in the last two steps. If a trajectory file is given as input, also take a look to the total number of most pupulated clusters to analyze and the path to the GROMACS binary file in the global properties of the configuration file. Make sure the binary path specified and the module loaded in the run file agree between them.

See [biobb documentation](https://mmb.irbbarcelona.org/biobb/documentation/source) for additional properties not included in the YAML configuration file.

To run in an HPC environment adapt the run_HPC.sl and input_HPC.yml scripts and send a job to the slurm queue:

```bash
sbatch run_HPC.sl
```

To run locally, modify run_local.sh and input_local.yml if needed:

```bash
./run_local.sh
```

The output will be generatedin the "working_dir_path" folder selected in the corresponding YAML input. The global log files will be in "/working_dir_path/log.out" and "/working_dir_path/log.err". Each successful step will have its log files and output in a separate folder inside "/working_dir_path".

## Description

This workflow has several steps. The input for the workflow can be either (1) a trajectory and topology or (2) a path to a folder containing representative structures in pdb format from an external clustering. In the former case the workflow will cluster the trajectory to find representative structures (steps 0-2), in the later case the workflow will directly use the representative structures for the cavity analysis and filtering (steps 3-4). The command line arguments can be used to provide some inputs and settings that will be prioritized over those in the YAML configuration file.

**Step 0 (A-C)**: Creation of index files that define groups of atoms used during the clustering step.

- **Step 0A**: Create initial index file with standard groups from structure file (e.g. System, Protein, Protein-H, C-alpha, Backbone, MainChain...).

- **Step 0B**: Addition of 'RmsdGroup'. Corresponds to the group of atoms that will be used to fit the trajectory (unless -nofit option is used - see biobb docs for gmx_cluster) and to do the calculation of the RMSD. Check the documentation of gmx select to see all the possible atom selections. Some examples: 

    - Centers of mass of residues 1 to 5 and 10: "res_com of resnr 1 to 5 10"
    - All atoms of a residue LIG within 0.5 nm of a protein (with a custom name): '"Close to protein" resname LIG and within 0.5 of group "Protein"'
    - All protein residues that have at least one atom within 0.5 nm of a residue
  LIG: "group "Protein" and same residue as within 0.5 of resname LIG"

- **Step 0C**: Addition of 'OutputGroup' corresponds to atoms that will be included in the output representative structures. Check the documentation of gmx select to see all the possible atom selections.

- **Step 1**: Clustering of the trajectory. A trajectory (accepted formats: xtc, trr, cpt, gro, g96, pdb or tng) and a topology (accepted formats: tpr, gro, g96, pdb or brk) are read.

- **Step 2**: From the pdb file with all the centroids the most populated ones are extracted. The number of extracted centroids is defined in the num_clusters keyword of the YAML configuration file.

- **Step 3**: Cavity analysis of the centroid structures using Fpocket.

- **Step 4**: Filtering of the cavities found according to the criteria defined in the properties of this step. This includes filtering according to minimum and maximum values for the score, druggability score and volume of the pocket.

- **Step 5**: Filtering of the cavities according to the distance from their centers of mass to the center of mass of a selection defined in the properties of this step. The center of mass of a pocket is computed using the corresponding pqr file from Fpocket.

Note that re-launching the workflow will skip the previously successful steps if restart is True and the output folder is the same. 





