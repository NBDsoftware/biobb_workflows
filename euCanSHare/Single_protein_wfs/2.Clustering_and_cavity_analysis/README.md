# Clustering and cavity analysis

![alt text](../../img/clust_cavity_analysis.png?raw=true)

## Quick installation and run

Go to workflow folder and install the conda environment (running in Nostrum's cluster use the already installed environments located in */shared/work/BiobbWorkflows/envs*)

```bash
conda env create -f environment.yml
conda activate biobb_sp_cavity_analysis
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
Specially important are: the filtering settings in the last two steps. If a trajectory file is given as input, also take a look to the total number of most populated clusters to analyze, the definition of the rmsd group to do the clustering and the path to the GROMACS binary file in the global properties of the configuration file. Make sure the binary path specified and the module loaded in the run file agree between them. If the "--prepare_traj" flag is used, take a look at the selection group definition to trim the trajectory and the start, end and dt values to sub-sample the trajectory if needed.

### Command line arguments

The command line arguments can be used to provide some inputs and settings that will be prioritized over those in the YAML configuration file.

```bash
python biobb_clustering_cavity_analysis.py --help
```

Specially important are: the configuration file path, the path to the external clustering results with representative structures in pdb format or the path to the trajectory and topology files. The trajectory file can be in any Amber or Gromacs compatible format  (xtc, trr, cpt, gro, g96, pdb or tng). The topology file should be a pdb file if the trajectory is in Amber format (mdcrd, crd, cdf, netcdf, restart, ncrestart, dcd, charmm, cor, mol2, trr, binpos, xtc, cif, arc, sqm, sdf, conflib) or any Gromacs-compatible format (tpr, gro, g96, pdb or brk) if the trajectory is in a Gromacs format.

## Description

This workflow has several steps. The input for the workflow can be either (1) a trajectory and topology or (2) a folder containing representative structures from an external clustering. In the former case the workflow will cluster the trajectory to find representative structures, in the later case the workflow will directly use the representative structures for the cavity analysis and filtering.

- **Step 0 **: Conversion of trajectory from Amber to Gromacs xtc format. This step is only activated if the "--prepare_traj" flag is used and the provided trajectory is in Amber format.

- **Step 1 (A-B)**: Creation of index file that will be used to select some atoms from the trajectory. This step can be used to delete waters, ions or any exotic atom from the trajectory. As the subsequent clutering step using gmx cluster might give problems when dealing with these. The step is only activated if the "--prepare_traj" flag is used.

- **Step 2 (A-B)**: Extraction of selected atoms from both the trajectory and the topology. This step is only activated if the "--prepare_traj" flag is used.

- **Step 3 (A-C)**: Creation of index files that define groups of atoms used during the clustering step.

- **Step 3A**: Create initial index file with standard groups from structure file (e.g. System, Protein, Protein-H, C-alpha, Backbone, MainChain...).

- **Step 3B**: Addition of 'RmsdGroup'. Corresponds to the group of atoms that will be used to fit the trajectory (unless -nofit option is used - see biobb docs for gmx_cluster) and to do the calculation of the RMSD. Check the documentation of gmx select to see all the possible atom selections. Some examples: 

    - Centers of mass of residues 1 to 5 and 10: "res_com of resnr 1 to 5 10"
    - All atoms of a residue LIG within 0.5 nm of a protein (with a custom name): '"Close to protein" resname LIG and within 0.5 of group "Protein"'
    - All protein residues that have at least one atom within 0.5 nm of a residue
  LIG: "group "Protein" and same residue as within 0.5 of resname LIG"

- **Step 3C**: Addition of 'OutputGroup' corresponds to atoms that will be included in the output representative structures. Check the documentation of gmx select to see all the possible atom selections.

- **Step 4**: Clustering of the trajectory.

- **Step 5**: From the pdb file with all the centroids the most populated ones are extracted. The number of extracted centroids is defined in the num_clusters keyword of the YAML configuration file.

- **Step 6**: Cavity analysis of the centroid structures using Fpocket.

- **Step 7**: Filtering of the cavities found according to the criteria defined in the properties of this step. This includes filtering according to minimum and maximum values for the score, druggability score and volume of the pocket.

- **Step 8**: Filtering of the cavities according to the distance from their centers of mass to the center of mass of a selection defined in the properties of this step. The center of mass of a pocket is computed using the corresponding pqr file from Fpocket.

Note that re-launching the workflow will skip the previously successful steps if restart is True and the output folder is the same. 





