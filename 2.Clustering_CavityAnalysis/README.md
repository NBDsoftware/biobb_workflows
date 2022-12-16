# Clustering and cavity analysis

## Quick installation and run

Go to workflow folder and install conda environment:

```bash
conda env create -f environment.yml
conda activate eucanshare_wf2
```

See options for worklow:

```bash
vi input.yml
python biobb_clustering_cavity_analysis.py --help
```

See [biobb documentation](https://mmb.irbbarcelona.org/biobb/documentation/source) for additional properties not included in input.yml.

Modify run_local.sh if needed and launch:

```bash
./run_local.sh
```

The output will be generated in the "/output" folder by default and the global log files will be in "/output/log.out" and "/output/log.err". Each successful step will have its log files and output in a separate folder inside "/output".

## Description

This workflow has several main sections, the workflow can be run until the end or until one of the sections (see --until command line option). The input can be either a trajectory and topology or a path to a folder containing representative structures in pdb format. In the former case the workflow will cluster the trajectory to find representative structures, in the later case the workflow will directly use the representative structures cavity analysis.

- **Section 1**: Creation of index files. 'FitGroup' corresponds to the group of atoms that will be used in the least squares fit and RMSD calculation. 'OutputGroup' corresponds to atoms that will be included in the output representative structures.

- **Section 2**: Clustering of trajectory. A trajectory (accepted formats: xtc, trr, cpt, gro, g96, pdb or tng) and a topology (accepted formats: tpr, gro, g96, pdb or brk) are read. From the pdb file with all the centroids the most populated ones are extracted. The number of extracted models is defined in the extract_models step.

- **Section 3**: Cavity analysis of a subset of the most populated models and a subsequent filtering of cavities.

Note that re-launching the workflow will skip the previously successful steps if restart is True. 

## Sequential run

```bash
conda activate eucanshare_wf2
```

It's a good idea to run the workflow sequentially to check the output of the different steps for a given trajectory and structure file. Modify 'input_structure_path' and 'input_traj_path' from the input.yml or give them through the corresponding command line arguments. The atom selection to align and compute the RMSD is the same and is defined by the 'FitGroup'. To align and compute the RMSD using different atom selections do the clustering externally (using ttclust/mdtraj or cpptraj). The atom selection that will be used when writing the output representative structures is defined by 'OutputGroup'. 

If no index file is needed, (if we want to align all the backbone atoms of the protein for example) then just use the corresponding groups in 'fit_selection' and 'output_selection'. If the index file is already created, then provide its path with '--ndx-file' command line option. If the index file should be created, then provide a suitable atom selection string in 'selection' from the input.yml and the corresponding group in 'fit_selection' and/or 'output_selection'.

In this example, we set the structure and trajectory paths to the output of the previous workflow.

- input_structure_path : /path/to/eucanshare_wfs/1.MD_setup_mutation/output/step23_dry/imaged_structure.gro
- input_traj_path : /path/to/eucanshare_wfs/1.MD_setup_mutation/output/step24_trjcat/all_trajectories.trr

Launch index file creation:

```bash
python biobb_clustering_cavity_analysis.py --config input.yml --until ndx
```

Then, cluster the trajectory:

```bash
python biobb_clustering_cavity_analysis.py --config input.yml --until cluster
```

After seeing the results, adjust the 'cutoff' and 'method' of the clustering if needed. Adjust the number of models to extract if needed. To re-launch this section simply erase the corresponding step folders and launch again (make sure 'restart' is True in the input.yml). Finally the cavity analysis and filtering:

```bash
python biobb_clustering_cavity_analysis.py --config input.yml --until all
```

Repeat the filtering with different values removing the filtering step folder and launching again.





