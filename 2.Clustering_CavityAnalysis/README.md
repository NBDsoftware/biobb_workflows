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

This workflow has several main sections, the workflow can be run until the end or until one of the sections (see --until command line option):

- **Section 1**: Clustering of trajectory. A trajectory (accepted formats: xtc, trr, cpt, gro, g96, pdb or tng) and a topology (accepted formats: tpr, gro, g96, pdb or brk) are read. The trajectory is clustered using the atoms defined in created index file (if --ndx command line option is used, see make_ndx step) or atoms defined in gmx_cluster step. From the pdb file with all the centroids the most populated ones are extracted. The number of extracted models is defined in the extract_models step.

- **Section 2**: Cavity analysis of a subset of the most populated models and a subsequent filtering of cavities.

Note that re-launching the workflow will skip the previously successful steps if restart is True. 

