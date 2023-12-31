# SuMD
Python code to run Supervised Molecular Dynamics (SuMD) simulations using GROMACS
**N.B.**: a version using ACEMD3 engine can be found at github.com/molecularmodelingsection/SuMD

SuMD is a Python code that can be utilized to perform Supervised Molecular Dynamics simulations. The algorithm is deeply explained is the works of Sabbadin et. al (2014)<sup>1</sup>, Cuzzolin et. al (2016)<sup>2</sup> and Salmaso et al. (2017)<sup>3</sup>.
1. Sabbadin, D.; Moro, S. Supervised Molecular Dynamics (SuMD) as a Helpful Tool to Depict GPCR-Ligand Recognition Pathway in a Nanosecond Time Scale. J. Chem. Inf. Model. 2014, 54, 372–376.
2. Cuzzolin, A.; Sturlese, M.; Deganutti, G.; Salmaso, V.; Sabbadin, D.; Ciancetta, A.; Moro, S. Deciphering the Complexity of Ligand-Protein Recognition Pathways Using Supervised Molecular Dynamics (SuMD) Simulations. J. Chem. Inf. Model. 2016, 56, 687–705.
3. Salmaso, V.; Sturlese, M.; Cuzzolin, A.; Moro, S. Exploring Protein-Peptide Recognition Pathways Using a Supervised Molecular Dynamics Approach. Structure. 2017, 25,655–662.e2.

SuMD simulations can be used to enhance the sampling along a certain CV. Given a CV and a target value, the method will only accept the short MD cycles that approach the CV to its target value. For example, starting from a pre-equilibrated system in which the ligand is placed far enough from the protein. The general idea is to place the ligand at a distance, at least, bigger than the PME cut-off from any protein atom. The distance should be set even depending on the hydrodynamic properties of the ligand and the complexity of the binding
event. As rule of thumb, we suggest to place the ligand in a random conformation, in a range of 30-70 Å. The python script currently takes as inputs the GROMACS topology of the system compressed in a zip file, the initial structure as a .gro file and the configuration as a YAML file. 

The **configuration file** (here named **“input.yml”**) organized in three major sections containing information about (i) the supervision protocol, (ii) the MD simulation settings, and (iii) the analysis (pending). **An explanation on each required parameter to run the simulation is provided in the "input.yml" file.**  

To reconstitute the right Python virtual environment:
- **conda create -f environment.yml**  

If the environment reconstruction is too slow or raises conflicts, try to reconstruct the environment manually (i.e. installing the packages one by one in a new environment) starting with biobb_gromacs.

To run the code:
1. open a terminal within the directory of interest
2. activate the right conda environment (**conda activate biobb_sumd**)
3. run the code (**python /PATH/TO/SUMD_DIRECTORY/run_local.sh**)

Limitations of the method:

The criteria used by the method to accept or reject new sampling together with the short nature of the MD simulations used makes the study of complex and slow processes difficult with Supervised MD. Imagine a ligand adsorption that requires a large conformational change in the protein. The method will accept short MD that lead to a reduction in the main CV, i.e. distance of the ligand to the center of the pocket for example. But this is completely agnostic to any large conformational change in the protein... In general, any process that involves more than one slow CV will not be suitable to be studied with SuMD. It should work best for simple FE landscapes or if one uses a multiple walker strategy, see next point. 

For more complicated processes involving several slow degrees of freedom one has to resort to more refined strategies to identify a good CV from the MD simulation and then try to enhance sampling through the identified CV. Some examples are: spectral gap optimization of parameters (SGOOP), variational approach to conformational dynamics (VAC or VAMPnets), time-structure Independent Component Analysis (tICA), time-lagged Autoencoders (such as the RAVE method)... See works of Tiwary et al. and Frank Noe at al.

Currently we are using a single walker approach. Multiple instances can be launched but this would lead to an independent-multiple walker. A possible work-around is to launch several short SuMD trajectories in an iterative fashion. Where the more advanced trajectory along the CV is used to spawn the subsequent SuMD trajectories in the next iteration.

