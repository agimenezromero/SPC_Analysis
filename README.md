# MagChain Analysis
Graphical interface in PyQt5 for the analysis of chain formation.

# Overview
The aim of this program is to simplify and automate the analysis of chain formation of superparamagnetic colloids under magnetic fields. With this graphical interface the analysis of different kinds of aggregation, including chain formation among others, can be done simply and easily, so anyone with basic or even null programming notions will be able to use all the software features. Although the software use is quite intuitive, a documentation has been written to assure its properly use, so we strongly recommend having a look on it, specially to understand the use of the > lateral_chains_finder function.


Table of contents
=================

<!--ts-->
   * [ChainAnalysis](#MagChain-Analysis)
   * [Table of contents](#table-of-contents)
   * [Installation](#installation)
   * [Usage](#usage)
     * [Loading files](#Loading-files)
       * [Pre analysis](#Pre-analysis)
       * [Post analysis](#Post-analysis)
      * [Analysing data](#Analysing-data)
        * [Cluster Finder](#Cluster-Finder)
        * [Linear Chains Finder](#Linear-Chains-Finder)
        * [Isolated Particle Finder](#Isolated-Particle-Finder)
        * [Radial Distribution Function](#Radial-Distribution-Function)
      * [Plotting analised data](#Ploting-analised-data)
      * [Customizing plots](#Customizing-plots)
   * [Examples](#Examples)
   * [Authors](#Authors)
   * [License](#License)
   * [Acknowledgments](#Acknowledgments)
<!--te-->

# Installation
1. Download ChainAnalysis.pyw file (no console python script) and Data folder and copy in the same location.
2. You must have the following python3 packages installed:
   - NumPy
   - SciPy
   - MdAnalysis
   - Matplotlib
   - Pandas
   - PyQt5
3. To run de program simply run the ChainAnalysis.pyw file by double clicking.

# Usage

## Loading Files

### Pre analysis
- In order to realize any analysis topology and trajectory files of the analysed system are required.
- The topology file is usually a single frame trajectory file (i.e. the initial frame) while trajectory file include the whole particle trajectories for all the frames.
- The admited formats are those admited by MdAnalysis to built the univers object, which include .pdb, .dcd, .xyz, etc. See [MdAnalysis documentation](https://www.mdanalysis.org/docs/) for further information.

### Post analysis
- After any analysis is done some .npy files are saved into self-created folders. You will have to load them to make the desired plots.
- Moreover, some .txt files are also created while any analysis executed, containing the results for the analysied data. In this way we let the user create more customized plots or even realise further analysis.
- More information about this files will be given in the following sections

## Analysing data
This software includes four different types of analysis: `Cluster Finder`, `Linear Chains Finder`, `Isolated Particles Finder` and `Radial Distribution Function (RDF)`. From each of these analysis we can get different information of the studied system.

### Cluster Finder
The Cluster Finder method combines MdAnalysis and SciPy libraries to find the number of aggregates of particles of selected type as well as the average length of this aggregates, that is to say the average number of particles in this aggregates. More specifically, this method makes use of the [scypi.cluster.hierarchy.linkage](https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html) method, which returns a hierarchical clustering encoded as a linkage matrix and then applyes the [scipy.cluster.hierarchy.fcluster](https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.cluster.hierarchy.fcluster.html) method to form the clusters from the linkage matrix.

Although this information may be irrelevant for some users it could be important for others. Now we'll explain in some lines how to use properly this analysis method in our software.

After loading the topology and trajectory files, we just have to set the parameters properly to perform our desired analysis. At the end of the analysis a .npy file is saved in the `Clusters` folder, inside the `Analysis files` folder, containing the number of selected particles that are found to be in a cluster over all the selected frames. The parameters needed are the following:

- `Radial cutoff`: Maximum interaction distance of the particles. Should be set accordingly to the RDF first minimum (as we are accounting for the firsts neighbours).
- `Initial frame`: Frame where the analysis will begin
- `Final frame`: Frame where thr analysis will finish
- `Time step`: Number of frames between consecutive frames in the analysis. In a good computer this method can take about 5 minutes to analyse 500 frames, so set the time step accordingly to this consideration.
- `Selection 1`: Type of particles to analyse. The selection syntax is the same used by [MdAnalysis](https://www.mdanalysis.org/docs/documentation_pages/selections.html) or [Visual Molecular Dynamics (VMD) software](https://www.ks.uiuc.edu/Research/vmd/current/ug/).
- `Filename`: The name of the .npy file that will be saved.

*The other paramaters are not used in this particular analysis

### Linear Chains Finder
The Linear Chains Finder method uses a usefull MDAnalysis geometric selection to decide whether selected particles are in a linear chain or not. The selection used is `cyzone`, with which a cylinder is formed arround every particle in the imput selection and search for neighbours. At the end of the analysis a .npy file is saved in the `Linear chains` folder, inside the `Analysis files` folder, containing the number of selected particles that are found to be in a chain over all the selected frames.

The parameters needed for this analysis are the following:

- `Radial cutoff`: It must be set to the value of the particles radius
- `Linear cutoff`: Maximum interaction distance of the particles in the z axis. Should be set accordingly to the RDF first minimum (as we are accounting for the firsts neighbours).
- `Initial frame`: Frame where the analysis will begin
- `Final frame`: Frame where thr analysis will finish
- `Time step`: Number of frames between consecutive frames in the analysis. In a good computer this method can take about 25 minutes to analyse 500 frames, so set the time step accordingly to this consideration.
- `Selection 1`: Type of particles to analyse. 
- `Filename`: The name of the .npy file that will be saved.

*The other paramaters are not used in this particular analysis

### Isolated Particle Finder
The Isolated Particle Finder method uses a spherical selection to determine if selected particles are isolated or not by counting the number of neighbours at a certain distance. With this information, it is clear that we can also know the number of general aggregated particles, being linear aggregates or other structures, as the laterally aggregated particles, by comparison with the results of the same temporal analysis done by the Linear Chains Finder (as if the particles are not linearly aggregated, they are laterally aggregated). At the end of the analysis two .npy files are saved in the `Isolated particles` and `General aggregates` folders, inside the `Analysis files` folder, containing the number of isolated particles and the number of particles in general aggregates, respectively over all the selected frames.

The parameters needed for this analysis are the following:

- `Radial cutoff`: Maximum interaction distance of the particles. Should be set accordingly to the RDF first minimum (as we are accounting for the firsts neighbours).
- `Initial frame`: Frame where the analysis will begin
- `Final frame`: Frame where thr analysis will finish
- `Time step`: Number of frames between consecutive frames in the analysis. In a good computer this method can take about 25 minutes to analyse 500 frames, so set the time step accordingly to this consideration.
- `Selection 1`: Type of particles to analyse. 
- `Filename`: The name of the .npy file that will be saved.

*The other paramaters are not used in this particular analysis

### Radial Distribution Function
The RDF method simply computes a radial distribution function using [MDAnalysis.analysis.rdf.InterRDF](https://www.mdanalysis.org/docs/documentation_pages/analysis/rdf.html).

The parameters needed for this analysis are the following:

- `Initial frame`: Frame where the analysis will begin
- `Final frame`: Frame where thr analysis will finish
- `Time step`: Number of frames between consecutive frames in the analysis.
- `r max`: Maximum radial distance
- `dr`: Radial distance step
- `Selection 1`: Type of group 1 particles.
- `Selection 2`: Type of group 2 particles.
- `Filename`: The name of the .npy file that will be saved.

## Plotting analised data
Once desired analysis are done is time to visualise the results. For this purpose a dropdown menu has been implemented in the graphical interface in order to let the user choose the results to plot. Nevertheless the .npy file saved by each analysis method, previously mentioned, must be loaded using the corresponding load button.

Once the data file is loaded and the corresponding option of the dropdown menu selected the only thing left to do is to push the plot button. Then, if all steps have been correctly made, there will appeaar a matplotlib.pyplot window with the analysis results.

In the following section we will focus in how to customize this plots.

## Customizing plots

