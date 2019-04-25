# ChainAnalysis
Graphical interface in PyQt5 for the analysis of chain formation.

# Overview
The aim of this program is to simplify and automate the analysis of chain formation of superparamagnetic colloids under magnetic fields. With this graphical interface the analysis of different kinds of aggregation, including chain formation among others, can be done simply and easily, so anyone with basic or even null programming notions will be able to use all the software features. Although the software use is quite intuitive, a documentation has been written to assure its properly use, so we strongly recommend having a look on it, specially to understand the use of the > lateral_chains_finder function.


Table of contents
=================

<!--ts-->
   * [ChainAnalysis](#ChaynAnalysis)
   * [Table of contents](#table-of-contents)
   * [Installation](#installation)
   * [Usage](#usage)
     * [Loading files](#Loading-files)
       * [Pre analysis](#Pre-analysis)
       * [Post analysis](#Post-analysis)
      * [Analysing data](#Analysing-data)
      * [Ploting analized data](#Ploting-analized-data)
      * [Customizing plots](#Customizing-plots)
   * [Tests](#tests)
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
- In order to realize any analysis topology and trajectory files of the analysed sistem are required.
- The topology file is usually a single frame trajectory file (i.e the initial frame) while trajectory file include the trajectories for all the frames.
- The admited formats are those admited by MdAnalysis to built the univers object, which include .pdb, .dcd, .xyz, etc. See [MdAnalysis documentation](https://www.mdanalysis.org/docs/) for further information.

