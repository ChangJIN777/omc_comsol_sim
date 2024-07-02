# omc-comsol
MATLAB scripts to automate COMSOL simulations for diamond optomechanical crystals. These simulations are based on the finite element method (FEM). Simulations can be run for optical and mechanical eigenmodes for nanobeam geometry, and mechanical bandstructure for unit cells

To learn more about OMC simulations in COMSOL, check the Loncar Group google drive, 02 Diamond Team/COMSOL-OMC tutorial.pptx

To edit the scripts - avoid editing the master branch unless it is to patch bugs that affects everyone using the scripts! Please pull a branch and edit within the branch.

## Getting started
These scripts can only be run by starting MATLAB with "COMSOL with MATLAB" - this will open an instance of COMSOL that is connected to the COMSOL server, allowing for COMSOL script commands to be interpreted

### Prerequisites
MATLAB (v2016 and above), COMSOL Multiphysics (v5.3a and above) <br />
For optical bandstructures: Python 3, anaconda or miniconda, meep, mpb, jupyter, matplotlib, pandas, mayavi

### Python and python packages installation instructions
MIT photonic bands simulation requires the installation of anaconda or miniconda:
https://www.anaconda.com/distribution/

Need to first install mpb and meep. Installation is really simple on Ubuntu (if anyone happens to use that). The website says Windows and Mac installation is a little bit more involved.
https://mpb.readthedocs.io/en/latest/Installation/
https://meep.readthedocs.io/en/latest/Installation/

### Insatallation in Ubuntu:
General dependencies: run in a terminal <br />
$ sudo apt-get install mpb h5utils <br />
$ sudo apt-get install meep h5utils <br />
<br />
Meep dependencies and virtual environment creation: <br />
This installs meeps dependencies into a python virtual environment to prevent conflicts with existing software, the environment must be activated in a terminal each time prior to running simulations. <br />
$ conda create -n mp -c conda-forge pymeep <br />
<br />
To activate the environent run: <br />
$ conda activate mp 

## Other dependencies:
mp is now a new environment with nothing but the necessary packages preinstalled, you will need to install everything necessary not included with anaconda while in the mp environment in a terminal while in the mp environment. <br />

$ conda install jupyter <br />
$ conda install matplotlib <br />
$ conda install pandas <br />
<br />
3D plotting functionality: Matplotlib is really good for figure quality 2D plotting but is not very well suited to 3D graphics the mayavi package is much better for this, follow the below installation instructions in a termial while in the mp environment and uncomment the initialization lines. <br />
$ pip install mayavi <br />
$ pip install PyQt5 <br />
<br />
Then to enable use within a jupyter notebook <br />
$ jupyter nbextension install --py mayavi --user <br />
$ jupyter nbextension enable --py mayavi --user <br />
<br />
### Relevant COMSOL settings
From the COMSOL GUI, open File > Preferences. In the dialog box:
- Select "LiveLink Connections", specify the correct MATLAB installation folder
- Select "Geometry", set Geometry representation in new models to "COMSOL kernel"

## Nanobeam simulations


### Overview of scripts
- **test_...** - scripts prefixed with test_ are front panel scripts with parameter specifications (in data structure P) and will run the RunNanobeamFEM function.
- **RunNanobeamFEM** - umbrella function containing sub-functions for each step in the COMSOL workflow
  - **CreateFileBase** - creates base file name given geometry parameters in P
  - **CreateNanobeamBlockTetGeom, CreateNanobeamGeom, CreateNanobeamGeom_asym** - generate hole dimensions for full geometries of block-tether nanobeam, symmetric nanobeam with holes, and asymmetric nanobeam with holes
  - **PlotDefectCells** - plots nanobeam geometry and hole defect parameters as function of mirror segment no.
  - **LoadMaterialParams** - loads predefined material properties into variable workspace
  - **RotateXtalTensor** - function to rotate elasticity and photoelastic tensors by given angle corresponding to crystal orientation to be simulated
  - **BuildNanobeamFEM** - builds geometry in COMSOL, generates selections for domains and boundaries
  - **SetupNanobeamFEM** - adds physics and defines boundary conditions in COMSOL for optical and/or mechanical simulations
  - **SolveNanobeamFEM** - adds eigensolver studies, meshes structure, and solves COMSOL model; plus postprocesses results to get resonance frequencies (optical and/or mechanical) and Q's (only for optical)
  - **CalcGOM** - calculates optomechanical coupling from optical and mechanical result dataset in COMSOL
  - **CalcStrCplSiV** - calculates strain coupling to SiV from strain results interpolated from mechanical simulations onto rectilinear grid, then processed in Matlab
  - **rotateStrainIntoSiV** - function in CalcStrCplSiV that does tensor rotation to rotate strain into SiV coordinate basis, with strain projected onto symmetry sectors if necessary
  - **PlotEy, PlotDispStr** - plots electric field (y-component), displacement and strain fields in nanobeam, using COMSOL functions
  - **PlotStrCplSiVXY, PlotStrCplSiVYZ** - plots strain coupling distribution, using postprocessing in Matlab

### Running the scripts
Modify parameters in test_... script, then run the script!

### Output from scripts
If the relevant options are enabled, a typical simulation outputs the following:
- COMSOL model file
- Matlab data structure (ds) containing parameter data structure (ds.P), optical and mechanical simulation parameters and results (ds.mfem, ds.ofem), coupling results (ds.cpl)
- Geometry plot generated by PlotDefectCells
- E-field plot generated by PlotEy
- Displacement and strain field plots generated by PlotDispStr
- Strain coupling distribution plots at max and at custom coordinate (if specified) generated by PlotStrCplSiVXY and PlotStrCplSiVYZ

## Mechanical unit cell simulations


### Overview of scripts
- **test_...** - scripts prefixed with test_ are front panel scripts with parameter specifications (in data structure P) and will run the SolveBands function.
  - **solveBands** - runs RunNanobeamBands for every specified symmetry combination, computes bandgaps, plots bandstructures
  - **RunNanobeamBands** - runs series of functions to set up, solve, and plot results from simulations
  - **Draw...** - creates geometry of unit cell in COMSOL
  - **RotateXtalTensor** - function to rotate elasticity and photoelastic tensors by given angle corresponding to crystal orientation to be simulated
  - **bndindex** - finds boundaries adjacent to given domain
  - **edgeindex** - finds edges adjacent to given boundary
  - **findGaps** - returns midgap frequencies and gap sizes for bands of each type of symmetry

### Running the scripts
Modify parameters in test_... script, then run the script!

### Output from scripts
If the relevant options are enabled, a typical simulation outputs the following:
- Matlab data structure containing bandstructure results for each symmetry, and complete bandgaps
- Bandstructure plots
- folder containing COMSOL model file, displacement and strain plots for unit cell at high symmetry points (Gamma- and X-points for 1D unit cells, plus M-point for 2D unit cells)

### Optical unit cell simulations
- **code/analyze_data.ipynb** - Jupyter Notebook used to visualize and plot data from already run simulations, can also run individual simulations.

- **code/analyze_sweep.py** - Capable of analyzing sweep data, currently only for convergence plot, will soon be implementing bandgap analysis of parameter sweeps.

- **code/photonic_crystal_waveguide.py** - Contains core code.

- **code/run_unit_cell.py** - Script used to initialize and run individual simulations or parameter sweeps, run this file in a terminal while in the mp environment.


## Acknowledgements
Original scripts written in COMSOL v3.5 from Michael Burek, who in turn got the scripts from Oskar Painter's group. Porting over to COMSOL v5.x done by Cleaven Chia. MIT Photonic Bands interface written by Graham Joe.
