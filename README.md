## Installation

### Within the Nextgen Frameowork
The first thing to do is install Nextgen. The next thing to do is ...  

### In stand alone version
Install any dependencies (see below). Then open the `run\_test.ipynb` and shift-click through the notebook.

----

# General Data Assimilation with BMI

Streamflow predictions derived from a hydrologic model are hampered by manysources of errors, including uncertainties in meteorological inputs and inadequaterepresentation of natural systems’ physical processes. One way to reduce the impact of theseuncertainties on models’ accuracy and improve streamflow predictions is to implement DataAssimilation (DA) methods. DA methods use streamflow observations to update models’ statevariables to improve streamflow predictions. This study aimed to build and demonstrate DA asa module ready for use in the Next Generation Water Resources Modeling Framework(Nextgen). All components of this module conform to the Basic Model Interface (BMI)standards to be compatible with the Nextgen framework. This includes a BMI-enabledpython-module that uses the Ensemble Kalman Filter (EnKF) to derive state transitionfunctions for other Nextgen modules. An additional BMI-enabled python-module was alsocreated to collect USGS streamflow observations to assimilate streamflow predictions. Thedeveloped EnKF and USGS modules were tested using streamflow simulations from theBMI-enabled Conceptual Functional Equivalent (CFE) version 2.1 of the National WaterModel (NWM). The existing BMI-enabled CFE was used to create a new version of theBMI-enabled CFE with added perturbations to produce an ensemble of streamflowsimulations required to run the EnKF method. The DA module was tested on synthetic dataand a real-world example at a USGS gauged basin. The outcomes showed that implementingthe EnKF DA method improved the CFE streamflow predictions. The difference between theassimilated streamflow predictions and USGS streamflow observations was smaller than thedifference between original CFE simulations and USGS observations. When compared to theUSGS streamflow observations, the assimilated CFE streamflow predictions improvedcorrelation coefficient (r), Kling-Gupta Efficiency (KGE), Nash-Sutcliffe Efficiency (NSE),and percent bias (PBIAS) statistics.

## Dependencies

There are two paths to run this software, and they have separate dependencies.  
There will be an environment.yml and a environment.txt to run this code in stand alone version with Python.

## Installation

Detailed instructions on how to install, configure, and get the project running.
This should be frequently tested to ensure reliability. Alternatively, link to
a separate [INSTALL](INSTALL.md) document.


## Usage

There will be a notebook which demonstrates how to run the code outside of the Nextgen Framework.  
There will be some code explaining how to run this within the Nextgen Framework

## How to test the software

If the software will include automated tests.

## Known issues

This repository has not code in it yet.

## Getting help

If you have questions, concerns, bug reports, etc, please file an issue in this repository's Issue Tracker.

## Getting involved


Please see our contribution instructions [CONTRIBUTING](CONTRIBUTING.md).


----

## Open source licensing info
1. [TERMS](TERMS.md)
2. [LICENSE](LICENSE)


----

## Credits and references

1. Deardorff E., A. Modaresi Rad, et al. (2022). National Water Center Innovators Program - Summer Institute, CUAHSI Technical Report, HydroShare, http://www.hydroshare.org/resource/096e7badabb44c9f8c29751098f83afa
2. Wolkeba F. T., K. O. Ekpetere, M. S. Abualqumboz , and Z. J. Butler (2022). Data assimilation of USGS streamflow in the Nextgen Framework
3. Abualqumboz et al. 2022. “Data Assimilation of USGS Streamflow Data in the Nextgen Framework”. AGU Fall Meeting




# Some old figures

![](./doc/readme_figs/nwm_vs_bmi_da.png)
This figure describes the difference between the naive NWM data assimilation method and this BMI-based general data assimilation framework.
![](./doc/readme_figs/buckets.png)
![](./doc/readme_figs/logic_soil.png)
![](./doc/readme_figs/soil_and_giuh.png)
