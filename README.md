# Paired Image X-fel Integrator

Integrate pump probe X-FEL data using nXDS (http://nxds.mpimf-heidelberg.mpg.de/). PIXI accepts sets of images acquired with crystals in the same position and integrates them on a per image basis using the same reflection profiles and geometric parameters for equivalent images. 


# Usage

pixi [OPTION] [DATASET 1 IMAGE 1] [DATASET 2 IMAGE 2] ...

PIXI will integrate an arbitrary number of datasets with parameters learned from the first dataset. Specify the first image in each dataset. Images must be enumerated with fixed with, zero padded integers __immediately__ prior to the file extension (/path/to/image/prefix00001.cbf). 

*   -x, --xdsin
XDS.INP file with default parameters to use during integration. Some parameters are overridden during the process, including but not limitted to those relating to spot profiling and parameter refinement. 

*   -g, --generate-xdsin
Optionally, use the generate_XDS.INP to automatically generate sane parameters for nXDS. This should be treated with extreme skepticism. nXDS and XDS share a lot of options but not all. Furthermore, this generate_XDS.INP can only do as well as your input file's header. 

*   -h, --help
Print a description of the arguments.

