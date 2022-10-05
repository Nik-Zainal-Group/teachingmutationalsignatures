# teachingmutationalsignatures


## Introduction

Material for teaching about  mutational signatures extraction and fit algorithms using the ```signature.tools.lib``` R package.

## Versions

1.0.0

- First six teaching sessions
- Functions for loading simulated datasets and testing accuracy of results


## Installation

To use this package you need to install the ```signature.tools.lib``` R package version 2.2.0 or higher, which is available on github [here](https://github.com/Nik-Zainal-Group/signature.tools.lib).

To install the ```teachingmutationalsignatures``` R package you need to have the ```devtools``` package installed.
Open R and set the working directory to the main package folder, and use the command ```devtools:install()```.

## Package functions and teaching sessions

The teaching material is organised into teaching sessions. Each session is organised into one single R file within the
```teaching_material``` folder, and typically will cover the analysis of a simulated dataset.

The ```teachingmutationalsignatures``` R package provides functions to access the simulated data and test how good
the obtained results are.

- **```fetchData(...)```**: returns the mutational catalogues of the samples from a given simulated dataset. 
The only parameter is the dataset name, such as ```SD001```, ```SD002```,... 
- **```checkPerformanceSignatures(...)```**: check a matrix of estimated signatures against the true signatures
of a given simulated dataset. Multiple estimated signatures can be compared using the ```checkPerformanceSignaturesList``` function.
- **```checkPerformanceExposures(...)```**: check a matrix of estimated exposures against the true exposures
of a given simulated dataset. Multiple estimated exposures can be compared using the ```checkPerformanceExposuresList``` function.
