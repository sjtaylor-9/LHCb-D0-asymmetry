# D0 production asymmetry calculator using K-pi
This project is aimed to calculating the production asymmetry of the D0 meson. This is done in differents regions of the phase space spanned by the transverse momentum and pseudorapidity of the D0 meson.

In this repository there are the necessary tools in order to:
 - Make a selection of the events given a certain criteria
 - Remove multiple candidates
 - Perform a global fit on the data using a simultaneous fit
 - Fit both using a binned and an unbinned approach
 - Create a uniform binning across the phase space
 - Perform a local fit in each of the phase space regions
 - Process the results and output them with relevant figures

Note that when performing the local fit some of the parameters have been fixed, using the values obtained in the global fit, in order to ensure convergence. Also note that if this program is to be used with a different set of data, the file *selection_of_events.py* will need to be modified, and other modifications may be required as well.

**Warnings:**
While running the code be aware that any change to one of the scripts can lead to a malfunction. In addition, make sure that the directories that will be generated while running the program don't already exist. If they do exist beforehand, this program might not work as intended.

## How to download
In order to download this package you can use the following commands in your terminal:
```
mkdir ProductionAsymmetry
cd ProductionAsymmetry
git init
git remote add -f origin https://github.com/mark1ry/D0_production_asymmetry
git pull origin master
```

## How to use
The different scripts can be run individually (note that a different set of arguments is required for each), or as a whole using the bash script *main.sh*.
In order to use *main.sh* 4 arguments are required. These are:
- The path where the output should be written
- The year the data to be used was taken [16, 17, 18]
- The size or amount of data to be used [small, medium, large, 1, 2, 3, 4, 5, 6, 7, 8]
- Whether a binned fit should be performed (otherwise unbinned fit) [y, Y, n, N]

Here is an example of how to call *main.sh*:
```
bash main.sh example 18 large y
```
This should produce the same output as shown in the folder *example* (still to be implemented).
## Credits
A large amount of the scripts uses or is inspired by the code written by Camille Jarvis-Stiggants and Michael England during their MPhys project.


**Author:** Marc Oriol Pérez (marc.oriolperez@student.manchester.ac.uk) / **Last modified:** 16th September 2023
