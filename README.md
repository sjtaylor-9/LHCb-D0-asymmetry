# Measuring the Production Asymmetry between $D^0$ and  $\bar{D}^0$ mesons due to proton-proton collisions at the LHCb
This project is aimed to calculating the production asymmetry of the $D^0$ meson. This is done in differents regions of the phase space spanned by the transverse momentum ($p_T$) and pseudorapidity ($\eta$) of the $D^0$ meson. This is the current version (7th December) of the code and the predominent focus has been to implement a binned fit in order to increase run time speeds by 2-3 orders of magnitude. This repository will be updated throughout the next 5 months. The project was done in collaboration with Laxman Seelan.

In this repository there are the necessary tools in order to:
 - Make a selection of the events given a certain criteria,
 - Remove multiple candidates,
 - Perform a global fit on the data using either a binned or unbinned simultaneous extended maximum likelihood fit,
 - Produce a global model using the global fit fit parameters ,
 - Create a uniform binning across the phase space - ($p_T, \eta$) and individual $p_T$, $\eta$ binning schemes,
 - Perform a local fit in each of the phase space bins, by using a simultaneous extended maximum likelihood fit and plotting,
 - Process the results and output them with relevant figures,
 - Calculate production Asymmetry: Integrated and Average over the binning scheme bins,
 - Plot the asymmetry in the bins of ($p_T, \eta$),
 - Plot Asymmetry against $p_T$ and $\eta$.

Note that if this program is to be used with a different set of data, the file *selection_of_events.py* will need to be modified, and other modifications may be required as well. If phase space is changed, the plotting of Asymmetry and create and apply binning scheme python files may need to be changed.

**Warnings:**
While running the code be aware that any change to one of the scripts can lead to a malfunction. In addition, make sure that the directories that will be generated while running the program don't already exist. If they do exist beforehand, this program might not work as intended. ```model_fitting``` produces a segmentation violation, however this bug does not affect the code.

## How to download
In order to download this package you can use the following commands in your terminal:
```
mkdir ProductionAsymmetry
cd ProductionAsymmetry
git init
git clone git@github.com:sjtaylor-9/LHCb-D0-asymmetry.git
git pull origin main
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
A large amount of the scripts uses or is inspired by the code written by Camille Jarvis-Stiggants and Michael England during their MPhys project and Marc Oriol Pérez in his summer internship.


**Authors:** Sam Taylor (samuel.taylor-9@student.manchester.ac.uk) and Laxman Seelan (laxman.seelan@student.manchester.ac.uk)/ **Last modified:** 7th December 2023
