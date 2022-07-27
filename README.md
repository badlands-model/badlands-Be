# 10Be Production & Transport from Source-to-Sink
## Credit: Carole Petit

Université Côte d'Azur, CNRS, Observatoire de la Côte d'Azur, IRD, Géoazur \
250 rue Albert Einstein, Sophia Antipolis 06560 Valbonne, France.

Contact: carole.petit@univ-cotedazur.fr
***

This repository contains a modified version of badlands code that simulates erosion, deposition and detrital 10Be transfer in source-to-sink sedimentary system.

The modified version of the code is available in the **badlands** folder, available pre- and post-processing functions to extract data from specific rock types and to compute their isotopic signatures are provided in the **badlands-companion** folder. Finally, a series of jupyter notebooks that illustrates how these new developments could be used are provided in the **Alps** folder.

## Installation

Easiest approach to install this version of badlands consists in downloading the repository from GitHub and use Anaconda. 

Once downloaded, go to the folder from the terminal and use the command below:

```bash
conda env create --file environment.yml
```

Once the environment has been succesfully installed, activate it:

```bash
conda activate badlandsBe
```

Then install the new developments by doing:

```bash
 cd badlands-companion
 python3 setup.py install --user
```

and 

```bash
 cd ../badlands
 python3 setup.py install --user
```

You are done!

## Testing the code on the Var River catchment (Southern French Alps) source-to-sink sedimentary system 

In a terminal, activate the `badlandsBe` conda environment:

```bash
conda activate badlandsBe
```

then launch the Jupyter environment:

```bash
jupyter notebook
```

In Jupyter, go to the `Alps` folder and run the provided notebooks.
