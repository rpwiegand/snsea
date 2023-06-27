# Simple Novelty Search Source Code

This code supports the experiments run in the following paper:

> Wiegand, R. P.  (2023).  "Preliminary Analysis of Simple Novelty Search".  *Evolutionary Computation*.

## What Is In This Project
The primary purpose of this repository is to provide the Python 3 source code for implementing the *Simple Novelty Search Evolutionary Algorithm*, as well as the *Population Based Simple Novelty Search Evolutionary Algorithm* discussed in the paper.  I've also tried to include as many of the `ini` files needed to replicate experiments in multiple papers as possible.

| Source File          | Description                                             |
| -------------------- | ------------------------------------------------------- |
| configReader.py      | A class to parse command-line and `ini` file parameters |
| pseudoRobot.py       | A class to to implement genotype-behavior mapping       |
| snseaBase.py         | The base class for all SNSEAs                           |
| snseaNoPopulation.py | The base population-less algorithm                      |
| snseaPopulation.py   | The algorithms that uses a population                   |

The SNSEA itself maintains only an individual and an archive.  New individuals are spawned from a single parent, and added to the archive if they meet the sparesness criterion.  A new parent is drawn randomly from the archive.  In the *population* version, a population is maintained *in addition to* the archive.  Parents are drawn from the population to produce *lambda* children and the *mu* "most sparse" become parents.  Sparseness can be computed based on the archive or the existing population (this is a parameter), and we also have the ability to select a percentage of random individuals from the archive.  The algorithms are described in detail in the paper.  

All metrics discussed in the paper are implemented in `snseaBase.py`.

In addition, the `ExperimentParams` directory contains several studies directories.  Each of those contains `ini` files with the parameters used to run that specific parameter set for that study.  The studies are:

| Study           | Description                                                                                        |
| --------------- | -------------------------------------------------------------------------------------------------- |
| ECJ-NoPop       | Experiments for the ECJ paper involving just the SNSEA                                             |
| ECJ-PopPop      | ECJ experiments with a population where survival is based on sparseness compared to the population |
| ECJ-RandSample  | ECJ experiments for archive visualization                                                          |
| flairs-paper    | Experiments for the FLAIRS 2020 paper                                                              |
| Misc            | Other ini files that I kept                                                                        |


## Executing the Code
If you want to run without a population, you will executate the `snseaNoPopulation.py` file.  Otherwise, you will run the `snseaPopulation.py` file.  The program takes as a parameter the `ini` file with the experiment parameters to use.  For example:

```
python3 snseaNoPopulation.py ExperimentParams/ECJ-NoPop/hamming-500-n10.ini
```


## Modules Needed
This code relies on some basic Python modules, which must be installed.  These include:

* ImportLib
* Matplotlib
* Numpy
* Pandas
* SafeConfigParser
* SciPy


