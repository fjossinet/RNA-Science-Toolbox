RNA Science Toolbox -- A virtual environment to do RNA science
===============================================================

The RNA Science Toolbox provides a Python API (PyRNA) to do RNA science. This virtual environment allows you to:
* write your own Python scripts based on PyRNA,
* deploy a Web server providing documentations and Web services.

The PyRNA API allows you to:

* parse RNA data from "classical" file formats (PDB, CT, FASTA, VIENNA,...) and convert them into easy-to-use and easy-to-analyse data structures:
    * [Pandas Series and DataFrames](http://pandas.pydata.org/pandas-docs/stable/dsintro.html)
    * "PyRNA objects": defined in the module pyrna.features.
* compute RNA data from RNA algorithms (see list below) and convert them into Pandas data structures and PyRNA objects,
* recover RNA data from public databases ([PDB](http://www.rcsb.org/pdb/home/home.do), [RFAM](http://rfam.sanger.ac.uk),...) and convert them into Pandas data structures and PyRNA objects,
* deploy some functionalities as REST Web services.

In its current state, PyRNA is able to handle an ever-increasing number of (mainly) RNA algorithms like:

* [The RNA Vienna Package](http://goo.gl/7frDgF)
* [CONTRAfold](http://goo.gl/4BCI7)
* [Infernal](http://goo.gl/SxLHJO)
* [Mlocarna](http://goo.gl/AIGKrl)
* [RNAMotif](http://goo.gl/MDdOQ2)
* [RNAVIEW](http://goo.gl/c5o19v)
* ....

You can [follow this project on twitter](https://twitter.com/RnaSciToolbox) to get updates as they happen.

#Quick Start

You can find all the details [here](http://jossinetlab.github.io/RNA-Science-Toolbox/)
