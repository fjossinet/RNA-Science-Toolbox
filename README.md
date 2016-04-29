RNA Science Toolbox
===================

The RNA Science Toolbox provides a Python API (PyRNA) to do RNA science. PyRNA allows you to:

* parse RNA data from "classical" file formats (PDB, CT, FASTA, VIENNA,...) and convert them into easy-to-use and easy-to-analyse data structures:
    * [Pandas Series and DataFrames](http://pandas.pydata.org/pandas-docs/stable/dsintro.html)
    * "PyRNA objects" (as defined in the module pyrna.features).
* compute RNA data from RNA algorithms (see below for details) and convert them into Pandas data structures and PyRNA objects,
* recover RNA data from public databases ([PDB](http://www.rcsb.org/pdb/home/home.do), [RFAM](http://rfam.sanger.ac.uk),...) and convert them into Pandas data structures and PyRNA objects,
* deploy some functionalities as REST Web services.

This project is related to the [DockeRNA project](https://github.com/fjossinet/DockeRNA) which provides Docker images containing the RNA algorithms you may need.

You can [follow this project on twitter](https://twitter.com/RnaSciToolbox) to get updates as they happen.

Installation
------------

To use the RNA Science Toolbox, you will need to go through several steps. But don't be afraid, each step is really easy to follow. We do suppose that you are using either MacOSX or Linux.

### Python environment

You need at first to have a Python distribution installed on your computer. If you don't have one, we recommend you a distribution like [Anaconda](https://www.continuum.io/why-anaconda).

### Fabric

You also need the tool [Fabric](http://www.fabfile.org). If you're using the [Anaconda distribution](https://www.continuum.io/why-anaconda), you can get it by typing:

    conda install fabric

### Python dependencies

Once done, download the RNA Science Toolbox and go into its directory. To install its Python dependencies, you can  use either the package manager conda (from the [Anaconda distribution](https://www.continuum.io/why-anaconda)) or pip. To use conda, type:

    fab

To use pip, type:

    fab install:manager=pip

### Docker

To install the RNA algorithms, you need first to get the tool Docker. You will find all the details [here](https://docs.docker.com/engine/installation/).

### RNA algorithms

Each Docker image available contains several algorithms:

 * [fjossinet/assemble2](https://hub.docker.com/r/fjossinet/assemble2/): provides [RNAVIEW](http://ndbserver.rutgers.edu/ndbmodule/services/download/rnaview.html), [Vienna RNA package](https://www.tbi.univie.ac.at/RNA/), [foldalign](http://rth.dk/resources/foldalign/), [LocARNA](http://rna.informatik.uni-freiburg.de/LocARNA/)
 * [fjossinet/rnaseq](https://hub.docker.com/r/fjossinet/rnaseq/): provides [SAMtools](http://samtools.sourceforge.net), [Tophat2](https://ccb.jhu.edu/software/tophat/), [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

Go to their Web page for more details about their installation and usage.

### Configure your PATH

In the configuration file of your shell (.bashrc, .zshrc,...), add the folowing lines:

    export TOOLBOX=THE_PATH_TO_YOUR_RNA_SCIENCE_TOOLBOX
    export PYTHONPATH=$PYTHONPATH:$TOOLBOX
    export PATH=$PATH:$TOOLBOX/files/scripts/python:$PATH

Restart your shell and type:

    pyrna_tests.py

Your RNA Science Toolbox is fully configured if you get something like:

<pre>
Recovering entry 1EHZ from Protein Databank...

## 3D annotation ##

List of base-pairs computed with RNAVIEW:

edge1 edge2 orientation  pos1  pos2
0      (     )           c     1    72
1      (     )           c     2    71
2      (     )           c     3    70
[...]
</pre>
