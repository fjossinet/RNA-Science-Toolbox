RNA Science Toolbox
===================

The RNA Science Toolbox provides a Python API (PyRNA) to do RNA science. PyRNA allows you to:

* parse RNA data from "classical" file formats (PDB, CT, FASTA, VIENNA,...) and convert them into easy-to-use and easy-to-analyse data structures:
    * [Pandas Series and DataFrames](http://pandas.pydata.org/pandas-docs/stable/dsintro.html)
    * "PyRNA objects": defined in the module pyrna.features.
* compute RNA data from RNA algorithms (see list below) and convert them into Pandas data structures and PyRNA objects,
* recover RNA data from public databases ([PDB](http://www.rcsb.org/pdb/home/home.do), [RFAM](http://rfam.sanger.ac.uk),...) and convert them into Pandas data structures and PyRNA objects,
* deploy some functionalities as REST Web services.

To be able to use RNA algorithms, you need to install them. This can be easily done by installing Docker images providing these algorithms fully installed and configured.

You can [follow this project on twitter](https://twitter.com/RnaSciToolbox) to get updates as they happen.

Installation
------------

To install the RNA Science Toolbox, you need to go through several steps. But don't be afraid, each step is really easy to follow. This installation workflow targets computers running under MacOSX or Linux.

### Python environment

You need at first to have a Python distribution installed on your computer. If you don't have one, we recommend you a distribution like [Anaconda](https://www.continuum.io/why-anaconda).

### Fabric

You need also to have the tool [Fabric](http://www.fabfile.org) installed. If you have installed a Python distribution like [Anaconda](https://www.continuum.io/why-anaconda), you can install Fabric by typing:

    conda install fabric

### Python dependencies

Once done, download the RNA Science Toolbox and go into its directory. Now you need to install its Python dependencies. If you have installed the [Anaconda](https://www.continuum.io/why-anaconda) distribution, type:


    fab

This will launch the package manager conda to install dependencies. Otherwise, to use a package manager like pip, type:

    fab install:manager=pip

### Docker

Now you need to install the algorithms. This will be done with Docker images providing these algorithms fully configured. You first need to install the tool Docker. You will find all the details [here](https://docs.docker.com/engine/installation/).

### RNA algorithms

Now you need to install the Docker images containing the algorithms you want to work with. Here is the list of the images available:

 * [fjossinet/assemble2](https://hub.docker.com/r/fjossinet/assemble2/): provides algorithms like [RNAVIEW](http://ndbserver.rutgers.edu/ndbmodule/services/download/rnaview.html), [Vienna RNA package](https://www.tbi.univie.ac.at/RNA/), [foldalign](http://rth.dk/resources/foldalign/), [LocARNA](http://rna.informatik.uni-freiburg.de/LocARNA/)
 * [fjossinet/rnaseq](https://hub.docker.com/r/fjossinet/rnaseq/): provides algorithms like [SAMtools](http://samtools.sourceforge.net), [Tophat2](https://ccb.jhu.edu/software/tophat/), [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

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
</pre>
