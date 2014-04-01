PyRNA -- Produce and Mine RNA Data with Python
=================================================

PyRNA allows you to:

* parse RNA data from "classical" file formats (PDB, CT, FASTA, VIENNA,...) and convert them into easy-to-use and easy-to-analyse data structures:
    * [Pandas Series and DataFrames](http://pandas.pydata.org/pandas-docs/stable/dsintro.html)
    * "PyRNA objects": defined in the module pyrna.features.
* compute RNA data from RNA algorithms (see list below) and convert them into Pandas data structures and PyRNA objects, 
* recover RNA data from public databases ([PDB](http://www.rcsb.org/pdb/home/home.do), [RFAM](http://rfam.sanger.ac.uk),...) and convert them into Pandas data structures and PyRNA objects,
* deploy some functionalities as REST Web services.

[Assemble2](http://www.bioinformatics.org/assemble/) is an example of Java graphical client consuming PyRNA Web services.

To learn more about PyRNA, check the [PyRNA Cookbook](http://goo.gl/q20VoF)

PyRNA has been designed to be used on a UNIX system (Linux, MacOSX,...). It has been developed with MacOSX and tested on MacOSX 10.8 and Scientific Linux 5.8 (python 2.7, gcc 4.1 or 4.2).

In its current state, PyRNA is able to handle the following algorithms:

* [Blast](ftp://ftp.ncbi.nlm.nih.gov/blast/)
* [Blastr](http://goo.gl/lKCR1u)
* [Bowtie](http://goo.gl/nmXKH)
* [Clustalw](http://goo.gl/Z9FRV)
* [CONTRAfold*](http://goo.gl/4BCI7)
* [Gotohscan](http://goo.gl/2atKpi) 
* [Infernal](http://goo.gl/SxLHJO)
* [Mlocarna](http://goo.gl/AIGKrl)
* [RNA Vienna Package (RNAfold*, RNAplot*)](http://goo.gl/7frDgF)
* [RNAMotif](http://goo.gl/MDdOQ2)
* [RNAVIEW*](http://goo.gl/c5o19v)
* [snoGPS](http://goo.gl/66pnrF)
* [SnoReport](http://goo.gl/pq3qXu)
* [Snoscan](http://goo.gl/P5EQiH)

The algorithms highlighted with a * can be used directly from my own server (arn-ibmc.in2p3.fr), without any installation. You just have to configure your Python scripts properly (see "Create secondary structures from algorithms" in [this IPython notebook](http://goo.gl/WHpfWh)).

You can [follow me on twitter](https://twitter.com/fjossinet) to get updates as they happen.

#Dependencies

To be able to use 99% of all PyRNA functionalities, you will only need a single Python dependency:

* [Pandas](http://pandas.pydata.org/): the data structures from pandas are used to describe RNA data,

If you want to handle NGS data, you will also need to install [pysam](https://code.google.com/p/pysam/), a Python interface for the SAM/BAM sequence alignment and mapping format.

A first option is to install and use a "scientific" Python distribution like [Canopy](https://www.enthought.com/products/canopy/) or [Anaconda](https://store.continuum.io/cshop/anaconda/). They provide access to numerous pre-installed scientific packages. 

A second option is to install the dependencies yourself for your Python distribution. I recommend you to first install the command [easy_install](http://pythonhosted.org/distribute/easy_install.html). Then you can type:

* sudo easy_install pandas
* sudo easy_install pysam

#Quick Start

* if you don't have it, install the command [easy_install](http://pythonhosted.org/distribute/easy_install.html)
* install the Python dependencies Pandas and pysam (if you plan to handle NGS data):
    * sudo easy_install pandas
    * sudo easy_install pysam
* in a UNIX shell, type:

        $ git clone https://github.com/fjossinet/PyRNA.git

* in the configuration file of your shell (.bashrc for example), add the following lines:

    * export PYRNA_HOME=the_full_path_of_your_pyrna_location
    * export PYTHONPATH=$PYTHONPATH:$PYRNA_HOME
    * export PATH=$PATH:$PYRNA_HOME/pyrna:$PYRNA_HOME/files/scripts/python:$PYRNA_HOME/files/scripts/shell:$PYRNA_HOME/files/scripts/grid_tasks:

* reload your shell,
* type:

        $ pyrna_tests.py

You should get an output like:

        Recovering entry 1EHZ from Protein Databank...

        ## 3D annotation ##

        List of base-pairs computed with RNAVIEW:

        edge1 edge2 orientation  pos1  pos2
        0      (     )           c     1    72
        1      (     )           c     2    71
        2      (     )           c     3    70
        3      (     )           c     4    69
        4      (     )           c     5    68
        5      (     )           c     6    67
        [...]

        ## 2D prediction ##

        RNA sequence from 1EHZ:

        GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA

        List of base-pairs computed with RNAfold (RNA Vienna Package):

        orientation edge1 edge2  pos1  pos2
        0            c     (     )    20    27
        1            c     (     )    19    28
        2            c     (     )    17    30
        3            c     (     )    16    31
        4            c     (     )    15    32
        5            c     (     )    14    33
        6            c     (     )    13    34
        [...]

#Local installation of the RNA algorithms

To increase the speed of your computations and/or to avoid to flood my own server, i recommend you to install the RNA algorithms locally, on your own computer. You can download and install the algorithms yourself or you can use the script "install_algorithms.sh" provided with PyRNA (and located in files/scripts/shell). This script should simplify the installation process, but this is without any guarantee. To be able to run this script efficiently, you will need to have compilation tools installed (gcc, g++ and make). You don't need to install all of the algorithms listed above. It will depend on the classes you will import from the module pyrna.computations in your scripts. Once PyRNA installed successfully (see the Quickstart section above):

* in a terminal, type:

        $ install_algorithms.sh your_RNA_algorithms_location

* once done, in the configuration file of your shell (.bashrc for example), add the following line:

    * source your_RNA_algorithms_location/setmyenv

* reload your shell.

#The PyRNA modules

* features.py: classes to model the RNA concepts used by PyRNA (a.k.a. "PyRNA objects"),
* parsers.py: functions to parse and convert RNA data,
* computations.py: classes to wrap RNA algorithms,
* db.py: classes to connect public databases,
* restserver.py: executable script to launch the embedded REST server

The modules task.py and glite.py allow to design, submit and manage grid jobs using the [gLite](http://glite.web.cern.ch/glite/) middleware. As a "basic user", you should not care about them.

#The PyRNA scripts

PyRNA provides several scripts in the files directory:

* scripts/python: examples of Python scripts using PyRNA
* scripts/shell: shell scripts  
* grid_tasks: Python scripts that parallelize some "large" tasks (genomic annotation for example). As a "basic user", you should not care about them.

#The PyRNA REST Server

PyRNA provides you the ability to deploy some functionalities as REST Web services. This allows you to install and configure PyRNA on a server of your local network and to allow other computers to connect it. You will need to install: 

* [MongoDB](http://www.mongodb.org/): a NoSQL database,
* [PyMongo](http://api.mongodb.org/python/current/): a Python dependency providing an easy way to connect MongoDB databases,
* [Tornado](http://www.tornadoweb.org/): a Python web framework and asynchronous networking library,
* [Flask](http://flask.pocoo.org/): a Python microframework,
* [ujson](https://pypi.python.org/pypi/ujson): a Python dependency to load/dump JSON data,

Once everything installed, launch your MongoDB and type: 

    $ ./restserver.py [-wh webserver_host (default: localhost)] [-wp webserver_port (default: 8080)] [-mh mongodb_host (default: localhost)] [-mp mongodb_port (default: 27017)] [-conf configuration_file] 

The configuration file allows you to define which algorithms are enabled for remote computations. A sample configuration is provided with PyRNA (named pyrna.conf).

Examples: 

    $ ./restserver.py

    $ ./restserver.py -wh my_host_name -wp 80

    $ ./restserver.py -wh my_host_name -wp 80 -conf $PYRNA_HOME/pyrna.conf
