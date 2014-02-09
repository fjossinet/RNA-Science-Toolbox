PyRNA -- Produce and Analyse RNA Data with Python
=================================================

PyRNA allows you to:

* parse RNA data from "classical" file formats (PDB, CT, FASTA, VIENNA,...) and convert them into easy-to-use and easy-to-analyse data structures (like [Pandas Series and DataFrames](http://pandas.pydata.org/pandas-docs/stable/dsintro.html) and PyRNA objects),
* compute RNA data from RNA algorithms ([Vienna RNA package](http://www.tbi.univie.ac.at/RNA/), [RNAVIEW](http://ndbserver.rutgers.edu/services/download/), [Contrafold](http://contra.stanford.edu/contrafold/),...),
* extract RNA data from public databases ([PDB](http://www.rcsb.org/pdb/home/home.do), [RFAM](http://rfam.sanger.ac.uk),...)
* export RNA data to the [JSON](http://www.json.org) format,
* submit and manage grid jobs using the [glite middleware](http://glite.cern.ch), 
* deploy some PyRNA functionalities as REST Web services.

To consume the REST Web Services, you can have a look on the code source of [very simple  REST clients](https://bitbucket.org/fjossinet/pyrna-rest-clients) or you can use a Java GUI like [Assemble2](http://www.bioinformatics.org/assemble/).

PyRNA has been designed to be used on a UNIX system (Linux, MacOSX,...). It has been developed with MacOSX and has been tested on MacOSX 10.8 and Scientific Linux 5.8 with python 2.7 and gcc 4.1 and 4.2.

#Dependencies

To be able to use 99% of all PyRNA functionalities, you will only need one Python dependency:

* [Pandas](http://pandas.pydata.org/): the data structures from pandas are used to describe RNA data,

If you want to handle NGS data, you will also need to install [pysam](https://code.google.com/p/pysam/), a Python interface for the SAM/BAM sequence alignment and mapping format.

A first option is to install and use a "scientific" Python distribution like [Canopy](https://www.enthought.com/products/canopy/) or [Anaconda](https://store.continuum.io/cshop/anaconda/). They provide access to numerous scientific packages that are already installed. 

A second option is to install the dependencies yourself for your Python distribution. We suggest you to install first the command [easy_install](http://pythonhosted.org/distribute/easy_install.html) for your Python distribution. Then everything will be easier to install:

* sudo easy_install pandas
* sudo easy_install pysam

You will also need a bunch of RNA algorithms installed and available in your PATH. You have two options:

* to install them yourself,
* to run the script "install_algorithms.sh" located in files/scripts/shell. To be able to run this script efficiently, you will need to have compilation tools installed (gcc and make).

PyRNA needs to find the following algorithms in your PATH:

* [Blast](ftp://ftp.ncbi.nlm.nih.gov/blast/)
* [Blastr](http://goo.gl/lKCR1u)
* [Bowtie](http://goo.gl/nmXKH)
* [Clustalw](http://goo.gl/Z9FRV)
* [CONTRAfold](http://goo.gl/4BCI7)
* [Gotohscan](http://goo.gl/2atKpi) 
* [Infernal](http://goo.gl/SxLHJO)
* [Mlocarna](http://goo.gl/AIGKrl)
* [RNA Vienna Package](http://goo.gl/7frDgF)
* [RNAMotif](http://goo.gl/MDdOQ2)
* [RNAVIEW](http://goo.gl/c5o19v)
* [snoGPS](http://goo.gl/66pnrF)
* [SnoReport](http://goo.gl/pq3qXu)
* [Snoscan](http://goo.gl/P5EQiH)

You don't need to install all of them. It will depend on the classes you will import from the module pyrna.computations in your scripts. 

#Quick Start

* if you don't have it, install the command [easy_install](http://pythonhosted.org/distribute/easy_install.html)
* install the Python dependencies Pandas and pysam:
    * sudo easy_install pandas
    * sudo easy_install pysam
* in a UNIX shell, type:

        $ hg clone https://bitbucket.org/fjossinet/pyrna

* in the configuration file of your shell (.bashrc for example), add the following lines:

    * export PYRNA_HOME=the_full_path_of_your_pyrna_location
    * export PYTHONPATH=$PYTHONPATH:$PYRNA_HOME
    * export PATH=$PATH:$PYRNA_HOME/pyrna:$PYRNA_HOME/files/scripts/python:$PYRNA_HOME/files/scripts/shell

* reload your shell,
* type:

        $ install_algorithms.sh your_RNA_algorithms_location

* in the configuration file of your shell (.bashrc for example), add the following lines:

    * source your_RNA_algorithms_location/setmyenv

* reload your shell,
* type:

        $ pyrna_tests.py

You should get an output like:

        RNA sequence from 1EHZ:

        GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA

        List of base-pairs computed from RNAfold (RNA Vienna Package):

              orientation edge1 edge2  pos1  pos2
        0            c     (     )    20    27
        1            c     (     )    19    28
        2            c     (     )    17    30
        3            c     (     )    16    31
        4            c     (     )    15    32
        [...]

#The PyRNA modules

* features.py: classes to model the RNA concepts used by PyRNA,
* parsers.py: functions to parse and convert RNA data,
* computations.py: classes to wrap RNA algorithms,
* db.py: classes to connect public databases,
* restserver.py: executable script to launch the embedded REST server, 
* task.py: a class allowing to design grid tasks,
* glite.py: executable file and functions to manage glite jobs.

#The PyRNA scripts

PyRNA provides several scripts in the files directory:

* scripts/python: examples of Python scripts using PyRNA
* scripts/shell: shell scripts  
* grid_tasks: these are Python scripts specializing the Task class (see the task.py module). These scripts parallelize some tasks and have been designed to be executed on a grid. 

#The PyRNA REST Server

PyRNA provides the ability to deploy some of its functionalities as REST webservices. You will need to install: 

* [MongoDB](http://www.mongodb.org/): a NoSQL database,
* [PyMongo](http://api.mongodb.org/python/current/): a Python dependency providing an easy way to connect MongoDB databases,
* [Flask](http://flask.pocoo.org/): a Python microframework,
* [ujson](https://pypi.python.org/pypi/ujson): a Python dependency to load/dump JSON data,
* [Scipy](http://www.scipy.org/): a Python-based ecosystem of open-source software for mathematics, science, and engineering

Once evrything installed, launch MongoDB, then type: 

    $ ./restserver.py [-wh webserver_host (default: localhost)] [-wp webserver_port (default: 8080)] [-mh mongodb_host (default: localhost)] [-mp mongodb_port (default: 27017)]

Examples: 

    $ ./restserver.py

    $ ./restserver.py -wh my_host_name -wp 80
