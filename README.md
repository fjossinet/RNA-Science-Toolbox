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

PyRNA has been designed to be used on a UNIX system (Linux, MacOSX,...).

In its current state, PyRNA is able to handle several algorithms:

* [Blast](ftp://ftp.ncbi.nlm.nih.gov/blast/)
* [Blastr](http://goo.gl/lKCR1u)
* [Bowtie](http://goo.gl/nmXKH)
* [Clustalw](http://goo.gl/Z9FRV)
* [CONTRAfold](http://goo.gl/4BCI7)
* [Gotohscan](http://goo.gl/2atKpi) 
* [Infernal](http://goo.gl/SxLHJO)
* [Mlocarna](http://goo.gl/AIGKrl)
* [RNA Vienna Package (RNAfold, RNAplot)](http://goo.gl/7frDgF)
* [RNAMotif](http://goo.gl/MDdOQ2)
* [RNAVIEW](http://goo.gl/c5o19v)
* [snoGPS](http://goo.gl/66pnrF)
* [SnoReport](http://goo.gl/pq3qXu)
* [Snoscan](http://goo.gl/P5EQiH)
* more to come....

You can [follow me on twitter](https://twitter.com/fjossinet) to get updates as they happen.

#Quick Start

* install [Vagrant](https://www.vagrantup.com/) and [Virtualbox](https://www.virtualbox.org/) on your computer. These tools will allow you to launch a fully configured Linux as a virtual machine. 

* in a UNIX shell, recover the last version of PyRNA by typing:

        $ git clone https://github.com/fjossinet/PyRNA.git

* in the PyRNA directory, type:

        $ vagrant up
        $ vagrant ssh        

This will launch and log you into a fully configured virtual machine.

* from the shell of the virtual machine, type:

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

Congratulations, you're ready to develop your own Python scripts with PyRNA. 

Since the PyRNA directory is a shared folder between your computer and the fully configured virtual machine, we recommended you to:

* keep your Python scripts in the directory files/scripts/python provided with PyRNA,
* run your Python scripts from the shell of the virtual machine.

#The PyRNA modules

* features.py: classes to model the RNA concepts used by PyRNA (a.k.a. "PyRNA objects"),
* parsers.py: functions to parse and convert RNA data,
* computations.py: classes to wrap RNA algorithms,
* db.py: classes to connect public databases,
* server.py: executable script to launch the embedded REST server

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

    $ ./server.py [-wh webserver_host (default: localhost)] [-wp webserver_port (default: 8080)] [-mh mongodb_host (default: localhost)] [-mp mongodb_port (default: 27017)] [-conf configuration_file] 

The configuration file allows you to define which algorithms are enabled for remote computations. A sample configuration is provided with PyRNA (named pyrna.conf).

Examples: 

    $ ./server.py

    $ ./server.py -wh my_host_name -wp 80

    $ ./server.py -wh my_host_name -wp 80 -conf $PYRNA_HOME/pyrna.conf
