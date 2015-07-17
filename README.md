RNA Science Toolbox -- A virtual environment to do RNA science
===============================================================

The RNA Science Toolbox provides a Python API (PyRNA) and a fully configured virtual machine (described in the Vagrantfile) to do RNA science. This virtual environment allows you to:
* write your own Python scripts based on PyRNA,
* deploy the PyRNA API as Web services to be consumed from graphical tool like [Assemble2](http://www.bioinformatics.org/assemble/).

The PyRNA API allows you to:

* parse RNA data from "classical" file formats (PDB, CT, FASTA, VIENNA,...) and convert them into easy-to-use and easy-to-analyse data structures:
    * [Pandas Series and DataFrames](http://pandas.pydata.org/pandas-docs/stable/dsintro.html)
    * "PyRNA objects": defined in the module pyrna.features.
* compute RNA data from RNA algorithms (see list below) and convert them into Pandas data structures and PyRNA objects,
* recover RNA data from public databases ([PDB](http://www.rcsb.org/pdb/home/home.do), [RFAM](http://rfam.sanger.ac.uk),...) and convert them into Pandas data structures and PyRNA objects,
* deploy some functionalities as REST Web services.

In its current state, PyRNA is able to handle several algorithms like:

* [The RNA Vienna Package](http://goo.gl/7frDgF)
* [Blastr](http://goo.gl/lKCR1u)
* [CONTRAfold](http://goo.gl/4BCI7)
* [Infernal](http://goo.gl/SxLHJO)
* [Mlocarna](http://goo.gl/AIGKrl)
* [RNAMotif](http://goo.gl/MDdOQ2)
* [RNAVIEW](http://goo.gl/c5o19v)
* ....

You can [follow this project on twitter](https://twitter.com/RnaSciToolbox) to get updates as they happen.

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

Congratulations, you're ready to .

Since the PyRNA directory is a shared folder between your computer and the fully configured virtual machine, we recommended you to:

* keep your Python scripts in the directory files/scripts/python provided with PyRNA,
* run your Python scripts from the shell of the virtual machine.

#The PyRNA Server

PyRNA provides you the ability to deploy its functionalities as REST Web services over a local network. From the command-line of the virtual machine, type:

    $ server.py

By default, the server runs on http://localhost:8080. Point a Web browser to this address to have a look at the website functionalities. The server provides graphical Web pages and Web services. The Web services can be consumed from dedicated tools like [Assemble2](http://www.bioinformatics.org/assemble/)
