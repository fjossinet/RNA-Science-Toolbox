import os, sys, re
from os.path import expanduser
from fabric.context_managers import cd
from fabric.contrib.console import confirm
from fabric.operations import prompt
from fabric.operations import run
from fabric.colors import green
from fabric.api import *

home = os.path.expanduser('~')
available_algorithms = ['rnaview', 'vienna', 'infernal', 'contrafold', \
                        'locarna', 'foldalign', 'trnascan-SE', 'blast', 'blastR', 'clustalw', \
                        'rnamotif', 'samtools', 'bowtie2', 'gotohscan', 'rnabob', 'bcheck', \
                        'snoreport', 'snoscan', 'snogps']

@task(default=True)
def basic_install():
    """
    Do a basic installation of the system. Not all the RNA algorithms will be installed.
    """
    update()
    python()
    algorithms(install=['rnaview', 'vienna', 'contrafold', 'foldalign', 'locarna'])
    website()
    PDB(limit=5000)
    RNA3DHub(limit=5000)

def full_install():
    """
    Do a full installation of the system. All the RNA algorithms will be installed.
    """
    if not confirm("This will install and configure all the dependencies. Do you wish to continue?") :
        print "Bye!"
        sys.exit()
    update()
    python()
    website()
    PDB()
    RNA3DHub()

def minimal_install():
    """
    Do a minimal installation of the system. No RNA algorithms will be installed, only the Python packages.
    """
    if not confirm("This will install and configure all the Python packages. Do you wish to continue?") :
        print "Bye!"
        sys.exit()
    python()

@task
def python(manager="conda"):
    """
    Install all the Python packages
    """
    print(green("Installing Python packages..."))
    if manager == "conda":
        local('conda config --set always_yes TRUE')
        local('conda install pandas')
        local('conda install pymongo')
        local('conda install pysam')
        local('conda install ujson')
    elif manager == "pip":
        local('pip install pandas')
        local('pip install pymongo')
        local('pip install pysam')
        local('pip install ujson')
    else:
        print "You need to install the following Python packages:"
        print ("pandas")
        print ("pymongo")
        print ("pysam")
        print ("ujson")

@task
def update():
    """
    Update the Ubuntu package list
    """
    print(green("Updating the package list..."))
    local('sudo apt-get update -qq')

@task
def algorithms(install=available_algorithms, list="False"):
    """
    Install all the RNA algorithms
    """
    if list == "True":
        print '\n'.join(available_algorithms)
        sys.exit()
    if isinstance(install, basestring):
        install = install.split(':')
    print(green("Installing the RNA algorithms..."))
    installation_directory = prompt('Where do you want to install your RNA algorithms?', default=os.path.join(home, 'algorithms'), validate=r'^'+expanduser('~')+'/.+/?$')
    if os.path.exists(installation_directory) and confirm("The directory %s already exist. Do you wish to continue?"%installation_directory) or not os.path.exists(installation_directory):
        if not os.path.exists(installation_directory):
            local('mkdir '+installation_directory)
        if not os.path.exists('%s/RNAVIEW/bin/rnaview'%installation_directory) and 'rnaview' in install:
            rnaview(installation_directory)
        if not os.path.exists('%s/ViennaRNA_2.1.8/bin/RNAfold'%installation_directory) and 'vienna' in install:
            vienna_rna_package(installation_directory)
        if not os.path.exists('%s/infernal_1.0.2/bin/cmsearch'%installation_directory) and 'infernal' in install :
            infernal(installation_directory)
        if not os.path.exists('%s/contrafold/src/contrafold'%installation_directory) and 'contrafold' in install:
            contrafold(installation_directory)
        if not os.path.exists('%s/locarna_1.8.1/bin/locarna'%installation_directory) and 'locarna' in install:
            locarna("%s/ViennaRNA"%installation_directory, installation_directory)
        if not os.path.exists('%s/foldalign.2.1.1/bin/foldalign'%installation_directory) and 'foldalign' in install:
            foldalign(installation_directory)
        if not os.path.exists('%s/tRNAscan-SE-1.3.1/bin/tRNAscan-SE'%installation_directory) and 'trnascan-SE' in install:
            trnaScanSE(installation_directory)
        if (not os.path.exists('%s/blast-2.2.26/bin/blastall'%installation_directory) or not os.path.exists('%s/ncbi-blast-2.2.31+/bin/blastp'%installation_directory)) and 'blast' in install:
            blast(installation_directory)
        if not os.path.exists('%s/blastR_package_V2.2/scripts/blastR.pl'%installation_directory) and 'blastR' in install:
            blastR(installation_directory)
        if not os.path.exists('%s/clustalw_2.1/bin/clustalw2'%installation_directory) and 'clustalw' in install:
            clustalw(installation_directory)
        if not os.path.exists('%s/rnamotif-3.0.7/src/rnamotif'%installation_directory) and 'rnamotif' in install:
            rnamotif(installation_directory)
        if not os.path.exists('%s/samtools-0.1.18/samtools'%installation_directory) and 'samtools' in install:
            samtools(installation_directory)
        if not os.path.exists('%s/bowtie2-2.2.2/bowtie2'%installation_directory) and 'bowtie2' in install:
            bowtie2(installation_directory)
        if not os.path.exists('%s/GotohScan_2.0-alpha/src/GotohScan2a'%installation_directory) and 'gotohscan' in install:
            gotohscan(installation_directory)
        if not os.path.exists('%s/rnabob-2.2.1/rnabob'%installation_directory) and 'rnabob' in install:
            rnabob(installation_directory)
        if not os.path.exists('%s/Bcheck-0.6/Bcheck'%installation_directory) and 'bcheck' in install:
            bcheck(installation_directory)
        if not os.path.exists('%s/SnoReport1.0/snoReport'%installation_directory) and 'snoreport' in install:
            snoreport(installation_directory)
        if not os.path.exists('%s/snoscan-0.9b/snoscan'%installation_directory) and 'snoscan' in install:
            snoscan(installation_directory)
        if not os.path.exists('%s/snoGPS-0.2/src/pseudoU_test'%installation_directory) and 'snogps' in install:
            snogps(installation_directory)
    print(green("The PATH variable has been updated in your .bashrc file"))

def bcheck(installation_directory = "%s/algorithms"%home):
    """
    Install BCheck
    """
    print(green("Installing BCheck..."))
    if not os.path.exists('%s/Bcheck/'%installation_directory):
        local('echo "export PATH=\$PATH:%s/Bcheck/" >> $HOME/.bashrc'%installation_directory)
    with lcd(installation_directory):
        #the Bcheck recovered from Dropbox has been modified to not check the rnabob and Infernal version
        local('wget -qO Bcheck-0.6.tar.gz http://dl.dropbox.com/u/3753967/algorithms/Bcheck-0.6.tar.gz')
        local('tar -xzf Bcheck-0.6.tar.gz')
        local('rm Bcheck-0.6.tar.gz')
        local('ln -sf %s/Bcheck-0.6 ./Bcheck'%installation_directory)
        local('echo "export PATH=\$PATH:%s/Bcheck/" >> $HOME/.bashrc'%installation_directory)

def blast(installation_directory = "%s/algorithms"%home):
    """
    Install Blast and NCBI-blast
    """
    print(green("Installing Blast..."))
    if not os.path.exists('%s/blast/bin'%installation_directory):
        local('echo "export PATH=\$PATH:%s/blast/bin" >> $HOME/.bashrc'%installation_directory)
    if not os.path.exists('%s/ncbi-blast/bin'%installation_directory):
        local('echo "export PATH=\$PATH:%s/ncbi-blast/bin" >> $HOME/.bashrc'%installation_directory)
    with lcd(installation_directory):
        local('wget -qO blast-2.2.26-x64-linux.tar.gz http://dl.dropbox.com/u/3753967/algorithms/blast-2.2.26-x64-linux.tar.gz')
        local('tar -xzf blast-2.2.26-x64-linux.tar.gz')
        local('rm blast-2.2.26-x64-linux.tar.gz')
        local('ln -sf %s/blast-2.2.26 ./blast'%installation_directory)
        local('wget -qO ncbi-blast-2.2.31+-x64-linux.tar.gz https://dl.dropboxusercontent.com/u/3753967/algorithms/ncbi-blast-2.2.31+-x64-linux.tar.gz')
        local('tar -xzf ncbi-blast-2.2.31+-x64-linux.tar.gz')
        local('rm ncbi-blast-2.2.31+-x64-linux.tar.gz')
        local('ln -sf %s/ncbi-blast-2.2.31+ ./ncbi-blast'%installation_directory)

def blastR(installation_directory = "%s/algorithms"%home):
    """
    Install BlastR
    """
    print(green("Installing BlastR..."))
    if not os.path.exists('%s/blastR/scripts/'%installation_directory):
        local('echo "export PATH=\$PATH:%s/blastR/scripts/" >> $HOME/.bashrc'%installation_directory)
    with lcd(installation_directory):
        local('wget -qO blastR_package_V2.2.tar.gz http://dl.dropboxusercontent.com/u/3753967/algorithms/blastR_package_V2.2.tar.gz')
        local('tar -xzf blastR_package_V2.2.tar.gz')
        local('rm blastR_package_V2.2.tar.gz')
        local('ln -sf %s/blastR_package_V2.2 ./blastR'%installation_directory)

def bowtie2(installation_directory = "%s/algorithms"%home):
    """
    Install Bowtie2
    """
    print(green("Installing Bowtie2..."))
    if not os.path.exists('%s/bowtie2/'%installation_directory):
        local('echo "export PATH=\$PATH:%s/bowtie2/" >> $HOME/.bashrc'%installation_directory)
    with lcd(installation_directory):
        local('wget -qO bowtie2-2.2.2.tar.gz https://dl.dropboxusercontent.com/u/3753967/algorithms/bowtie2-2.2.2.tar.gz')
        local('tar -xzf bowtie2-2.2.2.tar.gz')
        local('rm bowtie2-2.2.2.tar.gz')
        local('ln -sf %s/bowtie2-2.2.2 ./bowtie2'%installation_directory)
        with lcd("bowtie2-2.2.2"):
            local('make')

def clustalw(installation_directory = "%s/algorithms"%home):
    """
    Install Clustalw
    """
    print(green("Installing ClustalW..."))
    if not os.path.exists('%s/clustalw/bin'%installation_directory):
        local('echo "export PATH=\$PATH:%s/clustalw/bin" >> $HOME/.bashrc'%installation_directory)
    with lcd(installation_directory):
        local('wget -qO clustalw-2.1.tar.gz http://dl.dropbox.com/u/3753967/algorithms/clustalw-2.1.tar.gz')
        local('tar -xzf clustalw-2.1.tar.gz')
        local('rm clustalw-2.1.tar.gz')
        local('ln -sf %s/clustalw_2.1 ./clustalw'%installation_directory)
        with lcd('clustalw-2.1'):
            local('./configure --prefix=%s/clustalw_2.1'%installation_directory)
            local('make')
            local('make install')
        with lcd(installation_directory):
            local('rm -rf clustalw-2.1')

#problem to compile contrafold. For now, the dropbox provides a precompiled version for linux.
def contrafold(installation_directory = "%s/algorithms"%home):
    """
    Install Contrafold
    """
    print(green("Installing Contrafold..."))
    if not os.path.exists('%s/contrafold/src'%installation_directory):
        local('echo "export PATH=\$PATH:%s/contrafold/src" >> $HOME/.bashrc'%installation_directory)
    with lcd(installation_directory):
        local('wget -qO contrafold.tar.gz http://dl.dropbox.com/u/3753967/algorithms/contrafold.tar.gz')
        local('tar -xzf contrafold.tar.gz')
        local('rm contrafold.tar.gz')

def foldalign(installation_directory = "%s/algorithms"%home):
    """
    Install Foldalign
    """
    print(green("Installing Foldalign..."))
    if not os.path.exists('%s/foldalign/bin'%installation_directory):
        local('echo "export PATH=\$PATH:%s/foldalign/bin" >> $HOME/.bashrc'%installation_directory)
    with lcd(installation_directory):
        local('wget -qO foldalign.2.1.1.tar.gz http://dl.dropbox.com/u/3753967/algorithms/foldalign.2.1.1.tar.gz')
        local('tar -xzf foldalign.2.1.1.tar.gz')
        local('rm foldalign.2.1.1.tar.gz')
        local('ln -sf %s/foldalign.2.1.1 ./foldalign'%installation_directory)
        with lcd('foldalign.2.1.1'):
            local('make')

def gotohscan(installation_directory = "%s/algorithms"%home):
    """
    Install GotohScan
    """
    print(green("Installing GotohScan..."))
    if not os.path.exists('%s/GotohScan/src'%installation_directory):
        local('echo "export PATH=\$PATH:%s/GotohScan/src" >> $HOME/.bashrc'%installation_directory)
    with lcd(installation_directory):
        local('wget -qO GotohScan_2.0-alpha.tar.gz http://dl.dropbox.com/u/3753967/algorithms/GotohScan_2.0-alpha.tar.gz')
        local('tar -xzf GotohScan_2.0-alpha.tar.gz')
        local('rm GotohScan_2.0-alpha.tar.gz')
        local('ln -sf %s/GotohScan_2.0-alpha ./GotohScan'%installation_directory)
        with lcd('GotohScan_2.0-alpha'):
            local('make')

def infernal(installation_directory = "%s/algorithms"%home):
    """
    Install Infernal
    """
    print(green("Installing Infernal..."))
    if not os.path.exists('%s/infernal/bin'%installation_directory):
        local('echo "export PATH=\$PATH:%s/infernal/bin" >> $HOME/.bashrc'%installation_directory)
    with lcd(installation_directory):
        local('wget -qO infernal-1.0.2.tar.gz http://dl.dropbox.com/u/3753967/algorithms/infernal-1.0.2.tar.gz')
        local('tar -xzf infernal-1.0.2.tar.gz')
        local('rm infernal-1.0.2.tar.gz')
        local('ln -sf %s/infernal_1.0.2 ./infernal'%installation_directory)
        with lcd('infernal-1.0.2'):
            local('./configure --prefix="%s"/infernal_1.0.2'%installation_directory)
            local('make clean')
            local('make')
            local('make install')
        with lcd(installation_directory):
            local('rm -rf infernal-1.0.2')

def locarna(vrna_path, installation_directory = "%s/algorithms"%home):
    """
    Install Locarna
    """
    print(green("Installing Locarna..."))
    if not os.path.exists('%s/locarna/bin'%installation_directory):
        local('echo "export PATH=\$PATH:%s/locarna/bin" >> $HOME/.bashrc'%installation_directory)
    with lcd(installation_directory):
        local('wget -qO locarna-1.8.1.tar http://dl.dropbox.com/u/3753967/algorithms/locarna-1.8.1.tar')
        local('tar -xf locarna-1.8.1.tar')
        local('rm locarna-1.8.1.tar')
        local('ln -sf %s/locarna_1.8.1 ./locarna'%installation_directory)
        with lcd('locarna-1.8.1'):
            local('./configure --prefix=%s/locarna_1.8.1 --with-vrna=%s --without-perl --without-forester --without-kinfold'%(installation_directory, vrna_path))
            local('make clean')
            local('make')
            local('make install')
        with lcd(installation_directory):
            local('rm -rf locarna-1.8.1')

def rnamotif(installation_directory = "%s/algorithms"%home):
    """
    Install RNAMotif
    """
    print(green("Installing RNAMotif..."))
    if not os.path.exists('%s/rnamotif/src'%installation_directory):
        local('echo "export PATH=\$PATH:%s/rnamotif/src" >> $HOME/.bashrc'%installation_directory)
    local("sudo apt-get -y install flex")
    with lcd(installation_directory):
        local('wget -qO rnamotif-3.0.7.tar.gz https://dl.dropboxusercontent.com/u/3753967/algorithms/rnamotif-3.0.7.tar.gz')
        local('tar -xzf rnamotif-3.0.7.tar.gz')
        local('rm rnamotif-3.0.7.tar.gz')
        local('ln -sf %s/rnamotif-3.0.7 ./rnamotif'%installation_directory)
        with lcd('rnamotif-3.0.7'):
            local('make')

def rnabob(installation_directory = "%s/algorithms"%home):
    """
    Install RNABOB
    """
    print(green("Installing RNABOB..."))
    if not os.path.exists('%s/rnabob/'%installation_directory):
        local('echo "export PATH=\$PATH:%s/rnabob/" >> $HOME/.bashrc'%installation_directory)
    with lcd(installation_directory):
        local('wget -qO rnabob-2.2.1.tar.gz http://dl.dropbox.com/u/3753967/algorithms/rnabob.tar.gz')
        local('tar -xzf rnabob-2.2.1.tar.gz')
        local('rm rnabob-2.2.1.tar.gz')
        local('ln -sf %s/rnabob-2.2.1 ./rnabob'%installation_directory)
        with lcd('rnabob-2.2.1'):
            local('make')

def rnaview(installation_directory = "%s/algorithms"%home):
    """
    Install RNAView
    """
    print(green("Installing RNAView..."))
    if not os.path.exists('%s/RNAVIEW/bin'%installation_directory):
        local('echo "export RNAVIEW=%s/RNAVIEW/" >> $HOME/.bashrc'%installation_directory)
        local('echo "export PATH=\$PATH:\$RNAVIEW/bin" >> $HOME/.bashrc')
    with lcd(installation_directory):
        local('wget -qO RNAVIEW.tar.gz http://dl.dropbox.com/u/3753967/algorithms/RNAVIEW.tar.gz')
        local('tar -xzf RNAVIEW.tar.gz')
        local('rm RNAVIEW.tar.gz')
        with lcd('RNAVIEW'):
            local('make clean')
            local('make')

def samtools(installation_directory = "%s/algorithms"%home):
    """
    Install Samtools
    """
    print(green("Installing Samtool..."))
    if not os.path.exists('%s/samtools/'%installation_directory):
        local('echo "export PATH=\$PATH:%s/samtools/" >> $HOME/.bashrc'%installation_directory)
    with lcd(installation_directory):
        local('wget -qO samtools-0.1.18.tar http://dl.dropbox.com/u/3753967/algorithms/samtools-0.1.18.tar')
        local('tar -xf samtools-0.1.18.tar')
        local('rm samtools-0.1.18.tar')
        local('ln -sf %s/samtools-0.1.18 ./samtools'%installation_directory)
        with lcd('samtools-0.1.18'):
            local('make')

def snogps(installation_directory = "%s/algorithms"%home):
    """
    Install SnoGPS
    """
    print(green("Installing SnoGPS..."))
    if not os.path.exists('%s/snoGPS/'%installation_directory):
        local('echo "export PATH=\$PATH:%s/snoGPS/" >> $HOME/.bashrc'%installation_directory)
    with lcd(installation_directory):
        local('wget -qO snoGPS-0.2.tar.gz https://dl.dropboxusercontent.com/u/3753967/algorithms/snoGPS-0.2.tar.gz')
        local('tar -xzf snoGPS-0.2.tar.gz')
        local('rm snoGPS-0.2.tar.gz')
        local('ln -sf %s/snoGPS-0.2 ./snoGPS'%installation_directory)
        with lcd('snoGPS-0.2/src'):
            local('make')

def snoreport(installation_directory = "%s/algorithms"%home):
    """
    Install snoReport
    """
    print(green("Installing snoReport..."))
    if not os.path.exists('%s/SnoReport/'%installation_directory):
        local('echo "export PATH=\$PATH:%s/SnoReport/" >> $HOME/.bashrc'%installation_directory)
    with lcd(installation_directory):
        local('wget -qO SnoReport1.0.tgz http://dl.dropboxusercontent.com/u/3753967/algorithms/SnoReport1.0.tgz')
        local('tar -xzf SnoReport1.0.tgz')
        local('rm SnoReport1.0.tgz')
        local('ln -sf %s/SnoReport1.0 ./SnoReport'%installation_directory)
        with lcd('SnoReport1.0'):
            #SnoReport needs the Vienna package in the PATH. How to do it?? TO DO!!!!!
            local('./install.sh')

def snoscan(installation_directory = "%s/algorithms"%home):
    """
    Install Snoscan
    """
    print(green("Installing Snoscan..."))
    if not os.path.exists('%s/snoscan/'%installation_directory):
        local('echo "export PATH=\$PATH:%s/snoscan/" >> $HOME/.bashrc'%installation_directory)
    with lcd(installation_directory):
        local('wget -qO snoscan-0.9b.tar.gz https://dl.dropboxusercontent.com/u/3753967/algorithms/snoscan.tar.gz')
        local('tar -xzf snoscan-0.9b.tar.gz')
        local('rm snoscan-0.9b.tar.gz')
        local('ln -sf %s/snoscan-0.9b ./snoscan'%installation_directory)
        with lcd('snoscan-0.9b/squid-1.5j'):
            local('make')
        with lcd('snoscan-0.9b/'):
            local('make')

def trnaScanSE(installation_directory = "%s/algorithms"%home):
    """
    Install tRNAscan-SE
    """
    print(green("Installing tRNAscan-SE..."))
    if not os.path.exists('%s/tRNAscan-SE/'%installation_directory):
        local('echo "export PATH=\$PATH:%s/tRNAscan-SE/bin" >> $HOME/.bashrc'%installation_directory)
        local('echo "export PERL5LIB=%s/tRNAscan-SE/bin/:\$PERL5LIB" >> $HOME/.bashrc'%installation_directory)
    with lcd(installation_directory):
        local('wget -qO tRNAscan-SE-1.3.1.tar https://dl.dropboxusercontent.com/u/3753967/algorithms/tRNAscan-SE-1.3.1.tar')
        local('tar -xf tRNAscan-SE-1.3.1.tar')
        local('rm tRNAscan-SE-1.3.1.tar')
        local('ln -sf %s/tRNAscan-SE-1.3.1 ./tRNAscan-SE'%installation_directory)
        with lcd('tRNAscan-SE-1.3.1'):
            #edition of the Makefile
            local('sed "s,\$(HOME),%s/tRNAscan-SE-1.3.1," Makefile > Makefile2'%installation_directory)
            local('mv Makefile2 Makefile')
            local('make')
            local('make install')

def vienna_rna_package(installation_directory = "%s/algorithms"%expanduser("~")):
    """
    Install the Vienna RNA package
    """
    print(green("Installing Vienna RNA package..."))
    if not os.path.exists('%s/ViennaRNA/bin'%installation_directory):
        local('echo "export PATH=\$PATH:%s/ViennaRNA/bin" >> $HOME/.bashrc'%installation_directory)
        local('echo "export PATH=\$PATH:%s/ViennaRNA/share/ViennaRNA/bin/" >> $HOME/.bashrc'%installation_directory)
    with lcd(installation_directory):
        local('wget -qO ViennaRNA-2.1.8.tar.gz http://dl.dropbox.com/u/3753967/algorithms/ViennaRNA-2.1.8.tar.gz')
        local('tar -xzf ViennaRNA-2.1.8.tar.gz')
        local('rm ViennaRNA-2.1.8.tar.gz')
        local('ln -sf %s/ViennaRNA_2.1.8 ./ViennaRNA'%installation_directory)
        with lcd('ViennaRNA-2.1.8'):
            local('./configure --prefix="%s"/ViennaRNA_2.1.8'%installation_directory)
            local('make clean')
            local('make')
            local('make install')
        with lcd(installation_directory):
            local('rm -rf ViennaRNA-2.1.8')

@task
def PDB(limit = 5000):
    """
    Feed the database with PDB data
    """
    print(green("Feed the database with %i 3D junctions..."%limit))
    if not confirm("The current PDB database will be erased. Do you wish to continue?") :
        print "Bye!"
        sys.exit()
    local('mongo PDB --eval "db.dropDatabase()"')
    local("import_3Ds.py -annotate -l %i"%limit)

@task
def RNA3DHub(limit = 5000):
    """
    Feed the database with PDB data derived from the RNA3DHub website
    """
    print(green("Feed the database with %i 3D junctions derived from the RNA3DHub website..."%limit))
    if not confirm("The current RNA3DHub database will be erased. Do you wish to continue?") :
        print "Bye!"
        sys.exit()
    local('mongo RNA3DHub --eval "db.dropDatabase()"')
    local("import_3Ds.py -annotate -rna3dhub -l %i"%limit)

@task
def website():
    """
    Install the website
    """
    print(green("Installing a full web stack..."))

    print(green("Installing the Python packages..."))
    local('conda config --set always_yes TRUE')
    local('conda install ipython')
    local('conda install ipython-notebook')
    local('conda install tornado')

    print(green("Installing Node.js..."))
    local("sudo apt-get -y install nodejs npm")
    local("sudo ln -s `which nodejs` /usr/bin/node")

    print(green("Installing Bower..."))
    local("sudo npm install -g bower")

    print(green("Installing the website..."))
    local('cd /vagrant/website ; bower --config.interactive=false install')

@task
def backup():
    """
    Dump the database and store the data in the a folder /vagrant/backup accessible from the host.
    """
    print(green("The database is dumped in the folder /vagrant/backup..."))
    if not os.path.exists('/vagrant/backup'):
        local('mkdir /vagrant/backup')
    local('mongodump --out /vagrant/backup/')

@task
def restore():
    """
    Restore the database from a backup stored in the folder /vagrant/backup
    """
    print(green("The database is restored..."))
    if os.path.exists('/vagrant/backup'):
        local('mongorestore /vagrant/backup/')
