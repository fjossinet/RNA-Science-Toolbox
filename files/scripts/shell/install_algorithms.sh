#!/bin/bash

if [ $# -ne 1 ]
then
  echo "Usage: `basename $0` directory_name_with_absolute_path"
  exit -1
fi

DEST=$1

if [[ $DEST != /* ]]
then
	echo "Your installation directory has to be given with its full path"
	exit -1
fi

if [ -e $DEST ]
then
	read -p "$DEST already exists. Do you want to continue? [Y/n]" -n 1 -r
	if [[ $REPLY =~ ^[n]$ ]]
	then
		echo -e "\nInstallation aborted!!"
	    exit -1
	fi
else
	mkdir $DEST
	if [ ! -e $DEST ]
	then
		echo "I cannot create $DEST"
		exit -1	
	fi
fi

mkdir -p "$DEST/tmp"

export PYRNA_ALGORITHMS=$DEST

echo "export PYRNA_ALGORITHMS=$DEST" >> "$DEST/setmyenv"

export R_LIBS=$PYRNA_ALGORITHMS/R_LIBS

if [ ! -e $R_LIBS ]
then
	mkdir $R_LIBS
fi

#we add the user lib for R to the list of libpath of R
echo ".libPaths(c(\"$R_LIBS\",.libPaths()))" > $DEST/tmp/configure.R

R -q --save < $DEST/tmp/configure.R

echo "export R_LIBS=$PYRNA_ALGORITHMS/R_LIBS" >> "$DEST/setmyenv"

if [ "$(which make)" == "" ] || [ "$(which gcc)" == "" ]
then
    echo "You don't have make and gcc to compile the RNA algorithms."
        case "`uname`" in
            Darwin*)
                echo "You have to install the XCode Tools available as additional packages in your MacOSX install disc."
                ;;
            *)
                ;;
    esac
    echo "For any questions or remarks, please contact f.jossinet@ibmc-cnrs.unistra.fr"
    exit -1
fi

RESULT="NAME\t\tSUCCEED?\n"
RESULT=$RESULT"----\t\t-------\n"

CURRENTDIR=$(pwd)

if [ "$(which wget)" == "" ]
then
    echo "You need wget to download the RNA algorithms."
        case "`uname`" in
            Darwin*)
                echo "You can get it from MacPorts (http://www.macports.org/)."
                ;;
            *)
                ;;
        esac
    echo "For any questions or remarks, please contact f.jossinet@ibmc-cnrs.unistra.fr"
    exit -1
fi

FETCHER="wget -qO "

cd "$DEST"

################# INFERNAL ##########################

cd "$DEST"

if [ -f infernal/bin/cmsearch ]
then
    RESULT=$RESULT"Infernal\tALREADY INSTALLED\n"
else
	cd tmp
	echo "Installation of Infernal (http://infernal.janelia.org/)"
	echo "Please wait..."

	echo -e "\nInstallation of Infernal (http://infernal.janelia.org/)\n">> $DEST/installation.log

	$FETCHER infernal-1.0.2.tar.gz "http://dl.dropbox.com/u/3753967/algorithms/infernal-1.0.2.tar.gz" >> $DEST/installation.log 2>&1

	tar -xzf infernal-1.0.2.tar.gz
	cd infernal-1.0.2
	./configure --prefix="$DEST"/infernal >> $DEST/installation.log 2>&1
	make clean > /dev/null 2>&1
	make >> $DEST/installation.log 2>&1
	make install >> $DEST/installation.log 2>&1

	if [ -f $DEST/infernal/bin/cmsearch ]
		then
    		RESULT=$RESULT"Infernal\tYES\n"
		else
    		RESULT=$RESULT"Infernal\tNO\n"
	fi
	
	echo 'export PATH=$PATH:$PYRNA_ALGORITHMS/infernal/bin' >> "$DEST/setmyenv"
fi

################# RNAVIEW ##########################

cd "$DEST"

if [ -f RNAVIEW/bin/rnaview ]
then
    RESULT=$RESULT"RNAVIEW\t\tALREADY INSTALLED\n"
else

	echo "Installation of RNAVIEW (http://ndbserver.rutgers.edu/services/download/)"
	echo "Please wait..."

	echo -e "\nInstallation of RNAVIEW (http://ndbserver.rutgers.edu/services/download/)\n">> $DEST/installation.log

	$FETCHER RNAVIEW.tar.gz "http://dl.dropbox.com/u/3753967/algorithms/RNAVIEW.tar.gz" >> $DEST/installation.log 2>&1

	tar -xzf RNAVIEW.tar.gz
	cd RNAVIEW
	make clean > /dev/null 2>&1
	make >> $DEST/installation.log 2>&1

	if [ -f $DEST/RNAVIEW/bin/rnaview ]
		then
    		RESULT=$RESULT"RNAVIEW\t\tYES\n"
		else
    		RESULT=$RESULT"RNAVIEW\t\tNO\n"
	fi
	
	mv ../RNAVIEW.tar.gz "$DEST"/tmp

	echo 'export RNAVIEW=$PYRNA_ALGORITHMS/RNAVIEW/' >> "$DEST/setmyenv"
	echo 'export PATH=$PATH:$RNAVIEW/bin' >> "$DEST/setmyenv"
fi

################# VIENNA PACKAGE ##########################

cd "$DEST"

if [ -d ViennaRNA/bin ]
then
    RESULT=$RESULT"RNAdistance\tALREADY INSTALLED\n"
	RESULT=$RESULT"RNAfold\t\tALREADY INSTALLED\n"
	RESULT=$RESULT"RNAplot\t\tALREADY INSTALLED\n"
else
	cd tmp
	echo "Installation of the Vienna RNA Package (http://www.tbi.univie.ac.at/RNA/)"
	echo "Please wait..."

	echo -e "\nInstallation of the Vienna RNA Package (http://www.tbi.univie.ac.at/RNA/)\n">> $DEST/installation.log
	
	$FETCHER ViennaRNA-1.8.5.tar.gz "http://dl.dropbox.com/u/3753967/algorithms/ViennaRNA-1.8.5.tar.gz" >> $DEST/installation.log 2>&1

	tar -xzf ViennaRNA-1.8.5.tar.gz
	cd ViennaRNA-1.8.5
	./configure --prefix="$DEST"/ViennaRNA >> $DEST/installation.log 2>&1
	make clean > /dev/null 2>&1
    make >> $DEST/installation.log 2>&1
	make install >> $DEST/installation.log 2>&1
	
	mv Utils/b2ct $DEST/ViennaRNA/bin

	if [ -d $DEST/ViennaRNA/bin ]
	then
    	RESULT=$RESULT"RNAdistance\tYES\n"
		RESULT=$RESULT"RNAfold\t\tYES\n"
		RESULT=$RESULT"RNAplot\t\tYES\n"
	else
    	RESULT=$RESULT"RNAdistance\tNO\n"
		RESULT=$RESULT"RNAfold\t\tNO\n"
		RESULT=$RESULT"RNAplot\t\tNO\n"
	fi

	echo 'export PATH=$PATH:$PYRNA_ALGORITHMS/ViennaRNA/bin' >> "$DEST/setmyenv"
fi

################# CONTRAFOLD ##########################

cd "$DEST"

if [ -f contrafold/src/contrafold ]
then
    RESULT=$RESULT"Contrafold\tALREADY INSTALLED\n"
else
	
	echo "Installation of Contrafold (http://contra.stanford.edu/contrafold/)"
	echo "Please wait..."

	echo -e "\nInstallation of Contrafold (http://contra.stanford.edu/contrafold/)\n" >> $DEST/installation.log

	$FETCHER contrafold.tar.gz "http://dl.dropbox.com/u/3753967/algorithms/contrafold_v2_02.tar.gz" >> $DEST/installation.log 2>&1

	tar -xzf contrafold.tar.gz
	cd contrafold/src
	make >> $DEST/installation.log 2>&1

	if [ -f $DEST/contrafold/src/contrafold ]
	then
    	RESULT=$RESULT"Contrafold\tYES\n"
	else
    	RESULT=$RESULT"Contrafold\tNO\n"
	fi
	
	mv ../../contrafold.tar.gz "$DEST"/tmp

	echo 'export PATH=$PATH:$PYRNA_ALGORITHMS/contrafold/src/' >> "$DEST/setmyenv"

fi

################# LOCARNA ##########################

cd "$DEST"

if [ -d locarna/bin ]
then
    RESULT=$RESULT"Locarna\t\tALREADY INSTALLED\n"
else
	cd tmp
	
	echo "Installation of Locarna (http://www.bioinf.uni-freiburg.de/Software/LocARNA/)"
	echo "Please wait..."

	echo -e "\nInstallation of Locarna (http://www.bioinf.uni-freiburg.de/Software/LocARNA/)\n" >> $DEST/installation.log

	$FETCHER locarna-1.7.2.tar.gz "http://dl.dropbox.com/u/3753967/algorithms/locarna-1.7.2.tar.gz" >> $DEST/installation.log 2>&1

	tar -xzf locarna-1.7.2.tar.gz
	cd locarna-1.7.2/
	./configure --prefix="$DEST"/locarna --with-vrna="$DEST"/ViennaRNA/ --without-perl --without-forester --without-kinfold >> $DEST/installation.log 2>&1
	make >> $DEST/installation.log 2>&1
	make install >> $DEST/installation.log 2>&1

	if [ -d $DEST/locarna/bin ]
	then
    	RESULT=$RESULT"Locarna\t\tYES\n"
	else
    	RESULT=$RESULT"Locarna\t\tNO\n"
	fi

	echo 'export PATH=$PATH:$PYRNA_ALGORITHMS/locarna/bin/' >> "$DEST/setmyenv"
fi

################# FOLDALIGN #######################

cd "$DEST"

#the archive of foldalign contains a linux executable, so testing the executable existence is not the good way to test if already installed.
#for now, no test....

echo "Installation of Foldalign (http://foldalign.ku.dk/software/index.html)"
echo "Please wait..."

echo -e "\nInstallation of Foldalign (http://foldalign.ku.dk/software/index.htm)\n" >> $DEST/installation.log

cd tmp

$FETCHER foldalign.2.1.1.tar.gz "http://dl.dropbox.com/u/3753967/algorithms/foldalign.2.1.1.tar.gz" >> $DEST/installation.log 2>&1

tar -xzf foldalign.2.1.1.tar.gz
cd foldalign.2.1.1/
make >> $DEST/installation.log 2>&1

cd ..

mv foldalign.2.1.1 ../foldalign.2.1.1

RESULT=$RESULT"Foldalign\tYES\n"

echo 'export PATH=$PATH:$PYRNA_ALGORITHMS/foldalign.2.1.1/bin/' >> "$DEST/setmyenv"

################# BLAST ##########################

cd "$DEST"

if [ -d blast-2.2.26/bin ]
then
    RESULT=$RESULT"Blast\t\tALREADY INSTALLED\n"
else
	
	echo "Installation of Blast... (ftp://ftp.ncbi.nlm.nih.gov/blast/)"
	echo "Please wait..."

	echo -e "\nInstallation of Blast... (ftp://ftp.ncbi.nlm.nih.gov/blast/)\n" >> $DEST/installation.log

	case "`uname`" in
            Darwin*)
                $FETCHER blast-2.2.26-universal-macosx.tar.gz "http://dl.dropbox.com/u/3753967/algorithms/blast-2.2.26-universal-macosx.tar.gz" >> $DEST/installation.log 2>&1
                tar -xzvf blast-2.2.26-universal-macosx.tar.gz >> $DEST/installation.log 2>&1
                mv blast-2.2.26-universal-macosx.tar.gz tmp/
                ;;
            *)
                $FETCHER blast-2.2.26-x64-linux.tar.gz "http://dl.dropbox.com/u/3753967/algorithms/blast-2.2.26-x64-linux.tar.gz" >> $DEST/installation.log 2>&1
                tar -xzvf blast-2.2.26-x64-linux.tar.gz >> $DEST/installation.log
                mv blast-2.2.26-x64-linux.tar.gz tmp/
                ;;
    esac

	if [ -d $DEST/blast-2.2.26/bin ]
	then
    	RESULT=$RESULT"Blast\t\tYES\n"
	else
    	RESULT=$RESULT"Blast\t\tNO\n"
	fi

	echo 'export PATH=$PATH:$PYRNA_ALGORITHMS/blast-2.2.26/bin/' >> "$DEST/setmyenv"
fi


####################### BLASTR #########################

cd "$DEST"

if [ -d blastR_package_V2.2 ]
then
    RESULT=$RESULT"BlastR\t\tALREADY INSTALLED\n"
else
    
    echo "Installation of BlastR... (http://www.tcoffee.org/Projects/blastr/)"
    echo "Please wait..."

    echo -e "\nInstallation of BlastR... (http://www.tcoffee.org/Projects/blastr/)\n" >> $DEST/installation.log
   
    $FETCHER blastR_package_V2.2.tar.gz "http://dl.dropboxusercontent.com/u/3753967/algorithms/blastR_package_V2.2.tar.gz" >> $DEST/installation.log 2>&1
    tar -xzvf blastR_package_V2.2.tar.gz >> $DEST/installation.log 2>&1
    mv blastR_package_V2.2.tar.gz tmp/

    if [ -d $DEST/blastR_package_V2.2 ]
    then
        RESULT=$RESULT"BlastR\t\tYES\n"
    else
        RESULT=$RESULT"BlastR\t\tNO\n"
    fi

    echo 'export PATH=$PATH:$PYRNA_ALGORITHMS/blastR_package_V2.2/scripts/' >> "$DEST/setmyenv"
fi


################# CLUSTALW ##########################

cd "$DEST"

if [ -d clustalw/bin/ ]
then
    RESULT=$RESULT"Clustalw\tALREADY INSTALLED\n"
else

	cd tmp
	
	echo "Installation of Clustalw... (http://www.clustal.org/)"
	echo "Please wait..."

	echo -e "\nInstallation of Clustalw... (http://www.clustal.org/)\n" >> $DEST/installation.log
	
    $FETCHER clustalw-2.1.tar.gz "http://dl.dropbox.com/u/3753967/algorithms/clustalw-2.1.tar.gz" >> $DEST/installation.log 2>&1
    tar -xzvf clustalw-2.1.tar.gz >> $DEST/installation.log 2>&1
    cd clustalw-2.1
    ./configure --prefix="$DEST"/clustalw >> $DEST/installation.log 2>&1
	make >> $DEST/installation.log 2>&1
	make install >> $DEST/installation.log 2>&1
                
	if [ -d $DEST/clustalw/bin/ ]
	then
    	RESULT=$RESULT"Clustalw\tYES\n"
	else
    	RESULT=$RESULT"Clustalw\tNO\n"
	fi

	echo 'export PATH=$PATH:$PYRNA_ALGORITHMS/clustalw/bin/' >> "$DEST/setmyenv"

fi

################# SAMTOOLS ##########################

cd "$DEST"

if [ -f samtools/samtools ]
then
    RESULT=$RESULT"SAMtools\tALREADY INSTALLED\n"
else

	cd tmp
	
	echo "Installation of SAMtools... (http://samtools.sourceforge.net/)"
	echo "Please wait..."

	echo -e "\nInstallation of SAMtools... (http://samtools.sourceforge.net/)\n" >> $DEST/installation.log
	
    $FETCHER samtools-0.1.18.tar.bz2 "http://dl.dropbox.com/u/3753967/algorithms/samtools-0.1.18.tar.bz2" >> $DEST/installation.log 2>&1
    bunzip2 -d samtools-0.1.18.tar.bz2 >> $DEST/installation.log 2>&1
    tar -xvf samtools-0.1.18.tar >> $DEST/installation.log 2>&1
    cd samtools-0.1.18/
	make >> $DEST/installation.log 2>&1
	cd ..
	mv samtools-0.1.18 ../samtools
                
	if [ -f $DEST/samtools/samtools ]
	then
    	RESULT=$RESULT"SAMtools\tYES\n"
	else
    	RESULT=$RESULT"SAMtools\tNO\n"
	fi

	echo 'export PATH=$PATH:$PYRNA_ALGORITHMS/samtools/' >> "$DEST/setmyenv"

fi

################# GotohScan ##########################

cd "$DEST"

if [ -f GotohScan_2.0-alpha/src/GotohScan2a ]
then
    RESULT=$RESULT"GotohScan\tALREADY INSTALLED\n"
else

    cd tmp
    
    echo "Installation of GotohScan... (http://www.bioinf.uni-leipzig.de/Software/GotohScan/README)"
    echo "Please wait..."

    echo -e "\nInstallation of GotohScan... (http://www.bioinf.uni-leipzig.de/Software/GotohScan/README)\n" >> $DEST/installation.log
    
    $FETCHER GotohScan_2.0-alpha.tar.gz "http://dl.dropbox.com/u/3753967/algorithms/GotohScan_2.0-alpha.tar.gz" >> $DEST/installation.log 2>&1
    tar -xzvf GotohScan_2.0-alpha.tar.gz >> $DEST/installation.log 2>&1
    cd GotohScan_2.0-alpha/

    make >> $DEST/installation.log 2>&1
    cd ..
    mv GotohScan_2.0-alpha ../GotohScan_2.0-alpha

    if [ -f $DEST/GotohScan_2.0-alpha/src/GotohScan2a ]
    then
        RESULT=$RESULT"GotohScan\tYES\n"
    else
        RESULT=$RESULT"GotohScan\tNO\n"
    fi

    echo 'export PATH=$PATH:$PYRNA_ALGORITHMS/GotohScan_2.0-alpha/src/' >> "$DEST/setmyenv"

fi

################# BCheck ##########################

#Bcheck needs RNABOB and Infernal

cd "$DEST"

if [ -f Bcheck-0.6/Bcheck ]
then
    RESULT=$RESULT"BCheck\t\tALREADY INSTALLED\n"
else

    cd tmp

    echo "Installation of BCheck... (http://rna.tbi.univie.ac.at/bcheck/)"
    echo "Please wait..."

    echo -e "\nInstallation of BCheck... (http://rna.tbi.univie.ac.at/bcheck/)\n" >> $DEST/installation.log
    
    #the Bcheck recovered from Dropbox has been modified to not check the rnabob and Infernal version
    $FETCHER Bcheck-0.6.tar.gz "http://dl.dropbox.com/u/3753967/algorithms/Bcheck-0.6.tar.gz" >> $DEST/installation.log 2>&1
    tar -xzvf Bcheck-0.6.tar.gz >> $DEST/installation.log 2>&1

    mv Bcheck-0.6 ../Bcheck-0.6

    if [ -f $DEST/Bcheck-0.6/Bcheck ]
    then
        RESULT=$RESULT"Bcheck\t\tYES\n"
    else
        RESULT=$RESULT"Bcheck\t\tNO\n"
    fi

    echo 'export PATH=$PATH:$PYRNA_ALGORITHMS/Bcheck-0.6/' >> "$DEST/setmyenv"

fi

################# RNABOB ##########################

cd "$DEST"

if [ -f rnabob-2.2.1/rnabob ]
then
    RESULT=$RESULT"RNABOB\t\tALREADY INSTALLED\n"
else

    cd tmp
    
    echo "Installation of RNABOB... (http://selab.janelia.org/software.html)"
    echo "Please wait..."

    echo -e "\nInstallation of RNABOB... (http://selab.janelia.org/software.html)\n" >> $DEST/installation.log
    
    $FETCHER rnabob.tar.gz "http://dl.dropbox.com/u/3753967/algorithms/rnabob.tar.gz" >> $DEST/installation.log 2>&1
    tar -xzvf rnabob.tar.gz >> $DEST/installation.log 2>&1
    cd rnabob-2.2.1/
    make >> $DEST/installation.log 2>&1
    cd ..
    mv rnabob-2.2.1 ../rnabob-2.2.1

    if [ -f $DEST/rnabob-2.2.1/rnabob ]
    then
        RESULT=$RESULT"RNABOB\t\tYES\n"
    else
        RESULT=$RESULT"RNABOB\t\tNO\n"
    fi

    echo 'export PATH=$PATH:$PYRNA_ALGORITHMS/rnabob-2.2.1/' >> "$DEST/setmyenv"

fi

################# SNOREPORT ##########################

cd "$DEST"

if [ -f SnoReport1.0/snoReport ]
then
    RESULT=$RESULT"snoReport\tALREADY INSTALLED\n"
else

    cd tmp
    
    echo "Installation of snoReport... (http://www.bioinf.uni-leipzig.de/Software/snoReport/)"
    echo "Please wait..."

    echo -e "\nInstallation of snoReport... (http://www.bioinf.uni-leipzig.de/Software/snoReport/)\n" >> $DEST/installation.log
    
    $FETCHER SnoReport1.0.tgz "http://dl.dropboxusercontent.com/u/3753967/algorithms/SnoReport1.0.tgz" >> $DEST/installation.log 2>&1
    tar -xzvf SnoReport1.0.tgz >> $DEST/installation.log 2>&1
    cd SnoReport1.0/
    #SnoReport needs the Vienna package in the PATH. How to do it?? TO DO!!!!!
    ./install.sh >> $DEST/installation.log 2>&1
    cd ..
    mv SnoReport1.0 ../SnoReport1.0

    if [ -f $DEST/SnoReport1.0/snoReport ]
    then
        RESULT=$RESULT"snoReport\tYES\n"
    else
        RESULT=$RESULT"snoReport\tNO\n"
    fi

    echo 'export SNOREPORT=$PYRNA_ALGORITHMS/SnoReport1.0/models' >> "$DEST/setmyenv"
    echo 'export PATH=$PATH:$PYRNA_ALGORITHMS/SnoReport1.0/' >> "$DEST/setmyenv"

fi

################# SNOSCAN ##########################

cd "$DEST"

if [ -f snoscan-0.9b/snoscan ]
then
    RESULT=$RESULT"Snoscan\t\tALREADY INSTALLED\n"
else

    cd tmp
    
    echo "Installation of Snoscan... (http://lowelab.ucsc.edu/snoscan/)"
    echo "Please wait..."

    echo -e "\nInstallation of Snoscan... (http://lowelab.ucsc.edu/snoscan/)\n" >> $DEST/installation.log
    
    $FETCHER snoscan.tar.gz "https://dl.dropboxusercontent.com/u/3753967/algorithms/snoscan.tar.gz" >> $DEST/installation.log 2>&1
    tar -xzvf snoscan.tar.gz >> $DEST/installation.log 2>&1
    cd snoscan-0.9b/squid-1.5j
    make >> $DEST/installation.log 2>&1
    cd ..
    make >> $DEST/installation.log 2>&1
    cd ..
    mv snoscan-0.9b ../snoscan-0.9b

    if [ -f $DEST/snoscan-0.9b/snoscan ]
    then
        RESULT=$RESULT"Snoscan\t\tYES\n"
    else
        RESULT=$RESULT"Snoscan\t\tNO\n"
    fi

    echo 'export PATH=$PATH:$PYRNA_ALGORITHMS/snoscan-0.9b/' >> "$DEST/setmyenv"
fi

################# SNOGPS ##########################

cd "$DEST"

if [ -f snoGPS-0.2/src/pseudoU_test ]
then
    RESULT=$RESULT"SnoGPS\t\tALREADY INSTALLED\n"
else

    cd tmp
    
    echo "Installation of SnoGPS... (http://lowelab.ucsc.edu/snoGPS/)"
    echo "Please wait..."

    echo -e "\nInstallation of SnoGPS... (http://lowelab.ucsc.edu/snoGPS/)\n" >> $DEST/installation.log
    
    $FETCHER snoGPS-0.2.tar.gz "https://dl.dropboxusercontent.com/u/3753967/algorithms/snoGPS-0.2.tar.gz" >> $DEST/installation.log 2>&1
    tar -xzvf snoGPS-0.2.tar.gz >> $DEST/installation.log 2>&1
    cd snoGPS-0.2/src
    make >> $DEST/installation.log 2>&1
    cd ../..
    mv snoGPS-0.2 ../snoGPS-0.2

    if [ -f $DEST/snoGPS-0.2/src/pseudoU_test ]
    then
        RESULT=$RESULT"SnoGPS\t\tYES\n"
    else
        RESULT=$RESULT"SnoGPS\t\tNO\n"
    fi

    echo 'export PATH=$PATH:$PYRNA_ALGORITHMS/snoGPS-0.2/src' >> "$DEST/setmyenv"
fi

################# RNAMOTIF ##########################

cd "$DEST"

if [ -f rnamotif-3.0.7/src/rnamotif ]
then
    RESULT=$RESULT"RNAMotif\tALREADY INSTALLED\n"
else

    cd tmp
    
    echo "Installation of RNAMotif... (http://casegroup.rutgers.edu/casegr-sh-2.5.html)"
    echo "Please wait..."

    echo -e "\nInstallation of RNAMotif... (http://casegroup.rutgers.edu/casegr-sh-2.5.html)\n" >> $DEST/installation.log
    
    $FETCHER rnamotif-3.0.7.tar.gz "https://dl.dropboxusercontent.com/u/3753967/algorithms/rnamotif-3.0.7.tar.gz" >> $DEST/installation.log 2>&1
    tar -xzvf rnamotif-3.0.7.tar.gz >> $DEST/installation.log 2>&1
    cd rnamotif-3.0.7
    make >> $DEST/installation.log 2>&1
    cd ..
    mv rnamotif-3.0.7 ../rnamotif-3.0.7

    if [ -f $DEST/rnamotif-3.0.7/src/rnamotif ]
    then
        RESULT=$RESULT"RNAMotif\tYES\n"
    else
        RESULT=$RESULT"RNAMotif\tNO\n"
    fi

    echo 'export PATH=$PATH:$PYRNA_ALGORITHMS/rnamotif-3.0.7/src' >> "$DEST/setmyenv"
fi

################# MUSCLE ##########################

cd "$DEST"

if [ -f muscle3.8.31/src/muscle ]
then
    RESULT=$RESULT"MUSCLE\tALREADY INSTALLED\n"
else

    cd tmp
    
    echo "Installation of MUSCLE... (http://www.drive5.com/muscle/)"
    echo "Please wait..."

    echo -e "\nInstallation of MUSCLE... (http://www.drive5.com/muscle/)\n" >> $DEST/installation.log
    
    $FETCHER muscle3.8.31_src.tar.gz "https://dl.dropboxusercontent.com/u/3753967/algorithms/muscle3.8.31_src.tar.gz" >> $DEST/installation.log 2>&1
    tar -xzvf muscle3.8.31_src.tar.gz >> $DEST/installation.log 2>&1
    cd muscle3.8.31/src
    make >> $DEST/installation.log 2>&1
    cd ../..
    mv muscle3.8.31 ../muscle3.8.31

    if [ -f $DEST/muscle3.8.31/src/muscle ]
    then
        RESULT=$RESULT"MUSCLE\tYES\n"
    else
        RESULT=$RESULT"MUSCLE\tNO\n"
    fi

    echo 'export PATH=$PATH:$PYRNA_ALGORITHMS/muscle3.8.31/src' >> "$DEST/setmyenv"
fi

################# RNAZ ##########################

cd "$DEST"

if [ -f RNAz-2.1/rnaz/RNAz ]
then
    RESULT=$RESULT"RNAz\tALREADY INSTALLED\n"
else

	cd tmp
    
    echo "Installation of RNAz... (http://www.tbi.univie.ac.at/~wash/RNAz/)"
    echo "Please wait..."

    echo -e "\nInstallation of RNAz... (http://www.tbi.univie.ac.at/~wash/RNAz/)\n" >> $DEST/installation.log
    
    $FETCHER RNAz-2.1.tar.gz "https://dl.dropboxusercontent.com/u/3753967/algorithms/RNAz-2.1.tar.gz" >> $DEST/installation.log 2>&1
    tar -xzvf RNAz-2.1.tar.gz >> $DEST/installation.log 2>&1

    cd RNAz-2.1

    ./configure >> $DEST/installation.log 2>&1

    make >> $DEST/installation.log 2>&1

    cd ..

    mv RNAz-2.1 ../RNAz-2.1

    if [ -f $DEST/RNAz-2.1/rnaz/RNAz ]
    then
        RESULT=$RESULT"RNAz\tYES\n"
    else
        RESULT=$RESULT"RNAz\tNO\n"
    fi

    echo 'export PATH=$PATH:$PYRNA_ALGORITHMS/RNAz-2.1/rnaz/' >> "$DEST/setmyenv"
fi

######## END OF INSTALLATION #########################

cd "$DEST"

echo -e "\nResult of the RNA Algorithms installation:"
echo -e   "-------------------------------------------\n"
echo -e "The following algorithms have been installed in the directory $DEST :\n"
echo -e $RESULT
echo -e "IMPORTANT: now you have to add the following instruction in the configuration file of your shell: source $DEST/setmyenv\n"
echo "For any problem, check the installation.log file located in $DEST."

echo -e "\nFor any questions or remarks, please contact f.jossinet@ibmc-cnrs.unistra.fr"