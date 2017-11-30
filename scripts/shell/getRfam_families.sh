  #!/bin/bash

DOWNLOAD_DIR=$1

VERSION=$2

LINK="ftp://ftp.ebi.ac.uk/pub/databases/Rfam/$VERSION/database_files/"

cd $DOWNLOAD_DIR
wget -qO family.txt.gz $LINK/family.txt.gz
gzip -d family.txt.gz

#wget -qO pdb_rfam_reg.txt.gz $LINK/pdb_rfam_reg.txt.gz
#gzip -df pdb_rfam_reg.txt.gz
