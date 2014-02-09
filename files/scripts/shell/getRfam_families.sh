#!/bin/bash

DOWNLOAD_DIR=$1

VERSION=$2

LINK="ftp://ftp.sanger.ac.uk/pub/databases/Rfam/$VERSION/database_files/"

cd $DOWNLOAD_DIR
wget -qO rfam.txt.gz $LINK/rfam.txt.gz
gzip -d rfam.txt.gz

wget -qO pdb_rfam_reg.txt.gz $LINK/pdb_rfam_reg.txt.gz
gzip -df pdb_rfam_reg.txt.gz