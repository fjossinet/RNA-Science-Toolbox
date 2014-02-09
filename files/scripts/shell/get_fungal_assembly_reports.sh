#!/bin/bash

DOWNLOAD_DIR=$1

cd $1

wget --retr-symlinks -qr "ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/Eukaryotes/fungi/"