#!/bin/bash

DOWNLOAD_DIR=$1

DOWNLOAD_LINK=$2

cd $DOWNLOAD_DIR
wget -qO genome_entry.txt.gz $DOWNLOAD_LINK
gzip -d genome_entry.txt.gz