#!/bin/bash

DOWNLOAD_DIR=$1

DOWNLOAD_LINK=$2

cd $DOWNLOAD_DIR
wget -qO taxonomy.txt.gz $DOWNLOAD_LINK
gzip -d taxonomy.txt.gz