#!/bin/bash

DOWNLOAD_DIR=$1

DOWNLOAD_LINK=$2

cd $DOWNLOAD_DIR
wget -qO cdd.tar.gz $DOWNLOAD_LINK
tar-xzvf -d cdd.tar.gz