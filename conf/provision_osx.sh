#!/bin/bash
#launch this script at first to install what is necessary for OSX: Homebrew, Anaconda distribution, MongoDB, Fabric tool

echo "Installing HomeBrew"
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

echo "Installing Anaconda"
wget "https://3230d63b5fc54e62148e-c95ac804525aac4b6dba79b00b39d1d3.ssl.cf1.rackcdn.com/Anaconda2-4.0.0-MacOSX-x86_64.sh"
chmod u+x ./Anaconda2-4.0.0-MacOSX-x86_64.sh
./Anaconda2-4.0.0-MacOSX-x86_64.sh

echo "Installing Fabric"

echo "Installing MongoDB"
