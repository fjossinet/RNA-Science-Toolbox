#!/bin/bash

echo "[vagrant provisioning] Configuring the system..."
TIMEZONE="Europe/Paris" #you can change that to fit your needs
echo $TIMEZONE | sudo tee /etc/timezone
sudo dpkg-reconfigure --frontend noninteractive tzdata
echo "source /vagrant/conf/bash_login" >> /home/vagrant/.bash_login
sudo chown vagrant.vagrant /home/vagrant/.bash_login
echo "source /vagrant/conf/bashrc" >> /home/vagrant/.bashrc
sudo chown vagrant.vagrant /home/vagrant/.bash_login

echo "[vagrant provisioning] Installing vim..."
sudo apt-get -y install vim

echo "[vagrant provisioning] Installing MongoDB..."
sudo apt-get -y install mongodb-server

echo "[vagrant provisioning] Installing fabric..."
sudo apt-get -y install fabric

echo "[vagrant provisioning] Installing Miniconda..."
wget -q https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh
chmod u+x Miniconda-latest-Linux-x86_64.sh
./Miniconda-latest-Linux-x86_64.sh -b -p /home/vagrant/miniconda
rm Miniconda-latest-Linux-x86_64.sh
chown -Rf vagrant.vagrant /home/vagrant/miniconda
