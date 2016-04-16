#!/bin/bash

echo "[vagrant provisioning] Configuring the system..."
TIMEZONE="Europe/Paris" #you can change that to fit your needs
echo $TIMEZONE | sudo tee /etc/timezone
sudo dpkg-reconfigure --frontend noninteractive tzdata

if grep "source /vagrant/conf/bash_login" /home/vagrant/.bash_login > /dev/null
then
  echo ""
else
  echo "source /vagrant/conf/bash_login" >> /home/vagrant/.bash_login
  sudo chown vagrant.vagrant /home/vagrant/.bash_login
fi

if grep "source /vagrant/conf/bashrc" /home/vagrant/.bashrc > /dev/null
then
  echo ""
else
  echo "source /vagrant/conf/bashrc" >> /home/vagrant/.bashrc
  sudo chown vagrant.vagrant /home/vagrant/.bash_login
fi

if ! [ -x "$(command -v g++)" ]
then
  echo "[vagrant provisioning] Installing g++..."
  sudo apt-get -y install g++
fi

if ! [ -x "$(command -v vim)" ]
then
  echo "[vagrant provisioning] Installing vim..."
  sudo apt-get -y install vim
fi

if ! [ -x "$(command -v mongod)" ]
then
  echo "[vagrant provisioning] Installing MongoDB..."
  sudo apt-get -y install mongodb-server
fi

#to access the mongodb server from any interface
if grep "bind_ip = 127.0.0.1" /etc/mongodb.conf > /dev/null
then
  echo "[vagrant provisioning] Changing the MongoDB bind IP..."
  sudo sed -i 's/bind_ip = 127.0.0.1/bind_ip = 0.0.0.0/' /etc/mongodb.conf
  sudo service mongodb stop
  sudo service mongodb start
fi

if ! [ -x "$(command -v fab)" ]
then
  echo "[vagrant provisioning] Installing fabric..."
  sudo apt-get -y install fabric
fi

if ! [ -d "/home/vagrant/miniconda" ]
then
  echo "[vagrant provisioning] Installing Miniconda..."
  wget -q https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh
  chmod u+x Miniconda-latest-Linux-x86_64.sh
  ./Miniconda-latest-Linux-x86_64.sh -b -p /home/vagrant/miniconda
  rm Miniconda-latest-Linux-x86_64.sh
  chown -Rf vagrant.vagrant /home/vagrant/miniconda
fi
