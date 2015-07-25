#!/bin/bash

TIMEZONE="Europe/Paris" #you can change that to fit your needs

echo "[vagrant provisioning] Setting timezone..."
echo $TIMEZONE | sudo tee /etc/timezone
sudo dpkg-reconfigure --frontend noninteractive tzdata
