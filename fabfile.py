import os, sys
from os.path import expanduser
from fabric.context_managers import cd
from fabric.contrib.console import confirm
from fabric.operations import prompt
from fabric.operations import run
from fabric.colors import green
from fabric.api import *

def all():
    if not confirm("This will install and configure all the dependencies. Do you wish to continue?") :
        print "Bye!"
        sys.exit()
    update()
    algorithms()
    web_stack()
    website()

def update():
    print(green("Updating the package list..."))
    local('sudo apt-get update -qq')

def algorithms():
    print(green("Installation of the RNA algorithms..."))
    algorithms_home = prompt('Where do you want to install your RNA algorithms?', default=os.path.join(expanduser('~'), 'algorithms'), validate=r'^'+expanduser('~')+'/.+/?$')
    #first the dependencies needed for recompilation

def webstack():
    print(green("Installation of a full web stack..."))
    print(green("Node.js..."))
    local("sudo apt-get -y install nodejs npm")
    local("sudo ln -s `which nodejs` /usr/bin/node")
    print(green("Bower..."))
    local("sudo npm install -g bower")

def website():
    print(green("Installation of the dependencies for the website..."))
    local('cd website ; bower --config.interactive=false install')
