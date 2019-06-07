#!/bin/bash

INSTALL_DIRECTORY=$PWD
#compile for spin_pm1
cd $INSTALL_DIRECTORY/code_directACE_spin_pm1/
make
cp $INSTALL_DIRECTORY/code_directACE_spin_pm1/directACE_spin_pm1.out $INSTALL_DIRECTORY
cd $INSTALL_DIRECTORY

#compile for spin_01
cd $INSTALL_DIRECTORY/code_directACE_spin_01/
make
cp $INSTALL_DIRECTORY/code_directACE_spin_01/directACE_spin01.out $INSTALL_DIRECTORY
cd $INSTALL_DIRECTORY
