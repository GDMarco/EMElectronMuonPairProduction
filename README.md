# EMElectronMuonPairProduction

Plugin for the CRPropa simulation code to include the production of muon pairs from the interaction between highly-energetic electrons with background radiation fields.

## Structure
The header of the `EMElectronMuonPairProduction` class is defined in `EMElectronMuonPairProduction.h`. The C++ file `EMElectronMuonPairProduction.cc` contains the interaction class functions and the function to produce the secondaries. The file `EMElectronMuonPairProduction.i` is the SWIG interface for python usage. This structure follows the usual CRPropa modules for electromagnetic processes, e.g. see [CRPropa EMPairProduction](https://github.com/CRPropa/CRPropa3/blob/master/include/crpropa/module/EMPairProduction.h).
 
The python script `testPlugin.py` tests the installation (see below) of the example plugin as presented here.

## Interaction tables

An updated version of the [CRPropa3-data](https://github.com/CRPropa/CRPropa3-data) repository has been develped to include the interaction rates for this plugin, see [EMCascadePlugins-data](https://github.com/GDMarco/CRPropa3-data/tree/EMCascadePlugins-data). For the production of the secondaries, the inelasticieties are computed in the [MUNHECA code](https://github.com/afesmaeili/MUNHECA/tree/main).

# Installation
For the installation of the plugin you need a running CRPropa version (see [installation documentation](https://crpropa.github.io/CRPropa3/pages/Installation.html)).
This is done analogously to the installation of CRPropa. We recommend to activate the same virtual python environment that you use to run CRPropa.

First create a build folder within the plugin's directory and move there.

    mkdir build && cd build/

Now you can run cmake to configure your project:

    ccmake ..

At this step you have to set the installation path and the path to your swig interface of the current CRPropa installation (if it is not found by cmake).

After configuration (press c) and generation (press g) you can now build and install your plugin

    make install

## optional testing
Now you can run the python test script. 

    python ../testPlugin.py
