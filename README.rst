TS-PPTIS (Transition State PPTIS)
=================================

Version 2.0
-----------

New python implementation of TS-PPTIS for Gromacs 5.X and Plumed 2.X
Based on scripts by G. Saladino and J. Juraszek @2014-2016 

For reference see 

Juraszek J, Saladino G, van Erp TS, and Gervasio FL, "Efficient numerical reconstruction of protein folding kinetics with partial path sampling and pathlike variables." Phys Rev Lett. 2013 Mar 8;110(10):108106.

.. The package is now available on PyPI, to retrieve it just type ``pip install tspptisM`` or download it from here and install with ``python setup.py install``.

The code has been tested on Python 2.7 and Python 3.3

Dependencies
------------

- Numpy 1.11.0;
- MDTraj 1.8.0.

Utility Scripts Usage
---------------------

A number of utility scripts are provided in the utils folder. Here is an example of a typical workflow.

First run the setup script to prepare all the necessary files, by giving as input information on the system (topology, coordinates, etc...) and the location of your local gromacs installation::

    tsSetup.py -top topol.top -gro system.gro -mdp md.mdp -gmx /usr/bin/gmx
    -win windowsList.dat -xtc traj.xtc -col COLVAR

The script will automatically store the information in a file (default: init.info) that is required to run the next set of scripts.

A text file ``windowsList.dat`` containing information on the windows to set up needs to be included. The format should be::

    0.0,0.5,0.75
    0.65,1.1,1.5
    ...

with each line containing a floating point representing the position of the left, centre and right window bounds, separated by a comma.

With ``tsSetup.py`` the tspptis simulation is ready to run. 
To finalize a single run and setup the next iteration, ``tsFinalize.py`` and ``tsSetRun.py`` can be used for each window, all they need is the ``init.info`` and the tspptis window folder::

    tsFinalize.py -info init.info -fold pptis00
    tsSetRun.py -info init.info -fold pptis00

A standard analysis worfklow is available in ``tsAnalysis.py``, which requiresthe path to the tspptis windows folders, a free energy profile (in Plumed 2 format), and an estimate of the TS. If none provided, the FES maximum will be taken as TS::

    tsAnalysis.py -trans 1.55 -fold ./ -fes fes.dat 


Library Usage
-------------

Alternatively, you can use the library functions in a python script:, which gives you more control on the tspptis parameters::

    #Import the library
    import tspptis as tps

    #Initialise the system
    ts = tsSetup('./topol.top',
                 './system.gro',
                 './md.mdp',
                 gmx='/usr/bin/gmx')
        
    #Initialise a single window
    ts.initWindow('./pptis10',
                  [1.5,1.7,1.9],
                 './traj_fixed.xtc',
                  './COLVAR',
                   overwrite=True)

    #Setup the simulation for a specific window
    ts.setUpTPS('../testfiles/pptis10')

    #After the simulation terminates, finalize
    ts.finalizeTPS('../testfiles/pptis10')

    #Repeat the two steps above if neede, then proceed to the analysis
    tsa = tsAnalysis('./pptisoutput')

    #Get probabilities and crossings analysis at a given TS of 1.25
    tsa.getProbabilities()
    tsa.getCrossings(1.25)
    
    #Finally extract the Rates, after converting the fes if needed
    fesList=plumed2List('./outputs/fes.dat')
    tsa.getRates(fesList,valTS=1.25)

    #If desired, check the Memory Loss Assumption
    tsa.endPointVel()
    tsa.checkMLA()


Contributors
-------------

Federico Comitani (f.comitani@ucl.ac.uk)
Giulio Mattedi (g.mattedi.16@ucl.ac.uk)

.. When using TS-PPTIS,  please cite
.. G.Mattedi, F.Comitani ... F.L.Gervasio ...


