"""
TS-PPTIS v2.0

New implementation of TS-PPTIS for Gromacs 5.X and Plumed 2.X
Based on scripts by G. Saladino and J. Juraszek @2014-2016

For reference see

Juraszek J, Saladino G, van Erp TS, and Gervasio FL, "Efficient numerical reconstruction of protein folding kinetics with partial path sampling and pathlike variables." Phys Rev Lett. 2013 Mar 8;110(10):108106.

F. Comitani: f.comitani@ucl.ac.uk
G. Mattedi: giulio.mattedi.16@ucl.ac.uk
@2017-2018

"""


from __future__ import division

import sys
import os
import shutil # for removing directory tree
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
from tools import *

# Print debug info
debug = True

class tsSetup:
    """ Standard TS-PPTIS setup class. """

    def __init__(self, top, gro, traj, colvar, windows, par, ndx='', gmx='$GMX'):
        """Initialise TS-PPTIS.
        Args:
                top (string): path to topology file .top
                gro (string): path to structure file .gro
                traj (string): path to initial trajectory .trr/.xtc
                colvar (string): path to the colvar file of the input trajectory
                windows (string): path to text file containing information
                on the windows in the format left:center:right
                ndx (string, optional): path to groups definition file .ndx
                gmx (string, optional): path to the local gromacs executable.

        """
        if debug: print '# ' + timestamp()

        """Check and load trajectory data."""
        try:
            self.top = top
            self.gro = gro
            self.traj = traj
            self.trajData = md.load(traj, top=gro)
            print 'Topology and trajectory files:\tOK'
        except:
            sys.exit(
                'Error: invalid input topology/trajectory files ' + gro + ' ' + traj)

        """Check and load colvar file."""
        if os.path.isfile(colvar):
            self.colvar = colvar
            self.trajCV = parseColvar(colvar) # Might not be needed at this point. G.
        else:
            print 'COLVAR file:\t\t\tnot found'

        """Check and load windows file."""
        try:
            self.winList = parseWindows(windows)
            print 'PPTIS windows:\t\t\tOK'
        except:
            sys.exit('Error: invalid windows file ' + windows)

        """Check ndx file."""
        if os.path.isfile(ndx):
            self.ndx = ndx
            print "ndx file:\t\t\tOK"
        else:
            self.ndx = ''
            print 'nxd file:\t\t\tnot found'

        """Check for par file, if not foung generate from COLVAR"""
        if os.path.isfile(par):
            print "PAR file:\t\t\tOK"
        else:
            print "PAR file:\t\t\tnot found...generating it...",
            generatePar(colvar,par)
            print "OK"
        self.par = par

        """Check local gromacs installation."""
        self.gmx = findExe(gmx)
        if self.gmx != None:
            print 'Gromacs installation:\t\tOK'
        else:
            sys.exit('Error : invalid gmx path ' + gmx + '\n' +
                     'Make sure to have a working version of gromacs 5.X installed!')


    def initWindow(self, path, window, overwrite=False):
        """Initialize a window

        Args:
            path (string): path of the window directory
            window (list): interfaces of the window
            overwrite (bool): whether to overwrite existing folder

        """

        print "Initialising window", path

        # Check if folder exists and if overwriting is allowed
        if os.path.isdir(path):
            if overwrite:
                print "\tFolder exists, overwriting."
                shutil.rmtree(path)
            else:
                sys.exit('Refusing to overwrite directory')

        # Check length of interface list
        if len(window) != 3:
            sys.exit('Wrong number of elements as window interfaces')

        # Add trailing / to path if not provided
        if path[-1] != '/': path += '/'

        # Create the list of folders
        windowTree = [path, path+'data/', path+'run/', path+'temp/']
        for folder in windowTree:
            os.makedirs(folder)

        # In the data/ directory create symlinks to the initial
        # trajectory data
        os.symlink(os.path.abspath(self.traj),
                path+'data/00000.'+self.traj.split('.')[-1])
        os.symlink(os.path.abspath(self.colvar),
                path+'data/00000.cv')
        os.symlink(os.path.abspath(self.par),
                path+'data/00000.par')

        # Initialize a config file. Can be useful for storing paths and
        # various configurations. See the config file in the old implementation
        with open(path+'window.cfg','w') as handle:
            initText = '# %s\ninterfaces = %s\n' % (timestamp(), ':'.join(map(str,window)))
            handle.write(initText)


def testAll():

    # Test initialisation
    ts = tsSetup('../testfiles/topol.top',
                 '../testfiles/system.gro',
                 '../testfiles/traj_fixed.xtc',
                 '../testfiles/COLVAR',
                 '../testfiles/windows.dat',
                 '../testfiles/traj.par',
                 gmx='gmx')

    # Test extractFrame
#    print extractFrame([0,1.2],
#            trajFile='../testfiles/traj_fixed_skipped.xtc',
#            topFile='../testfiles/system.gro',
#            colvarFile='../testfiles/COLVAR',
#            outFile='../testfiles/out.pdb',trajStride=10,colvarStride=1)

    # Test initialising a window
    ts.initWindow('../testfiles/pptis10',[0,1,2], overwrite=True)

if __name__ == "__main__":

    print("Running test...\n")
    testAll()
    print("\nDone!")
