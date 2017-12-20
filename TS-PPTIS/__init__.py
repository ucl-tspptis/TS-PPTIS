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

        print SectionDelimiter("INITIALISATION")

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

        print "Initialising window:\t\t", path

        # Check if folder exists and if overwriting is allowed
        if os.path.isdir(path):
            if overwrite:
                print "Folder exists, overwriting."
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

    def setUpTPS(self,path):
        if path[-1] != '/': path+= '/'
        """Check if the folder is a ts-pptis window."""
        if os.path.isfile(path+'window.cfg'):
            print SectionDelimiter("SETUP RUN")
        else:
            sys.exit('Error: the folder does not seem to be a TS-PPTIS window')

        print 'Setting up run in:\t\t', path

        continuation = False

        tpsAccHandle = open(path+'tps_acc.log','a+')
        tpsAccLines = sum([1 for line in tpsAccHandle])

        if tpsAccLines > 1: continuation = True

        print 'First run:\t\t\t',not continuation

        config = parseConfig(path+'window.cfg')
        print 'Interfaces:\t\t\t', config['interfaces']


        if not continuation:
            pathLength = len(self.trajData)
            tpsAccHandle.write(('0     0000       -          initial    1 '
                                '{:>6} 1.0000   A  B  1   0.00       0     -      1 1 1 1'
                                '\n').format(pathLength))

            # TODO:
            # WHATEVER HAPPENS IF IT IS THE FIRST RUN
        else:
            # TODO:
            # WHATEVER HAPPENS IF IT IS NOT THE FIRST RUN
            pass
            # in both cases, roughly, **I believe**:
            # determining previous path length, selecting random frame,
            # determining the LPF of the selected frame (whether BW or FW, use PAR file)
            # setting tmap (maximum duration of current run) with
            # path_length/random(),
            # extracting random frame, randomising velocities, generating a reversed TPR
            # for the BW trajectory
        print "Path length:\t\t\t", pathLength

        tpsAccHandle.close()


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

    ts.setUpTPS('../testfiles/pptis10')

if __name__ == "__main__":

    print("Running test...\n")
    testAll()
    print("\nDone!")
