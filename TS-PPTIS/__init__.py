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
import subprocess
import numpy as np

import matplotlib.pyplot as plt

import mdtraj as md

from tools import *


class tsSetup:
    """ Standard TS-PPTIS setup class. """

    def __init__(self, top, gro, traj, colvar, windows, ndx='', gmx='$GMX'):
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
            self.trajCV = parseColvar(colvar)
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
        else:
            self.ndx = ''
            print 'nxd file:\t\t\tnot found'

        """Check local gromacs installation."""
        self.gmx = findExe(gmx)
        if self.gmx != None:
            print 'Gromacs installation:\t\tOK'
        else:
            sys.exit('Error : invalid gmx path ' + gmx + '\n' +
                     'Make sure to have a working version of gromacs 5.X installed!')

    def genConfig(self):
        """Generate a dictionary with paths and configs. DO WE REALLY WANT IT?"""
        self.config = {'top': self.top, 'gro': self.gro, 'traj': self.traj, 'colvar': self.colvar,
                       'winlist': self.winList, 'ndx': self.ndx, 'gmx': self.gmx}

        return self.config


    def findNearest(self, array, value):
        """ Find nearest value in an array."""
        idx = (np.abs(array-value)).argmin()
        return array[idx]

    def extractFrame(self, cvValue, trajFile='',  top='', colvarFile='', outFile='out.pdb', trajStride=1, colvarStride=1, save=True):
        """Extract a frame with a given CV value from a trajectory. If files are not specified
        it will load self.trajData e self.trajCV. The function looks for the frame with the closest CV values

        Args:
            cvValue (float/list/tuple): CV value or CV range
            trajFile (string): trajectory file.
            top (string): topology for mdtraj.
            colvarFile (string): COLVAR file
            outFile (string): output file name
            trajStride (int): trajectory stride
            colvarStride (int): COLVAR stride
            save (bool): Whether to save the frame or not

        Returns:
            frames (int): frame number
            frameCV (float): frame CV value

        """

        # if trajFile is not provided load traj from self.trajData
        if not trajFile:
            traj = self.trajData
        else:
            traj = md.load(trajFile,top=top)

        # if colvarFile is not provided load COLVAR from self.trajCV
        if not colvarFile:
            colvar = self.trajCV
        else:
            colvar = parseColvar(colvarFile)

        # Subsample trajectory or COLVAR depending on who has higher stride
        # in order to have (almost) equal number of frames and colvar lines

        # NEEDS MORE WORK
        if trajStride > colvarStride:
            stride = int(trajStride/colvarStride)
            colvar = colvar[::stride]
        else:
            stride = int( colvarStride/trajStride)
            traj = traj[::stride]

        # IF cvValue is a float look for the value, else for the inclusive range
        if type(cvValue) == float:
            framesCV = self.findNearest(colvar[:,1],cvValue)
            frames = np.where(colvar[:,1] == frameCV)[0]
        else:
            cvEntries = np.where(np.logical_and((colvar[:,1] >= cvValue[0]),(colvar[:,1] <= cvValue[1])))[0]
            framesCV = colvar[cvEntries,1]
            frames = cvEntries

        if save: traj[frames].save(outFile)

        return frames, framesCV


def testAll():

    ts = tsSetup('../testfiles/topol.top',
                 '../testfiles/system.gro',
                 '../testfiles/traj_fixed.xtc',
                 '../testfiles/COLVAR',
                 '../testfiles/windows.dat',
                 gmx='gmx')
    print ts.extractFrame([0,1.2], trajFile='../testfiles/traj_fixed_skipped.xtc',top='../testfiles/system.gro',
            outFile='../testfiles/out.pdb',trajStride=10,colvarStride=1)

if __name__ == "__main__":

    print("Running test...\n")
    testAll()
    print("\nDone!")
