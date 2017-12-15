"""
TS-PPTIS tools

A set of utility functions for TS-PPTIS

F. Comitani: f.comitani@ucl.ac.uk
G. Mattedi: giulio.mattedi.16@ucl.ac.uk
@2017-2018

"""

from __future__ import division

import sys
import os
import subprocess
import numpy as np


def parseWindows(nameFile):

    """Parse a text file containing information on the TS-PPTIS windows.

    Args:
        windows (string): path to  text file containing information
        on the windows in the format left:center:right

    Returns:
        outList (nested list, float): list containing left, center and right value
        of each window.
    """

    outList = []
    data =open(nameFile,"r")

    for line in data.readlines():
        outList.append([])
        for i in line.split(':'):
            outList[-1].append(float(i))
    data.close()

    return outList


def parseColvar(colvar):

    """Parse COLVAR files

    Args:
        colvar (string): COLVAR file

    Returns:

        outList (numpy array, float): list with time and colvar values
    """

    outList = []
    data = open(colvar, 'r')

    for line in data.readlines():
        if line[0] != '#':
            outList.append(
                    map(float,
                        filter(None, line.split(' '))))
    return np.array(outList)


def isExe(path):

    """Check if program exists and is executable.

    Args:
        path (string): path to the program.

    Returns:
        (bool): True or False if the program is executable.

    """
    return os.path.isfile(path) and os.access(path, os.X_OK)

def findExe(fileName):

    """Finds the full path of a program and check if it is executable.

    Args:
        fileName (string): name of the program or path pointing to it.

    Returns:
        (string): full path to the program if it exists..

    """

    """In case it's an alias."""
    if fileName.startswith('$'):
        p=subprocess.Popen('which '+fileName,shell=True,stdout=subprocess.PIPE)
        fileName=p.communicate()[0][:-1]

    fpath, fname = os.path.split(fileName)

    if fpath:
        if isExe(fileName): return fileName
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            fullPath = os.path.join(path, fileName)
            if isExe(fullPath):
                return fullPath

    return None






if __name__ == "__main__":

    pass

