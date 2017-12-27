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
import mdtraj as md
import datetime

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

    data.close()

    # Remove duplicate timesteps UNTESTED!
    u, uIndeces = np.unique([l[0] for l in outList], return_index=True)

    return np.array(outList[uIndeces])

def parseConfig(config):
    """Parse config file

    Args:
        config (string): config file

    Returns:
        configDict (dict): dictionary with config entries
    """

    data = open(config,'r')

    configDict = {}
    for line in data.readlines():
        if line[0] != '#':
            entry = [word.strip() for word in line.split('=')]
            configDict[entry[0]] = entry[1]

    data.close()

    return configDict



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
        (string): full path to the program if it exists.

    """

    #In case it's an alias
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


def findNearest(array, value):
    """Find nearest value in an array."""
    
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def timestamp():
    """Print timestamp in the YYYY:MM:DD HH:MM:SS format."""

    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def extractFrame(cvValue, trajFile, topFile, colvarFile, outFile='out.pdb', trajStride=1, colvarStride=1, save=True):
    """Extract a frame with a given CV value from a trajectory. If files are not specified
    it will load self.trajData e self.trajCV. The function looks for the frame with the closest CV values

    Args:
        cvValue (float/list/tuple): CV value or CV range
        trajFile (string): trajectory file name
        topFile (string): topology for mdtraj.
        colvarFile (string): COLVAR file name
        outFile (string, optional): output file name
        trajStride (int, optional): trajectory stride
        colvarStride (int, optional): COLVAR stride
        save (bool, optional): Whether to save the frame or not

    Returns:
        frames (numpy array): frame number
        frameCV (numpy array): frame CV value

    """

    # In theory traj is already loadd in trajData, am I right? So passing self we already 
    # have all the information... let's fix this if this is the case.

    traj = md.load(trajFile,top=topFile)
    colvar = parseColvar(colvarFile)

    # Subsample trajectory or COLVAR depending on who has highest stride
    # in order to have (almost) equal number of frames and colvar lines

    if trajStride=None or colvarStride=None:
        # Guess stride UNTESTED!
        ratio= traj.n_frames*1.0/len(colvar)
    else:
        ratio= trajStride*1.0/colvarStride

    if ratio<1: 
        traj= traj[::int(1.0/ratio)]
    else:
        colvar= colvar[::int(ratio)]
 
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

def generatePar(colvarFile, parFile, direction=''):
    """Generate and save a PAR file starting from a COLVAR file

    Args:
        colvarFile (string): COLVAR file name
        parFile (string): output PAR file name
        direction (list): list of directions for each frame.
            If not provided FW direction is assumed.

    """

    colvar =  parseColvar(colvarFile)

    # If directions not provided assume FW direction
    if not direction:
        direction = [1]*len(colvar)
    elif len(direction) != len(colvar):
        sys.exit("COLVAR length and direction list do not have the same length")

    # Write PAR
    with open(parFile,'w') as handle:
        for i in range(0, len(colvar)):
            handle.write('    '.join(map(str,np.append(colvar[i],direction[i])))+ '\n')

def SectionDelimiter(title, size=80, char='-'):
    """Prints the section title"""

    title = '[ '+title+' ]'
    titleLen = len(title)
    now = timestamp()
    if len(title) > size-len(now)-1: sys.exit('Title too long')

    return char*int(np.floor((size-titleLen)/2)) + title + char*int(np.ceil((size-titleLen)/2)- len(now)) + now


if __name__ == "__main__":

    pass

