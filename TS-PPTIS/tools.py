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
import subprocess as sub
import numpy as np
import mdtraj as md
import datetime


# Used as cdw in Popen
cwd = os.getcwd()


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
    data = open(nameFile, "r")

    for line in data.readlines():
        outList.append([])
        for i in line.split(':'):
            outList[-1].append(float(i))
    data.close()

    return outList


def parseTxt(colvar):
    """Parse COLVAR files

    Args:
        colvar (string): COLVAR file

    Returns:

        outList (numpy array, float): list with time and colvar values
    """

    outList = np.loadtxt(colvar)  # switched to np.loadtxt. More flexible

    # Remove duplicate timesteps UNTESTED!
    u, uIndeces = np.unique([l[0] for l in outList], return_index=True)

    return outList[uIndeces]


def parseConfig(config):
    """Parse config file

    Args:
        config (string): config file

    Returns:
        configDict (dict): dictionary with config entries
    """

    data = open(config, 'r')

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

    # In case it's an alias
    if fileName.startswith('$'):
        p = sub.Popen('which ' + fileName, shell=True, stdout=sub.PIPE)
        fileName = p.communicate()[0][:-1]

    fpath, fname = os.path.split(fileName)

    if fpath:
        if isExe(fileName):
            return fileName
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            fullPath = os.path.join(path, fileName)
            if isExe(fullPath):
                return fullPath

    return None


def findNearest(array, value):
    """Find nearest value in an array."""

    idx = (np.abs(array - value)).argmin()
    return array[idx]


def timestamp():
    """Print timestamp in the YYYY:MM:DD HH:MM:SS format."""

    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def extractFrame(cvValue, trajFile, topFile, colvarFile, outFile='out.pdb', trajStride=1, colvarStride=1, save=True):
    """Extract a frame with a given CV value from a trajectory. If files are not specified
    it will load self.trajData e self.trajCV. The function looks for the frame with the closest CV values

    Args:
        cvValue (float/list/tuple): CV value or CV range
        trajFile (string od mdtraj trajectory): trajectory file name or object
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

    if type(trajFile) == str:
        traj = md.load(trajFile, top=topFile)
    else:
        traj = trajFile

    colvar = parseTxt(colvarFile)

    # Subsample trajectory or COLVAR depending on who has highest stride
    # in order to have (almost) equal number of frames and colvar lines

    if trajStride == None or colvarStride == None:
        # Guess stride UNTESTED!
        ratio = traj.n_frames * 1.0 / len(colvar)
    else:
        ratio = trajStride * 1.0 / colvarStride

    if ratio < 1:
        traj = traj[::int(1.0 / ratio)]
    else:
        colvar = colvar[::int(ratio)]

    # IF cvValue is a float look for the value, else for the inclusive range
    if type(cvValue) == float:
        framesCV = findNearest(colvar[:, 1], cvValue)
        frames = np.where(colvar[:, 1] == framesCV)[0]
    else:
        cvEntries = np.where(np.logical_and(
            (colvar[:, 1] >= cvValue[0]), (colvar[:, 1] <= cvValue[1])))[0]
        framesCV = colvar[cvEntries, 1]
        frames = cvEntries

    if save:
        traj[frames].save(outFile)

    return frames, framesCV


def generatePar(colvarFile, parFile, direction=''):
    """Generate and save a PAR file starting from a COLVAR file or COLVAR data

    Args:
        colvarFile (string/iterable): COLVAR file name or nested list
        parFile (string): output PAR file name
        direction (list): list of directions for each frame.
            If not provided FW direction is assumed.

    """
    if type(colvarFile) == str:
        colvar = parseTxt(colvarFile)
    else:
        colvar = colvarFile

    # If directions not provided assume FW direction
    if len(direction) == 0:
        direction = [1] * len(colvar)
    elif len(direction) != len(colvar):
        sys.exit("COLVAR length and direction list do not have the same length")

    # Write PAR
    with open(parFile, 'w') as handle:
        for i in range(0, len(colvar)):
            handle.write(
                    '    '.join(
                        map(str, colvar[i]) + [str(int(direction[i]))]) + '\n')


def shootingPoint(parFile, interfaces):
    """ Define shooting point picking random CV value in the window range
    and extracting the closes frame in the CV space using the PAR file

    Args:
        parFile (string): PAR file of the trajectory
        interfaces (string): window in the X:Y:Z format

    Returns:
        frame (int): frame number
        CV values (float): cv value of the frame
        lpf (int): whether the frame is FW or BW
    """

    # Is redefining a function argument very bad coding? G.
    interfaces = np.array(map(float, interfaces.split(':')))[[0, 2]]

    # Read par file
    par = parseTxt(parFile)
    # Define shooting point CV value
    point = np.random.uniform(interfaces[0], interfaces[1])
    # Find closest frame
    frame = np.where(par[:, 1] == findNearest(par[:, 1], point))[0][0]

    # return frame number, CV value, lpf
    return frame, round(par[frame, 1], 3), int(par[frame, 2])


def sectionDelimiter(title, size=80, char='-'):
    """Prints the section title"""

    title = '[ ' + title + ' ]'
    titleLen = len(title)
    now = timestamp()
    if len(title) > size - len(now) - 1:
        sys.exit('Title too long')

    return char * int(np.floor((size - titleLen) / 2)) + title + char * int(np.ceil((size - titleLen) / 2) - len(now)) + now


def runGmx(cmd, logFile, logLine='', cwd=cwd):
    """ Run gromacs and save stdout and stderr to a logfile.
        Prints timestamp and custom message before the stdout/err.

    Args:
        cmd (string or iterable): command as string or argument iterable
        logFile (string): logfile
        logLine (string): debug message

    Returns:
        proc.returncode (int): process return code

    """
    if type(cmd) == str:
        cmd = cmd.split(' ')

    logHandle = open(logFile, 'a+')
    logHandle.write('\nDEBUG: %s\t%s\n' % (timestamp(), logLine))
    logHandle.flush()

    proc = sub.Popen(cmd, shell=False, stdout=logHandle,
                     stderr=logHandle, cwd=cwd)
    proc.wait()

    logHandle.close()

    if proc.returncode != 0:
        sys.exit('''\n\n!!!! GROMACS returned exit code: %d !!!!
Check the log file: \t%s
In section: \t\t%s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!''' % (proc.returncode, logFile, logLine))

    return proc.returncode


def tpsAccEntry(logFile, runNumber, trajLen, startingPoint, startingSide, bwAB, fwAB, leftCross, rightCross, crossCount, netCross):
    """ Add tps_acc.log/tps_rej.log entry for the combined BW+FW trjajectory.

    Args:
        logFile (string): file path
        runNumber (int): progressive run number. Initial is 0
        trajLen (float): trajectory length in ps
        startingPoint (float): starting point in CV space
        startingSide (int): whether the starting point is
                            before/after the central interface (-1/1)
        bwAB (string): which side the trajectory starts (A/B)
        fwAB (string): which side the trajectory end (A/B)
        leftCross (int): number of right to left crossings
        rightCross (int): number of left to right crossings
        crossCount (int): total number of crossings
        netCross (int): net number of crossings (rightCross - leftCross)

    """
    with open(logFile, 'a+') as handle:
        handle.write('{:10d} {:10d} {:10.3f} {:5d} {:>5s} {:>5s} {:5d} {:5d} {:5d} {:5d}\n'.format(
            runNumber,
            trajLen,
            startingPoint,
            startingSide,
            bwAB, fwAB,
            leftCross, rightCross,
            crossCount, netCross))

def parseTpsAcc(logFile):
    """ Parse tps_acc.log/tps_rej.log file.

    Args:
        logFile (string): file path


    Returns:
        parsed (nested list): tps_acc.log elements for each line

    """

    with open(logFile,'r') as handle:
        tpsAcc = handle.readlines()

    slices = (( 0,10),(11,21),(22,32),
              (33,38),(39,44),(45,50),
              (51,56),(57,62),(63,68),
              (69,74))
    parsed = []

    for entry in tpsAcc:
        line = []

        for s in slices:
            try:
                line.append(float(entry[s[0]:s[1]]))
            except:
                line.append(entry[s[0]:s[1]].replace(' ',''))

        parsed.append(line)
    return parsed

# --------------- GRO PARSING AND INVERTING -------------


def parseGro(groFile, slices=None):
    """ Returns a nested list with the gro fields.

        Args:
            gro (str): file name
            slices (nested list, optional): list of fields' limits (semiopen intervals)

        Returns:
            grofields (nested list): gro file fields

        The first two lines (title and atom number) and
        the last (box size) will be added verbatim.

        Assumes the standard number of decimals (3 for position,
        4 for velocity). See http://manual.gromacs.org/current/online/gro.html
        """

    with open(groFile, 'r') as handle:
        gro=handle.readlines()

    # If none given, use standard gro file structure
    if slices == None:
        slices=[[0, 5],
                  [5, 10],
                  [10, 15],
                  [15, 20],
                  [20, 28],
                  [28, 36],
                  [36, 44],
                  [44, 52],
                  [52, 60],
                  [60, 68]]

    # Copy verbatim the fist two line
    grofields=gro[:2]

    for line in gro[2:-1]:
        linefields=[]
        for s in slices:
            # For every line and every field try converting to float.
            # If it fails convert to string
            field=line[s[0]:s[1]]
            try:
                field=float(field)
            except:
                field=str(field).strip()
            linefields.append(field)
        grofields.append(linefields)

    # Append last line
    grofields.append(gro[-1])
    return grofields


def formatGro(grofields):
    """ Formats a parsed GRO file for writing

        Args:
            grofields (nested list): returned by parseGro

        Returns:
            text (string): formatted GRO

        Uses standard GRO precision.
    """
    text=''.join(map(str, grofields[:2]))
    for line in grofields[2:-1]:
        text += '%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n' % (
            line[0], line[1], line[2],
            line[3], line[4], line[5], line[6],
            line[7], line[8], line[9])

    text += grofields[-1]

    return text


def invertGro(grofields):
    """ Inverts the velocity of the atoms

        Args:
            grofields (nested list): returned by parseGro

        Returns:
            grofields (nested list): modified list

    """

    for i in range(2, len(grofields) - 1):
        grofields[i][7:]=[-1 * field for field in grofields[i][7:]]
    return grofields
# -------------------------------------------------------


if __name__ == "__main__":

    pass


###### In the following are the analysis tools, remove this line when done####


def calcR(posTS, ratesFile='rates.dat', crossFile='crossings.dat', debug=False):
    """ Calculates the R component of the rate constants with an iterative approach
    described by Juraszek et al. 2013

    Args:
        posTS (int): position of the Transition State along the path X-axis
        ratesFile (string, optional): path to the file containing the probabilities of
            crossing windows
        crossFile (string, optional): path to the file containing information on the
            crossing events
        debug (bool, optional): activate/deactivate the debug option

    Returns:
        R (float): final approximated value of R, if debug option is activated, returns
            vector with the value at each iteration
    """
    # NOTE FC: This implementation HEAVILY relies on files as formatted by Giorgio's script
    # it should be adequately adapted if we decide to output that information in a different
    # format

    data=readFile(ratesFile)
    lambdas, pmm, pmp, ppm, ppp=[], [], [], [], []
    numWindows=len(data)

    for line in data:
        if len(line) > 23:
            off=1
        else:
            off=0
        lambdas.append(float(line[2]))
        pmm.append(float(line[11 + off]))
        pmp.append(float(line[12 + off]))
        ppm.append(float(line[13 + off]))
        ppp.append(float(line[14 + off]))

    cross=readFile(crossFile)
    vel, we, pc, nc, ends=[], [], [], [], []

    for line in cross:
        vel.append(float(line[1]))
        we.append(int(line[4]))
        pc.append(int(line[2]))
        nc.append(int(line[3]))
        ends.append(line[5])

    ends=np.asarray(ends)
    Rtst=0.5 * np.sum(we) / np.sum([w / v for w, v in zip(we, vel)])

    pcw=np.sum([p * w for p, w in zip(pc, we)])
    ncw=np.sum([n * w for n, w in zip(nc, we)])

    pcw_ends=np.sum(np.ma.masked_where(ends == '-', we))
    ncw_ends=np.sum(np.ma.masked_where(ends == '+', we))

    p0p=pcw_ends * 1.0 / pcw
    p0n=ncw_ends * 1.0 / ncw

    diffL=99999
    lambda0=0
    for i in range(len(lambdas)):
        if np.abs(posTS - lambdas[i]) < diffL:
            diffL=np.abs(posTS - lambdas[i])
            lambda0=i

    A=[1, ppm[lambda0 + 1]]
    Ab=[1, pmp[lambda0 - 1]]
    R=[]

    for i in range(2, numWindows):
        m=i - 2
        if lambda0 + m + 1 > len(ppm) - 1 and lambda0 - m - 1 < 0:
            break

        AiNom=pmm[lambda0 + m] * ppm[lambda0 + m + 1] * A[m] * A[m + 1]
        AiDen=(pmp[lambda0 + m] * pmm[lambda0 + m + 1] + pmm[lambda0 + m] * ppm[lambda0 + m + 1])\
            * A[m] - pmp[lambda0 + m] * pmm[lambda0 + m + 1] * A[m + 1]
        if AiDen == 0:
            break
        if lambda0+m+1 <= len(ppm)-1:
            A.append(AiNom/AiDen)
        else:
            A.append(A[-1])

        AbiNom=ppp[lambda0 - m] * pmp[lambda0 - m - 1] * Ab[m] * Ab[m + 1]
        AbiDen=(ppm[lambda0 - m] * ppp[lambda0 - m - 1] + ppp[lambda0 - m] * pmp[lambda0 - m - 1])\
            * Ab[m] - ppm[lambda0 - m] * ppp[lambda0 - m - 1] * Ab[m + 1]
        if AbiDen == 0:
            break
        if lambda0-m-1 >= 0:
            Ab.append(AbiNom/AbiDen)
        else:
            Ab.append(Ab[-1])

    for m in range(0, len(A)):
        RiNom=0.5 * Rtst * \
            (p0n * ppm[lambda0] + p0p * pmp[lambda0]) * A[m] * Ab[m]
        RiDen=ppm[lambda0] * A[m] + pmp[lambda0] * Ab[m] + \
            (1 - ppm[lambda0] - pmp[lambda0]) * A[m] * Ab[m]
        R.append(RiNom / RiDen)

    if debug == True:
        return R
    return R[-1]
