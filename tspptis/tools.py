"""
TS-PPTIS tools

A set of utility functions for TS-PPTIS

F. Comitani: f.comitani@ucl.ac.uk
G. Mattedi: giulio.mattedi.16@ucl.ac.uk
@2017-2018

"""

from __future__ import division, print_function

import sys
import os
import subprocess as sub
import numpy as np
import mdtraj as md
import datetime
import re

# Used as cdw in Popen
cwd = os.getcwd()

def sigmoid(x,ref=0.0,beta=1.0):
    """A simple sigmoid function.

    Args:
        x (float):  input value to transform
        ref (float, optional): horizontal offset parameter
        beta (float, optional): function slope parameter

    Returns:
        (float): value of the sigmoid function at the input point
    """
   
    if beta*(x-ref)>np.log(9999): return 1
    return 1-1/(1+np.exp(beta*(x-ref)))


def convInt(text):
    """Takes an input and returns it as itself or as int if it's a digit.

    Args:
        text (string or number): input text 

    Returns:
        (string or int): the input converted to int if a digit
    """

    return int(text) if text.isdigit() else text

def natural_keys(text):
    """Help function to sort in human friendly order.

    Args:
        text (string or number): input text 

    Returns:
        (list,tring or int): list with the components of the splitted input
             converted to int if digits
    """
    return [ convInt(c) for c in re.split('(\d+)', text) ]

def getNumRows(nameFile):
    """Count the number of rows in a file. Returns 0 if the file does not exist.

    Args:
        nameFile (string): path to  text file 

    Returns:
        nrows (int) number of rows in the file
    """
    
    try:
        data = open(nameFile, 'r')
        nrows=len(data.readlines())
    except:
        nrows=0
    
    return nrows


def printLine(fieldName, field, end='\n'):
    print('{:<25} |{:s}{}'.format(fieldName, ' '*5, field), end=end)


def parseWindows(nameFile):
    """Parse a text file containing information on the TS-PPTIS windows.

    Args:
        nameFile (string): path to  text file containing information
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
    with open(colvar) as handle:
	raw = [map(float,line.strip().split())
		for line in handle.readlines() if line[0] != '#' and line != '\n']

    # Assume that the fist line is the maximum length a row can have
    outList = np.empty(shape=(len(raw), len(raw[0])))
    outList[:] = np.nan

    # Fill the nan matrix
    for i in range(outList.shape[0]):
	outList[i, :len(raw[i]):] = raw[i]


    # Remove duplicate timesteps UNTESTED!
    if len(outList.shape) > 1: # Avoid error if outList has 1 row
        u, uIndeces = np.unique(outList[:,0], return_index=True)
    else:
        uIndeces = 0

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
            try:
                configDict[entry[0]] = float(entry[1])
            except:
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

    # Handle case where multiple frames are selected
    if len(frames) > 1: frames = frames[:1]

    if save:
        traj[frames].save(outFile)

    return frames, framesCV


def shootingPoint(colvarFile, interfaces):
    """ Define shooting point picking random CV value in the window range
    and extracting the closest frame in the CV space using the colvar file

    Args:
        colvarFile (string): colvar file of the trajectory
        interfaces (string): window in the X:Y:Z format

    Returns:
        frame (int): frame number
        CV values (float): cv value of the frame
        lpf (int): whether the frame is FW or BW
    """

    # Is redefining a function argument very bad coding? G.
    interfaces = np.array(map(float, interfaces.split(':')))[[0, 2]]

    # Read colvar file
    colvar = parseTxt(colvarFile)
    # Define shooting point CV value
    point = np.random.uniform(interfaces[0], interfaces[1])
    # Find closest frame
    frame = np.where(colvar[:, 1] == findNearest(colvar[:, 1], point))[0][0]

    # return frame number, CV value, lpf
    return frame, round(colvar[frame, 1], 3), int(colvar[frame,0] >= 0)


def sectionDelimiter(title, size=80, char='_'):
    """Prints the section title"""

    title = '[ ' + title + ' ]'
    titleLen = len(title)
    now = timestamp()

    return '\n' + char * int(np.floor((size - titleLen) / 2)) +\
            title + char * int(np.ceil((size - titleLen) / 2) - len(now)) +\
            now + '\n'

def throwError(text, char='!', exitCode=1):
    """Thow error message and quit"""

    textLen = len(text)
    title = '[ ERROR ]'
    size = max(len(text)*2,len(title)*2)
    line = char * int(np.floor((size - textLen) / 2)) + title + char * int(np.ceil((size - textLen) / 2))

    print('\n\n\n%s\n\n%s\n\n%s' % (
            line,
            text,
            line))
    sys.exit(exitCode)


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
    header = not os.path.isfile(logFile) # whether to write header
    with open(logFile, 'a') as handle:
        if header:
            handle.write('# RUN_NUM LENGTH START_CV SIDE BW_END FW_END +CROSS -CROSS TOT_CROSS NET_CROSS\n')
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

    for entry in tpsAcc[1:]:   # ignore header
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


def setTmax(mdpFile, tmax, timestep):
    """ Sets nsteps to desider number in mdp file

        Args:
            mdpFile (string): mdp file
            tmax (float): max time in ps
            timestep (float): time step in ps

    """

    with open(mdpFile,'r') as handle:
        mdp = handle.readlines()

    with open(mdpFile,'w') as handle:
        for i in range(len(mdp)):
            if 'nsteps' in mdp[i][:6]:
                mdp[i] = 'nsteps = %d\n' % int(tmax/timestep)

            handle.write(mdp[i])


def getLambda(pathToFile):
    """Extracts the lambda value of a window from the window.cfg file in its 
    pptis output folder.
        
    Args:
        pathToFile (string): path to the folder containing th window.cfg file

    Returns:
        (float): the value of the central window
    """

    winFile=open(pathToFile+'/window.cfg',"r")
    for line in winFile.readlines():
        if line.startswith('interfaces'):
            return np.float(line.split(':')[1])


def getWeightTraj(pathToFile, index):
    """Extracts the weights from the rejected trajectory file, to be used
    when calculating velocities in crossing events.
        
    Args:
        pathToFile (string): path to the tps_rej.dat file
        index (string or int): trajectory index from which to count the weights

    Returns:
        (int): weight for the selected trajectory
    """

    trajFile=open(pathToFile, "r")
    counter=0

    for line in trajFile.readlines():
        if str(index) in line:
            counter+=1
    trajFile.close()

    return counter+1 #not sure why +1, double check!


def analyzeCross(fileName, target):
    """Reads an  ouput .info file and extracts information
    on the crossing events.

    Args:
        fileName (string): path to the input .info file        
        target (float): the position of the windows to analyze in the CV space
    
    Returns:
        info (dictionary): a ditctionary containing information on velocities
            and number of crossing events to be used in the calculation of the
            raates.
    """

    #Note FC: we need to decide if and what we want to log.

    #NOTE FC: (SUM of + - cross) is already contained in tps_acc.log, the only 
    #other info that's missing is the average velocity. For the moment let's leave it
    #like this, but let's consider moving everything into tps_acc.log

    cv,vel,crEvent,side=[],[],[],[]
    
    data = open(fileName,"r")
    for line in data.readlines():
        if line[0]!='#' and line!='\n':
            read=line.split()
            cv.append(np.float(read[1])) #not needed
            #clean this 
            if len(read)>4:
                crEvent.append(np.int(read[4]))
            side.append(np.int(read[3])) #not needed
            if len(read)>5:
                vel.append(np.float(read[-1]))
    data.close()
 
    vlist = []  #only for logging purposes
    cross = False    


######### Old version, please remove
#None of this seems to be needed if not for logging...
    n = 0
    average = 0
    p0p,p0m = 0,0
    posCross,negCross = 0,0
      
    for i in range(1,len(cv)):
    
        crossnow = False
        hit=''
        if  cv[i-1] < target  and  cv[i] >= target :
            cross = True
            direction = 1
            crossnow = True
        elif  cv[i-1] > target  and  cv[i] <= target :
            cross = True
            direction = -1
            crossnow = True
      
        if  cross and direction > 0 and cv[i] >= target+1: 
            hit= '+'
            p0p+=1
        elif  cross and direction < 0 and cv[i] <= target-1: 
            hit= '-'
            p0m+=1
                
        if  crossnow:
            average += np.abs(cv[i-1] - cv[i]) 
            n+=1 
            vlist.append(cv[i] - cv[i-1])  #only for logging
            if direction>0:
                posCross += 1
            elif direction<0:
                negCross += 1
           
        if cv[i] >= target+1: #why + and - 1???
            end = '+'
        elif cv[i] <= target-1:
            end = '-'
#########################################


    if len(vel) > 0:
        avrVel=np.mean(vel)
    else: avrVel = 0 

    if p0p>1: p0p=1
    if p0m>1: p0m=1
       
    info={}
    info['vel']=avrVel
    info['p0p']=p0p #only for logging
    info['p0m']=p0m #only for logging
    info['nrPos']=np.sum([1 if c>0 else 0 for c in crEvent])
    info['nrNeg']=np.sum([1 if c<0 else 0 for c in crEvent])
    if side[-1]>0:
        info['end']='+'   #we can cange this to a smarter 0/1 later
    else:
        info['end']='-'
        
    return info


def calcR(posTS, crossInfo, probInfo, debug=False):
    """ Calculates the R component of the rate constants with an iterative approach
    described by Juraszek et al. 2013

    Args:
        posTS (int): position of the Transition State along the path X-axis
        crossInfo (list, mixed): list containing information on the crossing 
            events on the TS
        probInfo (list, mixed): list containing information on the crossing 
            probabilities for each window events
        debug (bool, optional): activate/deactivate the debug option

    Returns:
        R (float): final approximated value of R, if debug option is activated, returns
            vector with the value at each iteration
    """

    lambdas, pmm, pmp, ppm, ppp=[], [], [], [], []
    numWindows=len(probInfo)

    for line in probInfo:
        lambdas.append(float(line[0]))
        pmm.append(float(line[1]))
        pmp.append(float(line[2]))
        ppm.append(float(line[3]))
        ppp.append(float(line[4]))

    vel, we, pc, nc, ends=[], [], [], [], []

    for line in crossInfo:
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

    lambda0=getClosestLambda(posTS,lambdas)
   
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
        if AiDen == 0 or np.isinf(AiNom):
            break
        if np.isnan(AiNom) or np.isnan(AiDen): #clean this mess...
           A.append(A[-1])
        if lambda0+m+1 <= len(ppm)-1:
            A.append(AiNom/AiDen)
        else:
            A.append(A[-1])

        AbiNom=ppp[lambda0 - m] * pmp[lambda0 - m - 1] * Ab[m] * Ab[m + 1]
        AbiDen=(ppm[lambda0 - m] * ppp[lambda0 - m - 1] + ppp[lambda0 - m] * pmp[lambda0 - m - 1])\
            * Ab[m] - ppm[lambda0 - m] * ppp[lambda0 - m - 1] * Ab[m + 1]
        if AbiDen == 0 or np.isinf(AbiNom):
            Ab.append(Ab[-1]) # to ensure they are the same length
            break
        if np.isnan(AbiNom) or np.isnan(AbiDen):
            Ab.append(Ab[-1])
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

def plumed2List(fileName):
    """Read a monodimensional free energy profile fes.dat file in Plumed 2 format
    and convert it to a format readable by getRates.

    Args:
         fileNam (string)e: path to Plumed 2 free energy file.

    Returns:
         fesList (nested list, float): fes information in the appropriate format.
    """

    fesList=[]
    f1,f2=[],[]
    for line in open(fileName,'r'):
        line=line.split()
        if line[0]!='#!' and len(line)>0:
            f1.append(float(line[0])),f2.append(float(line[1])) 
    
    fesList.append(f1), fesList.append(f2)

    return fesList

def getClosestLambda(posTS,lambdas):
    """Find the closest matching transition state lambda value from a list.

    Args:
         posTS (float): position of the transition state.
         lambdas (list, float): list of all available windows positions.

    Returns:
         lambda 0 (float): the closest matching TS lambda.
    """


    diffL=99999
    lambda0=0
    for i in range(len(lambdas)):
        if np.abs(posTS - lambdas[i]) < diffL:
            diffL=np.abs(posTS - lambdas[i])
            lambda0=i

    return lambda0

##################################################################################


if __name__ == "__main__":

    pass
