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
import shutil  # for removing directory tree
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
from tools import *


class tsSetup:
    """ Standard TS-PPTIS setup class. """

    def __init__(self, top, gro, mdp, ndx='', gmx='$GMX'):
        """Initialise TS-PPTIS setup.
        Args:
                top (string): path to topology file .top
                gro (string): path to structure file .gro
                mdp (string): mdp file for the MDs
                ndx (string, optional): path to groups definition file .ndx
                gmx (string, optional): path to the local gromacs executable.

        """

        print sectionDelimiter("INITIALISATION")

        """Check local gromacs installation."""
        self.gmx = findExe(gmx)
        if self.gmx != None:
            print 'Gromacs installation:\t\tOK'
        else:
            sys.exit('Error : invalid gmx path ' + gmx + '\n' +
                     'Make sure to have a working version of gromacs 5.X installed!')

        for label, filePath in zip(['top','gro','mdp', 'ndx'],[top, gro, mdp, ndx]):
            if os.path.isfile(filePath):
                print "%s file:\t\t\tOK" % label
            else:
                if label != 'ndx':
                    sys.exit('%s file not found: %s' (label, filePath))
                else:
                    print '%s file:\t\t\tnot found' % label

        self.gro = gro
        self.top = top
        self.mdp = mdp
        self.ndx = ndx


    def initWindow(self, path, window, traj, colvar, overwrite=False):
        """Initialize a window

        Args:
            path (string): path of the window directory
            traj (string): path to initial trajectory .trr/.xtc
            colvar (string): path to the colvar file of the input trajectory
            window (list): interfaces of the window
            overwrite (bool): whether to overwrite existing folder

        """

        print sectionDelimiter("CREATING WINDOW")

        # Absolute path
        path = os.path.abspath(path)

        # Add trailing / to path if not provided
        if path[-1] != '/':
            path += '/'

        if overwrite:
            print "Initialising window [O]:\t", path
        else:
            print "Initialising window:\t\t", path

        # Check if folder exists and if overwriting is allowed
        if os.path.isdir(path):
            if overwrite:
                shutil.rmtree(path)
            else:
                sys.exit('Refusing to overwrite directory')


        """Check and load trajectory data."""
        if os.path.isfile(traj):
            self.traj = traj
            self.trajData = md.load(traj, top=self.gro)
            print 'Trajectory_file:\t\tOK'
        else:
            sys.exit('Trajectory file not found: '+ traj)

        """Check and load colvar file."""
        if os.path.isfile(colvar):
            self.colvar = colvar
            print 'COLVAR file:\t\t\tOK'
        else:
            sys.exit('COLVAR file not found: ' + colvar)
        # Check length of interface list
        if len(window) != 3:
            sys.exit('Wrong number of elements as window interfaces')

        # Create a handy dictionary with the path folders
        pathTree = {'path' : path,
                      'data' : path + 'data/',
                      'run'  : path + 'run/',
                      'temp' : path + 'temp/'}
        for folder in pathTree.values(): os.makedirs(folder)

        # In the data/ directory create symlinks to the initial
        # trajectory data
        os.symlink(os.path.abspath(self.traj),
                   pathTree['data']+'00000.' + self.traj.split('.')[-1])
        os.symlink(os.path.abspath(self.colvar),
                   pathTree['data']+'00000.cv')

        # Copy mdp file for generating the final TPR files.
        # The step number will be set to TMAX
        shutil.copy(os.path.abspath(self.mdp),
                   pathTree['temp']+'md.mdp')

        # 1. Copy the plumed config file from the script directory to the
        # window directory.
        # *** WILL PROBABLY BREAK IF TS-PPTIS IS LOADED AS A LIBRARY FROM
        # A DIFFERENT WORKING DIRECTORY ***
        # 2. Replace placeholders with window limits

        with open('plumed.dat', 'r') as handle:
            committorText = handle.read()

        committorText = committorText.replace('__LL__', str(
            window[0])).replace('__UL__', str(window[2]))

        with open(path + 'run/plumed_fw.dat', 'w') as handle:
            handle.write(committorText.replace('__COLVAR__', path+'run/COLVAR_FW'))

        with open(path + 'run/plumed_bw.dat', 'w') as handle:
            handle.write(committorText.replace('__COLVAR__', path+'run/COLVAR_BW'))


        # Determine XTC, TRR and COLVAR stride and timestep. WIll be written in window.cfg

        config = {}
        config['interfaces'] = ':'.join(map(str, window))

        with open(self.mdp,'r') as handle:
            for line in handle.readlines():
                if 'nstxout-compressed' in line:
                    config['xtc_stride'] = int(line.split()[2])
                elif 'nstxout' in line:
                    config['trr_stride'] = int(line.split()[2])
                elif 'dt' in line:
                    config['timestep']   = float(line.split()[2])

        with open(path+'run/plumed_bw.dat','r') as handle:
            for line in handle.readlines():
                if 'print' in line.lower():
                    for arg in line.split():
                        if 'stride' in arg.lower():
                            config['colvar_stride'] = int(arg.split('=')[1])

        # Complain if not determined
        if not 'xtc_stride' in config:
            sys.exit('XTC stride not found in ' + self.mdp)
        elif not 'timestep' in config:
            sys.exit('Timestep not found in ' + self.mdp)
        elif not 'colvar_stride' in config:
            sys.exit('COLVAR stride not found in ' + path+'run/plumed_bw.dat')

        print 'TRR stride:\t\t\t',config['xtc_stride']
        print 'COLVAR stride:\t\t\t', config['colvar_stride']

        # Initialize a config file. Can be useful for storing paths and
        # various configurations. See the config file in the old implementation
        with open(path + 'window.cfg', 'w') as handle:
            # write timestamp, interfaces, trrStride and colvarStride
            initText = '#' + timestamp() + '\n'
            for key in config.keys():
                initText += '{:20s} = {:20s}\n'.format(key, str(config[key]))
            handle.write(initText)

        # Update 17/01/18: structure differs.
        tpsAccEntry(path+'tps_acc.log', 0, len(self.trajData), 0, 0, 'A', 'B', 0, 0, 0, 0)

        # Hint at gromacs command for running the dynamics

        print '\n *** Example gromacs commands for BW and FW replicas: ***\n'

        for r in ['bw','fw']:
            print '%s mdrun -s %s.tpr -plumed %s.dat -deffnm %s ; ' % (
                    self.gmx,
                    pathTree['run'] + r,
                    pathTree['run'] + 'plumed_' + r,
                    pathTree['run'] + r)

    def setUpTPS(self, path):
        """Setup a TPS run

        Args:
            path (string): path of the window directory

        """

        # Absolute path
        path = os.path.abspath(path)

        # Add trailing / to path if not provided
        if path[-1] != '/':
            path += '/'

        # Create a handy dictionary with the path folders
        pathTree = {'path' : path,
                      'data' : path + 'data/',
                      'run'  : path + 'run/',
                      'temp' : path + 'temp/'}

        # Gromacs log file. Passed to runGmx
        gmxLog = path + 'gmx.log'

        # Determine whether the folder is a window by the presence of window.cfg
        if os.path.isfile(path + 'window.cfg'):
            print sectionDelimiter("SETUP RUN")
        else:
            sys.exit('Error: the folder does not seem to be a TS-PPTIS window')

        print 'Setting up run in:\t\t', path

        continuation = False
        # Open tps_acc.log, which holds info about accepted trajectories
        tpsAcc = parseTpsAcc(path+'tps_acc.log')
        # Number of accepted trajectories
        runNumber = len(tpsAcc)


        if runNumber > 1:
            continuation = True  # The first is the initial so, > 1 is continuation

        print 'First run:\t\t\t', not continuation

        config = parseConfig(path + 'window.cfg')
        print 'Interfaces:\t\t\t', config['interfaces']

        # Delete everything in the temp/ and run/ subdirectory. Keep plumed files
        try:
            for fileName in os.listdir(pathTree['temp']):
                if fileName != 'md.mdp':
                    os.remove(pathTree['temp']+fileName)

            for fileName in os.listdir(pathTree['run']):
                if not fileName.startswith('plumed'):
                    os.remove(pathTree['run']+fileName)
        except:
            pass

        # ------------ RECOVER PREVIOUS TRAJECTORY ----------

        prevRun = pathTree['data'] + '%05d' % (runNumber - 1)
        print 'Source trajectory data:\t\t', prevRun

        # load previous path:
        if os.path.isfile(prevRun + '.xtc'):
            prevTrajExt = '.xtc'
        elif os.path.isfile(prevRun + '.trr'):
            prevTrajExt = '.trr'
        else:
            sys.exit('Trajectory file not found:', prevRun + ' [xtc/trr]')


        # Get number of frames of the initial trajectory.
        prevTraj = md.load(prevRun + prevTrajExt, top=self.gro)
        pathLength = len(prevTraj)
        print "Path length:\t\t\t%d (%.1f ps)" % (pathLength,
                    pathLength*config['timestep']*config['trr_stride'])


        # Define TMAX
        tmax = (pathLength*config['timestep']*config['trr_stride'])/np.random.random()

        print 'TMAX:\t\t\t\t%.3f ps' % tmax

        # Write TMAX in mdp file:
        setTmax(pathTree['temp']+'md.mdp', tmax, config['timestep'])

        # Define shooting point and dump gro file
        point = shootingPoint(prevRun + '.cv', config['interfaces'])

        print 'Shooting point:\t\t\t', point[1]
        print 'Shooting frame:\t\t\t', point[0]
        print 'Shooting frame LPF:\t\t', point[2]
        extractFrame(point[1], prevTraj, self.gro,
                     prevRun + '.cv', pathTree['temp'] + 'frame.gro',
                     trajStride=config['xtc_stride'],
                     colvarStride=config['colvar_stride'])

        print '\nInitialising FW replica velocities...\t\t',

        # Generate tpr for velocity generation
        cmd = '%s grompp -c %s -f %s -p %s -maxwarn 1 -o %s -po %s' % (
            self.gmx, path + 'temp/frame.gro', './invert.mdp',
            self.top, path + 'temp/genvel.tpr', path + 'temp/mdout.mdp')

        runGmx(cmd, gmxLog, 'Generating tpr for velocity generation')

        cmd = '%s mdrun -s %s -deffnm %s' % (
            self.gmx, path + 'temp/genvel.tpr', path + 'temp/genvel')

        runGmx(cmd, gmxLog, 'Running 1 step')

        print 'Done'

        # Invert the velocities
        print 'Inverting velocities for the BW replica...\t',

        with open(pathTree['temp'] + 'genvel_inverted.gro', 'w') as handle:
            handle.write(
                formatGro(
                    invertGro(
                        parseGro(pathTree['temp'] + 'genvel.gro'
                                 ))))

        print 'Done'



        # Generating TPR files for FW and BW replicas

        print 'Generating TPR files for FW and BW replicas...\t',

        cmd = '%s grompp -c %s -f %s -p %s -maxwarn 1 -o %s -po %s' % (
            self.gmx, pathTree['temp'] + 'genvel.gro', pathTree['temp']+'md.mdp',
            self.top, pathTree['temp'] + 'fw.tpr', path + 'temp/mdout.mdp')

        # Use ndx if specified
        if self.ndx != '': cmd += '-n ' + self.ndx

        runGmx(cmd, gmxLog, 'Generating TPR file for FW replica')

        cmd = '%s grompp -c %s -f %s -p %s -maxwarn 1 -o %s -po %s' % (
            self.gmx, pathTree['temp'] + 'genvel_inverted.gro', pathTree['temp']+'md.mdp',
            self.top, pathTree['temp'] + 'bw.tpr', path + 'temp/mdout.mdp')

        # Use ndx if specified
        if self.ndx != '': cmd += '-n ' + self.ndx

        runGmx(cmd, gmxLog, 'Generating TPR file for BW replica')

        print 'Done'

        # Moving the TPR file in the run/ subdir
        os.rename(pathTree['temp'] + 'fw.tpr', pathTree['run'] + 'fw.tpr')
        os.rename(pathTree['temp'] + 'bw.tpr', pathTree['run'] + 'bw.tpr')


    def finalizeTPS(self, path):
        """Setup finalize a TPS run

        Args:
            path (string): path of the window directory

        """
# tps.info structure in TSPPTIS 1:
#      TIME      BW:    TIS CV     SIDE CROSS NR. STOP     FW:    TIS CV     SIDE CROSS NR. STOP
#     0.000      BW   6.062549591    1    0    0    0      FW   6.062549591    1    0    0    0      [INIT]
#     1.000      BW   6.338336945    1    0    0    0      FW   6.277128220    1    0    0    0      STOP    0

        # Absolute path
        path = os.path.abspath(path)

        # Add trailing /
        if path[-1] != '/':
            path += '/'

        # Create a handy dictionary with the path folders
        pathTree = {'path' : path,
                      'data' : path + 'data/',
                      'run'  : path + 'run/',
                      'temp' : path + 'temp/'}

        # Determine whether the folder is a window by the presence of window.cfg
        if os.path.isfile(path + 'window.cfg'):
            print sectionDelimiter("FINALIZING")
        else:
            sys.exit('Error: the folder does not seem to be a TS-PPTIS window')

        print 'Finalizing:\t\t', path

        config = parseConfig(path + 'window.cfg')

        window = map(float, config['interfaces'].split(':'))

        # Get run number by finding highest numbered trajectory in data/ dir (+1):
        runNumber = np.max([int(f.split('.')[0]) for f in os.listdir(pathTree['data']) if f.endswith('.cv')]) + 1

        print 'Run number:\t\t\t', runNumber

        print 'Interfaces:\t\t\t', config['interfaces']

        # Inverting BW replica and joining trajectories...
        # Follow the TSPPTIS 1 convention of getting frame 0 from the FW replica
        replTraj = [md.load(pathTree['run'] + 'bw.xtc', top=self.gro)[:0:-1],
                    md.load(pathTree['run'] + 'fw.xtc', top=self.gro)]  # *** CHANGE WITH TRAJFILE NAME ***


        endPoint = []
        jointColvar = []
        for i, repl in enumerate(('BW', 'FW')):

            # Load replica colvar
            replColvar = parseTxt(pathTree['run'] + 'COLVAR_' + repl)

            print "%s path length:\t\t\t%d (%.1f ps)" % (repl,len(replTraj[i]),
                        len(replTraj[i])*config['timestep']*config['trr_stride'])

            # Invert colvar if BW
            if repl == 'BW':
                replColvar = replColvar[:0:-1]
                replColvar[:, 0] = -replColvar[:, 0]
            jointColvar.append(replColvar)

        replTraj = md.join(replTraj)
        replTraj.save(pathTree['run'] + 'fulltraj.trr')
        pathLength = len(replTraj)
        jointColvar = np.concatenate([jointColvar[0], jointColvar[1]])

        # map CV values to -1/1 depending on which side of the central interface they lie
        jointSide = map(int, np.sign(jointColvar[:, 1] - window[1]))

        # iterate i,i+1 pairs of CV points to count transitions
        crossHist = np.zeros([len(jointSide) - 1])
        for i in range(len(jointSide) - 1):
            seq = (jointSide[i], jointSide[i + 1])
            if seq == (-1, 1):
                crossHist[i] = 1
            elif seq == (1, -1):
                crossHist[i] = -1

        crossCount = np.sum(crossHist == 1), np.sum(crossHist == -1)

        endPoint = ['A' if jointColvar[:,1][i] <= window[0] else 'B' if jointColvar[:,1][i] >= window[2] else ' ' for i in (0, -1)]

        # Accept if crossings = 0 and the trajectories reached the external interfaces (a.k.a. they were not killed by tmax
        # Or by any other external reason
        accepted = np.logical_and(
                        np.sum(crossCount) > 0,
                        ' ' not in endPoint)

        print 'Crossings (+/-):\t\t%d, %d' % (crossCount[0], crossCount[1])
        print 'Start/end side:\t\t\t%s -> %s' % (endPoint[0], endPoint[1]),

        if ' ' in endPoint:
            print '[TMAX REACHED]'
        else:
            print

        print 'Accepted\t\t\t%s' % accepted

        # Write tps.info in the run directory. Structure differs from TS-PPTIS 1
        tpsInfo = '# {:>8s} {:>8s} {:>10s} {:>8s} {:>8s}  {:s}\n'.format(
            'TIME', 'LPF', 'TIS CV', 'SIDE', 'CROSS', 'CROSS_SPEED')

        for i in range(len(jointColvar)):

            cross, crossSpeed = '', ''

            if i < len(crossHist):
                cross = str(int(crossHist[i]))
                # Calculate absolute crossing speed
                if crossHist[i] != 0:
                   crossSpeed =  str(
                           np.abs(
                               round(
                           (jointColvar[i+1, 1] - jointColvar[i, 1]) / (jointColvar[i+1, 0] - jointColvar[i, 0]),3)))


            tpsInfo += '{:10.3f} {:8d} {:10.3f} {:8d} {:>8s} {:>8s}\n'.format(jointColvar[i, 0],
                                                                       jointColvar[i,
                                                                                   0] >= 0,
                                                                       jointColvar[i, 1],
                                                                       jointSide[i],
                                                                       cross,
                                                                       crossSpeed)
        tpsInfo += '''
# Timestamp:\t\t%s
# Run number:\t\t%d
# Total crossings:\t%d
# Net crossing:\t\t%d
# Accept:\t\t%d
''' % (timestamp(), runNumber, np.sum(crossCount), crossCount[0] - crossCount[1], int(accepted))

        startFrame = np.where(jointColvar == 0)[0][0]

        # If accepted copy data to data/ directory:
        if accepted:
            # move traj
            shutil.move(pathTree['run'] + 'fulltraj.trr',
                    pathTree['data'] + '%05d.trr' % runNumber)

            # write trajectory info file
            with open(pathTree['data'] + '%05d.info' % runNumber, 'w') as handle:
                handle.write(tpsInfo)

            # write .cv file:
            with open(pathTree['data'] + '%05d.cv' % runNumber, 'w') as handle:
                for line in jointColvar:
                    handle.write(
                            ' '.join(
                                ['{:8.3f}'.format(field) for field in line]) + '\n')

# update tsp_acc.log
# ORIGINAL
# 0     0000       -          initial    1    501 1.0000   A  B  1   0.00       0     -      1 1 1 1
# NEW
# runNumber lenTraj startCV startSide startAB endAB +cross -cross |+cross|+|-cross| (+cross)-(-cross)

            tpsAccEntry(path+'tps_acc.log',runNumber,
                        len(jointColvar),
                        jointColvar[startFrame][1],
                        jointSide[startFrame],
                        endPoint[0], endPoint[1],
                        crossCount[0], crossCount[1],
                        np.sum(crossCount), crossCount[0] - crossCount[1])
        else:

            tpsAccEntry(path+'tps_rej.log',runNumber,
                        len(jointColvar),
                        jointColvar[startFrame][1],
                        jointSide[startFrame],
                        endPoint[0], endPoint[1],
                        crossCount[0], crossCount[1],
                        np.sum(crossCount), crossCount[0] - crossCount[1])

            # write trajectory info file
            with open(pathTree['data'] + 'rej_%05d.info' % runNumber, 'w') as handle:
                handle.write(tpsInfo)

class tsAnalysis:
    """ TS-PPTIS analysis class.

    It allows to calculate windows crossing probabilities and extract rates
    from a set of gromacs simulations.
    """

    def __init__(self, units='kJ/mol'):
        """Initialise the analysis.

        Args:

            ...

            units (string, optional): energy units of the input free energy

        """

        # FC: I guess here we will load anything that's needed to calculate
        # rates etc... (e.g. trajectories)

        if units not in ['kJ/mol', 'kcal/mol', 'kT']:
            print 'Warning:  unrecognised energy units, assuming kJ/mol'

        self.beta = 1 / 2.479
        self.crossInfo=[]
    
    def getCrossings(self, folderName='pptis00', rejectedName='pptis00/tps_rej.dat'):
        """Reads par files produced by ts-pptis and extracts information on crossing events to be used by getRates.

        Args:

            ...
            
            folderName (string, optional): path to the pptis output files #at the moment for one window only
            rejectedName (string, optional): path to the tps_rej.dat file
        
        """

        # FC: for the moment being, this function reads par files as in Giorgio's implementation
        # it will need to be adapted to read from our new output format.
        # another issue is that it works for one file at a time, I'm not sure how the crossings for each window
        # can be put together. 

        # WARNING: UNTESTED


        ###find the value of lambda0, for the moment old plumed format, alternatively we can implement it as input
        plumed=open(folderName+'/plumed.dat',"r")
        for line in plumed.readlines():
            if line[0]!='TIS':
                target=np.float(line[-1])
                break
        plumed.close()	   

        pp,pm = 0,0
        velSum, weightsSum = 0,0 #needed only for logging if we decide to keep it

        #load the par files, to be adapted 
        listPar=[]
        for fi in os.listdir(folderName+'/tps_data/'):
            if fi.endswith(".par"):
                listPar.append(fi)
        listSorted=sorted(listPar, key=natural_keys)

        #the output will be saved in a dictionary, but we can keep this as a logfile
        output=open(folderName+'/crossings.dat', 'w')                
        for fi in listSorted:
            crossData = analyzeCross(folderName+'/tps_data/'+fi, target) #see tools
            pp += crossData['p0p']
            pm += crossData['p0m']  
            weight = getWeightTraj(rejectedName, fi[:4]) #see tools
            velSum += crossData['vel']*weight
            weightsSum += weight
        
        #LOGGING and OUTPUTTING...                
        if crossData['vel'] > 0:    #but why only positive vel? 
            output.write('{:s}'.format(fi[:4]) +'\t'+ '{:.4f}'.format(crossData['vel'])+'\t' + '{:d}'.format(crossData['nrPos']) + '\t' +\
                    '{:d}'.format(crossData['nrNeg']) + '\t' +'{:d}'.format(weight)+ '\t' +crossData['end'] +'\n')  
            self.crossInfo.append([fi[:4],crossData['vel'],crossData['nrPos'],crossData['nrNeg'],weight,crossData['end']]) #a bit redundant at the moment
        output.close()

    
    def getRates(self, fes, Astate=-1, Bstate=-1, Acorr=0, Bcorr=0, indexTS=None, error=None, ratesFile='rates.dat', printFile=False):
        """Reads the free energy surface FES, TS-PPTIS crossing probabilities
        and ouputs, calculate the rate constants and print them to screen and/or to file.

        Args:

            fes (list (floats)): list containing the X and Y values of the calculated
                free energy
            ##We need to decide which format we want... for now list
            Astate (int, optional): index of the A state along the FES, if none provided
                assume minimum free energy point
            Bstate (int, optional): index of the B state along the FES, if none provided
                assume last point
            Acorr (float, optional): free energy correction to the A state
                (e.g. from external potentials)
            Bstate (int, optional): free energy correction to the B state
                (e.g. from external potentials)
            indexTS (int, optional): index of the TS in the FES vector provided,
                if none provided automatically look for the highest point in the FES
            error (float, optional): free energy error to calculate the rates range,
                if none provided automatically assume 1 kT
            ratesFile (string, optional): path to the file containing the probabilities of
             crossing windows
            printFile (bool, optional): activate/deactivate printing to file

        """

        # NOTE FC: This implementation HEAVILY relies on files as formatted by Giorgio's script
        # it should be adequately adapted if we decide to output that information in a different
        # format

        if Bstate == -1:
            As = np.argmin(fes[1])
        else:
            As = Astate

        if Bstate == -1:
            Bs = len(fes)-1
        else:
            Bs = Bstate

        if indexTS == None or indexTS > np.max(fes[0]) or indexTS < np.min(fes[0]):
            iTS = np.argmax([y for x, y in fes])
        else:
            iTS = indexTS
        # TS=np.argmax(val[:int(len(val)/4)])

        if error == None:
            error = 1 / self.beta

        offFES = [f - fes[1][TS] for f in fes[1]]

        norm = 0
        for i in range(As, iTS):
            norm += 0.5 * (np.exp(-beta * offFES[i + 1]) + np.exp(-beta * \
                offFES[i])) * (fes[0][i + 1] - fes[0][i])
        PA = np.exp(-beta * (offFES[iTS] + Acorr)) / norm

        PAlow = np.exp(-beta * (offFES[TS] + Acorr - error)) / norm
        PAupp = np.exp(-beta * (offFES[TS] + Acorr + error)) / norm

        norm = 0
        for i in range(iTS, Bs - 1):
            norm += 0.5 * (np.exp(-beta * offFES[i + 1]) + np.exp(-beta * \
                 offFES[i])) * (fes[0][i + 1] - fes[0][i])
        PB = np.exp(-beta * (offFES[iTS] + Bcorr)) / norm

        PBlow = np.exp(-beta * (offFES[TS] + Bcorr - error)) / norm
        PBupp = np.exp(-beta * (offFES[TS] + Bcorr + error)) / norm

        R = calcR(fes[0][iTS], ratesFile=ratesFile, crossInfo=self.crossInfo)

        kAB = PA * R * 1e12
        kABlow = PAlow * R * 1e12
        kABupp = PAupp * R * 1e12

        kBA = PB * R * 1e12
        kBAlow = PBlow * R * 1e12
        kBAupp = PBupp * R * 1e12

        # FIX ALL THIS PRINTING
        # Add also the times tau...
        if printFile == True:
            f = open('RatesOutput.dat', 'w')
            f.write("Rates in s^-1")
            f.write("%.3e" % kAB, "%.3e" % kBA)
            f.write("\nRates low")
            f.write("%.3e" % kABlow, "%.3e" % kBAlow)
            f.write("\nRates upp")
            f.write("%.3e" % kABupp, "%.3e" % kBAupp)
            f.close()

        print "\nRates in s^-1"
        print "%.3e" % kAB, "%.3e" % kBA
        print "\nRates low"
        print "%.3e" % kABlow, "%.3e" % kBAlow
        print "\nRates upp"
        print "%.3e" % kABupp, "%.3e" % kBAupp


def testAll():
    """ Runs a standard set of commands to test the correct functioning of TS-PPTIS. """

    # Test initialisation
    ts = tsSetup('../testfiles/topol.top',
                 '../testfiles/system.gro',
                 '../testfiles/md.mdp',
#                 gmx='/usr/bin/gmx')
                  gmx='/usr/local/gromacs/bin/gmx')

    ts.initWindow('../testfiles/pptis10',
                  [0.55,1,1.25],
                  '../testfiles/traj_fixed.xtc',
                  '../testfiles/COLVAR',
                  overwrite=True)

    ts.setUpTPS('../testfiles/pptis10')

    ts.finalizeTPS('../testfiles/pptis10')

if __name__ == "__main__":

    print("Running test...\n")
    testAll()
    print("\nDone!")
