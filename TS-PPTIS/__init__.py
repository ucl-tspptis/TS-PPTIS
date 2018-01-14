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

    def __init__(self, top, gro, traj, colvar, windows, par, mdp, ndx='', gmx='$GMX'):
        """Initialise TS-PPTIS setup.
        Args:
                top (string): path to topology file .top
                gro (string): path to structure file .gro
                traj (string): path to initial trajectory .trr/.xtc
                colvar (string): path to the colvar file of the input trajectory
                windows (string): path to text file containing information
                on the windows in the format left:center:right
                par (string): path to file containing frame time, CV value and
                   forward/backward direction information.
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
            # Might not be needed at this point. G.
            self.trajCV = parseTxt(colvar)
            print 'COLVAR file:\t\t\tOK'
        else:
            self.colvar = None
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

        """Check mdp file."""
        if os.path.isfile(mdp):
            self.mdp = mdp
            print "mdp file:\t\t\tOK"
        else:
            self.mdp = ''
            print 'mdp file:\t\t\tnot found'

        """Check for par file, if not foung generate from COLVAR"""
        if os.path.isfile(par):
            print "PAR file:\t\t\tOK"
        else:
            print "PAR file:\t\t\tnot found, it will be generated...",
            if self.colvar != None:
                generatePar(self.colvar, par)
            print "OK"
        self.par = par

    def initWindow(self, path, window, overwrite=False):
        """Initialize a window

        Args:
            path (string): path of the window directory
            window (list): interfaces of the window
            overwrite (bool): whether to overwrite existing folder

        """
        # Absolute path
        path = os.path.abspath(path)

        # Add trailing / to path if not provided
        if path[-1] != '/':
            path += '/'

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

        # Create the list of folders
        windowTree = [path, path + 'data/', path + 'run/', path + 'temp/']
        for folder in windowTree:
            os.makedirs(folder)

        # In the data/ directory create symlinks to the initial
        # trajectory data
        os.symlink(os.path.abspath(self.traj),
                   path + 'data/00000.' + self.traj.split('.')[-1])
        os.symlink(os.path.abspath(self.colvar),
                   path + 'data/00000.cv')
        os.symlink(os.path.abspath(self.par),
                   path + 'data/00000.par')

        # 1. Copy the plumed config file from the script directory to the
        # window directory.
        # *** WILL PROBABLY BREAK IF TS-PPTIS IS LOADED AS A LIBRARY FROM
        # A DIFFERENT WORKING DIRECTORY ***
        # 2. Replace placeholders with window limits

        shutil.copy('./plumed.dat', path)

        with open(path+'plumed.dat','r') as handle:
            committorText = handle.read()

        committorText = committorText.replace('__LL__',str(window[0])).replace('__UL__',str(window[2]))

        with open(path+'run/plumed_fw.dat','w') as handle:
            handle.write(committorText.replace('__COLVAR__','COLVAR_FW'))

        with open(path+'run/plumed_bw.dat','w') as handle:
            handle.write(committorText.replace('__COLVAR__','COLVAR_BW'))

        # Initialize a config file. Can be useful for storing paths and
        # various configurations. See the config file in the old implementation
        with open(path + 'window.cfg', 'w') as handle:
            initText = '# %s\ninterfaces = %s\n' % (
                timestamp(), ':'.join(map(str, window)))
            handle.write(initText)

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
        tpsAccHandle = open(path + 'tps_acc.log', 'a+')
        # Number of accepted trajectories
        tpsAccLines = sum([1 for line in tpsAccHandle])

        if tpsAccLines > 1:
            continuation = True  # The first is the initial so, > 1 is continuation

        print 'First run:\t\t\t', not continuation

        config = parseConfig(path + 'window.cfg')
        print 'Interfaces:\t\t\t', config['interfaces']

        # Delete everything in the temp/ subdirectory
        try:
            shutil.rmtree(path + 'temp/')
            os.makedirs(path + 'temp/')
        except:
            pass

        if not continuation:
            # Get number of frames of the initial trajectory.
            pathLength = len(self.trajData)
            print "Path length:\t\t\t", pathLength

            # Write first line of tps_acc.log. Using same structure as previous implementation
            # The if avoids situations where the run is set up twice without running and
            # the file ends up having two lines about the initial traj.
            if tpsAccLines == 0:
                tpsAccHandle.write(('0     0000       -          initial    1 '
                                    '{:>6} 1.0000   A  B  1   0.00       0     -      1 1 1 1'
                                    '\n').format(pathLength))
            # Define shooting point and dump gro file
            point = shootingPoint(self.par, config['interfaces'])
            print 'Shooting point:\t\t\t', point[1]
            print 'Shooting frame:\t\t\t', point[0]
            print 'Shooting frame LPF:\t\t', point[2]
            extractFrame(point[1], self.traj, self.gro,
                         self.colvar, path + 'temp/frame.gro')

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

            with open(path + 'temp/genvel_inverted.gro', 'w') as handle:
                handle.write(
                    formatGro(
                        invertGro(
                            parseGro(path + 'temp/genvel.gro'
                                     ))))

            print 'Done'

            # Generating TPR files for FW and BW replicas

            print 'Generating TPR files for FW and BW replicas...\t',

            cmd = '%s grompp -c %s -f %s -p %s -maxwarn 1 -o %s -po %s' % (
                self.gmx, path + 'temp/genvel.gro', self.mdp,
                self.top, path + 'temp/fw.tpr', path + 'temp/mdout.mdp')

            runGmx(cmd, gmxLog, 'Generating TPR file for FW replica')

            cmd = '%s grompp -c %s -f %s -p %s -maxwarn 1 -o %s -po %s' % (
                self.gmx, path + 'temp/genvel_inverted.gro', self.mdp,
                self.top, path + 'temp/bw.tpr', path + 'temp/mdout.mdp')

            runGmx(cmd, gmxLog, 'Generating TPR file for BW replica')

            print 'Done'

            # Moving the TPR file in the run/ subdir
            os.rename(path + 'temp/fw.tpr', path + 'run/fw.tpr')
            os.rename(path + 'temp/bw.tpr', path + 'run/bw.tpr')

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

        tpsAccHandle.close()



    def finalizeTPS(self,path):
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

        print 'Finalizing:\t\t',path

        # Determine whether the folder is a window by the presence of window.cfg
        if os.path.isfile(path + 'window.cfg'):
            print sectionDelimiter("FINALIZING")
        else:
            sys.exit('Error: the folder does not seem to be a TS-PPTIS window')

        config = parseConfig(path+'window.cfg')

        window = map(float,config['interfaces'].split(':'))

        print 'Interfaces:\t\t\t',config['interfaces']

        
        for repl in ('BW', 'FW'):

            # Load replica colvar
            replColvar = parseTxt(path+'run/COLVAR_'+repl)

            print '%s Path length:\t\t\t%.2f ns' % (repl,replColvar[-1,0] / 1000)

            # map CV values to -1/1 depending on which side of the central interface they lie
            replSide = map(int,np.sign(replColvar[:,1] - window[1]))
            # iterate i,i+1 pairs of CV points to count transitions
            crossCount = [0,0]
            for i in range(len(replSide)-1):
                seq = (replSide[i], replSide[i+1])
                if seq == (-1,1): crossCount[1] += 1
                if seq == (-1,1): crossCount[0] += 1

            print '%s crossings (+/-):\t\t%d, %d' % (repl,crossCount[0], crossCount[1])

            replEndPoint = 'A' if replColvar[-1][1] <= window[0] else 'B'

            print '%s end point:\t\t\t%s' % (repl,replEndPoint)


            # TO DO: reverse BW trajectory


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

    def getRates(fes, Astate=0, Bstate=-1, Acorr=0, Bcorr=0, indexTS=None, error=None, ratesFile='rates.dat', crossFile='crossings.dat', printFile=False):
        """Reads the free energy surface FES, TS-PPTIS crossing probabilities
        and ouputs, calculate the rate constants and print them to screen and/or to file.

        Args:

            fes (list (floats)): list containing the X and Y values of the calculated
                free energy
            ##We need to decide which format we want... for now list
            Astate (int, optional): index of the A state along the FES, if none provided
                assume first point
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
            crossFile (string, optional): path to the file containing information on the
             crossing events
            printFile (bool, optional): activate/deactivate printing to file

        """

        # NOTE FC: This implementation HEAVILY relies on files as formatted by Giorgio's script
        # it should be adequately adapted if we decide to output that information in a different
        # format

        As = Astate
        if Bstate == -1:
            Bs = len(fes)
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
            norm += 0.5 * \
                (np.exp(-beta * offFES[i + 1]) + exp(-beta *
                                                     offFES[i])) * (fes[1][i + 1] - fes[1][i])
        PA = np.exp(-beta * (offFES[iTS] + Acorr)) / norm

        PAlow = np.exp(-beta * (offFES[TS] + Acorr - error)) / norm
        PAupp = np.exp(-beta * (offFES[TS] + Acorr + error)) / norm

        norm = 0
        for i in range(iTS, Bs - 1):
            norm += 0.5 * \
                (np.exp(-beta * offFES[i + 1]) + exp(-beta *
                                                     offFES[i])) * (fes[1][i + 1] - fes[1][i])
        PB = np.exp(-beta * (offFES[iTS] + Bcorr)) / norm

        PBlow = np.exp(-beta * (offFES[TS] + Bcorr - error)) / norm
        PBupp = np.exp(-beta * (offFES[TS] + Bcorr + error)) / norm

        R = calcR(fes[0][iTS], ratesFile=ratesFile, crossFile=crossFile)

        kAB = PA * R * 1e12
        kABlow = PAlow * R * 1e12
        kABupp = PAupp * R * 1e12

        kBA = PB * R * 1e12
        kBAlow = PBlow * R * 1e12
        kBAupp = PBupp * R * 1e12

        # FIX ALL THIS PRINTING
        # Add also the times tau...
        if printFile == True:
            f = open('RatesOuptut.dat', 'w')
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
                 '../testfiles/traj_fixed.xtc',
                 '../testfiles/COLVAR',
                 '../testfiles/windows.dat',
                 '../testfiles/traj.par',
                 '../testfiles/md.mdp',
                 gmx='/usr/bin/gmx')

    # Test extractFrame
#    print extractFrame([0,1.2],
#            trajFile='../testfiles/traj_fixed_skipped.xtc',
#            topFile='../testfiles/system.gro',
#            colvarFile='../testfiles/COLVAR',
#            outFile='../testfiles/out.pdb',trajStride=10,colvarStride=1)

    # Test initialising a window
    # ts.initWindow('../testfiles/pptis10',[0.75,1,1.25], overwrite=True)

    # ts.setUpTPS('../testfiles/pptis10')

    ts.finalizeTPS('../testfiles/pptis10')


if __name__ == "__main__":

    print("Running test...\n")
    testAll()
    print("\nDone!")
