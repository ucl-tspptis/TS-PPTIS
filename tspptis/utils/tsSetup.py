"""
TS-PPTIS v2.0 Setup Utility

F. Comitani: f.comitani@ucl.ac.uk
G. Mattedi: giulio.mattedi.16@ucl.ac.uk
@2017-2018

"""


############WARNING UNTESTED#################

#import __init__ as tsp
import tspptis as tsp
import argparse
import numpy as np

if __name__ == "__main__":
    """Run a standard TS-PPTIS setup sequence parsing inputs from command line."""


    """Parse command line input."""

    parser = argparse.ArgumentParser(description="TS-PPTIS v2.0 Setup Utility",
        formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=35))
    parser.add_argument("-top",help="system topology (gromacs format)", default="./topol.top")
    parser.add_argument("-gro",help="system coordinates (gromacs format)", default="./system.gro")
    parser.add_argument("-mdp",help="simulation parameters file (gromacs format)", default="./md.mdp")
    parser.add_argument("-ndx",help="index file (gromacs format, optional)", default="")
    parser.add_argument("-gmx",help="gromacs executable", default="/usr/bin/gmx")
    parser.add_argument("-fold",help="pptis windows output folder", default="./")
    parser.add_argument("-win",help="windows list file", default="./winList.dat")
    parser.add_argument("-xtc",help="reference trajectory file (gromacs format)", default="./traj.xtc")
    parser.add_argument("-col",help="CVs colvar file (Plumed 2 format)", default="./COLVAR")
    parser.add_argument("-info",help="file where input data will be stored", default="./init.info")
    args = parser.parse_args()

    """Save input to a file."""

    # Not sure if this is a good idea, but it will save time
    # when setupping and finilizing the windows iteratively
    info=open(args.info,'w')
    info.write(args.top+'\n'+\
               args.gro+'\n'+\
               args.mdp+'\n'+\
               args.ndx+'\n'+\
               args.gmx+'\n')

    """Initialize tsSetup."""

    print args.gmx

    ts = tsp.tsSetup(top=args.top,
                 gro=args.gro,
                 mdp=args.mdp,
                 ndx=args.ndx,
                 gmx=args.gmx)

    """Automatically set up windows from a file."""
 
    i=0
    for line in open(args.win,'r'):
        i+=1
        line=line.split(',')
        line=[np.float(l) for l in line]
        ts.initWindow(args.fold+'/pptis'+'{:<02}'.format(i),
                      line,
                      args.xtc,
                      args.col,
                      overwrite=True)

        """ Setup the simulation files for the specific window."""

        ts.setUpRun(args.fold+'/pptis'+'{:<02}'.format(i))

