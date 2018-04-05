"""
TS-PPTIS v2.0 Analysis Utility

F. Comitani: f.comitani@ucl.ac.uk
G. Mattedi: giulio.mattedi.16@ucl.ac.uk
@2017-2018

"""


############WARNING UNTESTED#################

from __future__ import print_function
import tspptis as tsp
import argparse

if __name__ == "__main__":
    """Run a standard TS-PPTIS analysis sequence parsing inputs from command line."""


    """Parse command line input."""

    parser = argparse.ArgumentParser(description="TS-PPTIS v2.0 Analysis Utility",
        formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=35))
    parser.add_argument("folders",nargs='+',help="pptis windows folders")
    parser.add_argument("-trans",type=float,help="position of the Transition State in the CV space", default=None)
    parser.add_argument("-fes",type=str,help="free energy profile file (Plumed 2 format)" , default="./fes.dat")
    parser.add_argument("-error",type=float,help="error on the free energy. Default is kbT" , default=None)
    args = parser.parse_args()

    """Initialize tsAnalysis."""

    print('\n- Initalising...')

    tsa = tsp.tsAnalysis(args.folders)

    """Extract Probabilities."""

    print('- Calculating window probabilities...')

    tsa.getProbabilities()

    """Change Fes Format."""

    print('- Converting FES...')

    fesList=tsp.plumed2List(args.fes)

    """Find TS if none provided."""

    if args.trans==None:
        iTS = np.argmax(fesList[1])
        args.trans=fesList[0][iTS]
        print('- Finding Transition State...',args.trans)

    """Extract Crossings."""

    print('- Extracting crossings...')

    tsa.getCrossings(args.trans)

    """Extract Rates."""

    print('- Calculating rates...\n')

    tsa.getRates(fesList, valTS=args.trans, error=args.error)

    """Check Memory Loss Assumption."""

    #tsa.endPointVel()
    #tsa.checkMLA(plot=True)

