"""
TS-PPTIS v2.0 Setup Run Utility

F. Comitani: f.comitani@ucl.ac.uk
G. Mattedi: giulio.mattedi.16@ucl.ac.uk
@2017-2018

"""


############WARNING UNTESTED#################

#import __init__ as tsp
import tspptis as tsp
import argparse

if __name__ == "__main__":
    """Run a standard TS-PPTIS single run setup sequence parsing inputs from command line."""


    """Parse command line input."""
    parser = argparse.ArgumentParser(description="TS-PPTIS v2.0 Setup Run Utility",
        formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=35))
    parser.add_argument("-info",help="file containing all the setup info", default="./init.info")
    parser.add_argument("fold",help="pptis window output folder to setup")
    args = parser.parse_args()

    """Reinitialize tsSetup."""

    info=[line.strip() for line in open(args.info,'r').readlines()]
    ts = tsp.tsSetup(top=info[0],
                 gro=info[1],
                 mdp=info[2],
                 ndx=info[3],
                 gmx=info[4])

    """Setup the run."""

    ts.setUpRun(args.fold)

