"""
TS-PPTIS v2.0 Finalize Run Utility

F. Comitani: f.comitani@ucl.ac.uk
G. Mattedi: giulio.mattedi.16@ucl.ac.uk
@2017-2018

"""


############WARNING UNTESTED#################

#import __init__ as tsp
import tspptis as tsp
import argparse

if __name__ == "__main__":
    """Run a standard TS-PPTIS single run finalizing sequence parsing inputs from command line."""


    """Parse command line input."""

    parser = argparse.ArgumentParser(description="TS-PPTIS v2.0 Finalize Run Utility",
        formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=35))
    parser.add_argument("-info",help="file containing all the setup info", default="./init.info")
    parser.add_argument("fold",help="pptis window output folder to setup")
    args = parser.parse_args()

    """Re-initialize tsSetup."""
    info=open(args.info,'r') 
    ts = tsp.tsSetup(info[0],
                 info[1],
                 info[2],
                 info[3])

    """Finalize the run."""

    ts.finalizeRun(args.fold)

