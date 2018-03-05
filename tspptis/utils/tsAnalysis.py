"""
TS-PPTIS v2.0 Analysis Utility

F. Comitani: f.comitani@ucl.ac.uk
G. Mattedi: giulio.mattedi.16@ucl.ac.uk
@2017-2018

"""


############WARNING UNTESTED#################

#import __init__ as tsp
import tspptis as tsp
import argparse

if __name__ == "__main__":
    "Run a standard TS-PPTIS analysis sequence parsing inputs from command line."


    #Parse command line input
    parser = argparse.ArgumentParser(description="TS-PPTIS v2.0 Analysis Utility",
        formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=35))
    parser.add_argument("trans",type=float,help="position of the Transition State in the CV space")
    parser.add_argument("-fold",help="pptis windows output folder", default="./")
    parser.add_argument("-fes",help="free energy profile file (Plumed 2 format)" , default="./fes.dat")
    args = parser.parse_args()

    # Initialize tsAnalysis
    tsa = tsp.tsAnalysis(args.fold)
    tsa.getProbabilities()
    tsa.getCrossings(args.trans)

    fesLists=tsp.plumed2List(arg.fes)
    tsa.getRates(fesList) 

    #Check memory loss assumption
    tsa.endPointVel()

