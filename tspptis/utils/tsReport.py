from __future__ import division
import numpy as np
import glob
import os.path
import re
import argparse

def natural_sort(l):
    ''' Natural sorting '''
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)


if __name__ == "__main__":
    """Monitor PPTIS windows"""


    """Parse command line input."""

    parser = argparse.ArgumentParser(description="TS-PPTIS v2.0 Window Monitoring Utility",
        formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=35))
    parser.add_argument("fold",nargs='+',help="pptis windows")

    args = parser.parse_args()

#    folderList = natural_sort(
#        glob.glob(args.fold+'/pptis*/'))

    folderList = natural_sort(args.fold)



    header = '{:15s} {:20s} {:>6s} {:>6s} {:>6s} {:>6s} {:>6s} {:>6s} {:>6s} {:>10s}{:>10s}'.format('WINDOW','LAMBDA','TOT','ACC','REJ','AA','AB','BA','BB','ACC_LEN','REJ_LEN')
    print header
    print '_'*len(header) + '\n'

    totals = []
    printTotals = False # set to True if at least one pptis folder is found

    for folder in folderList:
        if not os.path.isfile(folder+'/window.cfg'): continue
        printTotals = True
        with open(folder+'/window.cfg','r') as handle:
            window = [line.split('=')[1].strip()
                        for line in handle if 'interfaces' in line][0]
        with open(folder+'/tps_acc.log','r') as handle:
            tps_acc = [filter(None,line.split(' '))
                        for line in handle if line[0] != '#'][1:]
            acc_length = sum([float(line[1]) for line in tps_acc]) / 1000
            acc = len(tps_acc)

        if os.path.isfile(folder+'/tps_rej.log'):
            with open(folder+'/tps_rej.log','r') as handle:
                tps_rej = [filter(None,line.split(' '))
                        for line in handle if line[0] != '#']
                rej_length = sum([float(line[1]) for line in tps_rej]) / 1000
                rej = len(tps_rej)
        else:
            rej_length = 0
            rej = 0

        count = (
            sum([1 for line in tps_acc if line[4] == 'A' and line[5] == 'A']),
            sum([1 for line in tps_acc if line[4] == 'A' and line[5] == 'B']),
            sum([1 for line in tps_acc if line[4] == 'B' and line[5] == 'A']),
            sum([1 for line in tps_acc if line[4] == 'B' and line[5] == 'B'])
            )
        name = filter(None, folder.split('/'))[-1]

        print '{:15s} {:20s} {:6d} {:6d} {:6d} {:6d} {:6d} {:6d} {:6d} {:10.3f}{:10.3f}'.format(
                    name,
                    window,
                    acc+rej,
                    acc,
                    rej,
                    count[0],count[1],count[2],count[3],
                    acc_length,
                    rej_length)
        totals.append([acc+rej,acc,rej,count[0],count[1],count[2],count[3],acc_length,rej_length])

    if printTotals:
        print '_'*len(header)
        totals = np.sum(totals,axis=0, dtype=float)
        print '{:15s} {:20s} {:6} {:6} {:6} {:6} {:6} {:6} {:6} {:10.3f}{:10.3f}'.format(
                'TOTAL','',
                int(totals[0]),
                int(totals[1]),
                int(totals[2]),
                int(totals[3]),
                int(totals[4]),
                int(totals[5]),
                int(totals[6]),
                totals[7],
                totals[8]
                )
