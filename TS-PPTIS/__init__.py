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


import sys, os, subprocess
import numpy as np

import matplotlib.pyplot as plt

import mdtraj as md

from tools import *

class tsSetup:
	""" Standard TS-PPTIS setup class. """

	def __init__(self, top, gro, traj, windows, ndx='', gmx='$GMX'):

		"""Initialise TS-PPTIS.
		Args:
			top (string): path to topology file .top
			gro (string): path to structure file .gro
			traj (string): path to initial trajectory .trr/.xtc
			windows (string): path to text file containing information 
				on the windows in the format left:center:right
			ndx (string, optional): path to groups definition file .ndx
			gmx (string, optional): path to the local gromacs executable.

	 	"""

		"""Check and load trajectory data."""
		try:
			self.trajData=md.load(traj, top=top)
			print 'Topology and trajectory files:\tOK'
		except:
			sys.exit('Error: invalid input topology/trajectory files ' + top + ' ' + traj)
		
		"""Check and load windows file."""	
		try:
			self.winList=parseWindows(windows)
			print 'PPTIS windows:\t\t\tOK'
		except:
			sys.exit('Error: invalid windows file ' + windows)
		
		"""Check ndx file."""
		if os.path.isfile(ndx): 
			self.ndx=ndx
		else:
			self.ndx=''
			print 'nxd file:\t\t\tnot found' 	


		"""Check local gromacs installation."""
		self.gmx=findExe(gmx)
		if self.gmx!=None:
			print 'Gromacs installation:\t\tOK'
	 	else:
			sys.exit('Error : invalid gmx path ' + gmx+'\n'+\
				'Make sure to have a working version of gromacs 5.X installed!')

			

def testAll():
	
	ts=tsSetup('../testfiles/topology.top',
			'../testfiles/structure.gro',
			'../testfiles/traj.xtc',
			'../testfiles/windows.dat',
			gmx='gmx')

if __name__ == "__main__":

	print("Running test...\n")
	testAll()
	print("\nDone!")
