# script to change sign in the off-diagonal parts of self-energy
# Vladislav Pokorny; 2015-2016; pokornyv@fzu.cz

import scipy as sp
from scipy.interpolate import InterpolatedUnivariateSpline
from sys import argv,exit

try:
	NBand    = int(argv[1])
	NAtom    = int(argv[2])
	infile   = str(argv[3])
except IndexError:
	print 'Usage: python changesgnoff.py [NBand (int)] [NAtom (int)] [filename (str)]'
	exit()

Data_F = sp.loadtxt(infile)
NMats = len(Data_F)/(NAtom*NBand**2)
print '{0: 3d} Matsubaras on input.'.format(NMats)

def WriteBlock(block_F,NMats,f):
	for i in range(NMats):
		f.write('{0: .8f}\t{1: .8f}\t{2: .8f}\n'\
		.format(float(block_F[i,0]),float(block_F[i,1]),float(block_F[i,2])))

f = open(infile+'_sgnoff','w')

k = 0
for atom in range(NAtom):
	for i in range(NBand):
		for j in range(NBand):
			block_F = Data_F[k*NMats:(k+1)*NMats]
			if any([i==0 and j==1,i==1 and j==0]):
				block_F[:,1] = -block_F[:,1]
				block_F[:,2] = -block_F[:,2]
			WriteBlock(block_F,NMats,f)
			k = k+1
f.close()

