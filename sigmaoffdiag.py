# includes off-diagonal static self-energy into the sigma file
# Vladislav Pokorny; 2015-2016; pokornyv@fzu.cz

import scipy as sp
from scipy.interpolate import InterpolatedUnivariateSpline
from sys import argv,exit

try:
	NBand     = int(argv[1])
	NAtom     = int(argv[2])
	OffSigma1 = float(argv[3])
	OffSigma2 = float(argv[4])
	infile    = str(argv[5])
except IndexError:
	print 'Usage: python sigmaoffdiag.py [NBand (int)] [NAtom (int)] [SE1 (f)], [SE2 (f)] [infile (str)]'
	exit()

Data_F = sp.loadtxt(infile)

N     = len(Data_F)
NMats = N/(NBand**2*NAtom)

def WriteBlock(block_F,NMats,f):
	for i in range(NMats):
		f.write('{0: .8f}\t{1: .8f}\t{2: .8f}\n'\
		.format(float(block_F[i,0]),float(block_F[i,1]),float(block_F[i,2])))

f = open(infile+'_off','w')

k = 0
for atom in range(NAtom):
	for i in range(NBand):
		for j in range(NBand):
			block_F = Data_F[k*NMats:(k+1)*NMats]
			if any([i==0 and j==1,i==1 and j==0]):
				for w in range(NMats):
					block_F[w][1] = OffSigma1
					block_F[w][2] = 0.0
			if any([i==2 and j==3,i==3 and j==2]):
				for w in range(NMats):
					block_F[w][1] = OffSigma2
					block_F[w][2] = 0.0
			WriteBlock(block_F,NMats,f)
			k = k+1
f.close()

