# script to change beta and/or number of Matsubaras in sigma file
# Vladislav Pokorny; 2015-2016; pokornyv@fzu.cz

import scipy as sp
from scipy.interpolate import InterpolatedUnivariateSpline
from sys import argv,exit

def OmegaN(n,beta):
	return (2.0*n+1.0)*sp.pi/beta

def WriteBlock(MatsFreq_F,ReSE_F,ImSE_F,f):
	for i in range(len(MatsFreq_F)):
		f.write('{0: .8f}\t{1: .8f}\t{2: .8f}\n'\
		.format(float(MatsFreq_F[i]),float(ReSE_F[i]),float(ImSE_F[i])))

NBand     = int(argv[1])
NAtom     = int(argv[2])
NMats_out = int(argv[3])
beta_out  = float(argv[4])
infile    = str(argv[5])

Data_F        = sp.loadtxt(infile)
NMats_in      = len(Data_F)/(NBand**2*NAtom)
MatsFreq_in_F = Data_F[0:NMats_in,0]
MaxMats       = MatsFreq_in_F[-1]

ReSE_out_F     = sp.zeros(NMats_out)
ImSE_out_F     = sp.zeros(NMats_out)

beta_in = sp.around(2.0*sp.pi/(MatsFreq_in_F[1]-MatsFreq_in_F[0]),5)

print '# Number of Matsubaras on input: {0: 3d}, output: {1: 3d}'.format(NMats_in,NMats_out)
print '# Energy cutoff:                 {0: .3f}'.format(MaxMats)
print '# Inverse temperature on input:  {0: .3f}, output: {1: .3f}'.format(beta_in,beta_out)

if NMats_out > NMats_in: print '# Extending frequency range, additional frequencies \
will be filled with the last data point.'
elif NMats_out < NMats_in: print '# Reducing number of frequencies.'
else: print '# Number of frequencies unchanged.'

f = open(infile+'b'+str(beta_out)+'NM'+str(NMats_out),'w')

MatsFreq_out_F = OmegaN(sp.array(range(NMats_out)),beta_out)
k = 0
for atom in range(NAtom):
	for band1 in range(NBand):
		for band2 in range(NBand):
			ReSE_F = Data_F[k*NMats_in:(k+1)*NMats_in,1]
			ImSE_F = Data_F[k*NMats_in:(k+1)*NMats_in,2]
			ReTail = ReSE_F[-1]
			ImTail = ImSE_F[-1]
			RealSE = InterpolatedUnivariateSpline(MatsFreq_in_F,ReSE_F)
			ImagSE = InterpolatedUnivariateSpline(MatsFreq_in_F,ImSE_F)
			for l in range(NMats_out):
				w = MatsFreq_out_F[l]
				if w <= MaxMats:
					ReSE_out_F[l] = RealSE(MatsFreq_out_F[l])
					ImSE_out_F[l] = ImagSE(MatsFreq_out_F[l])
				else:
					ReSE_out_F[l] = ReSE_F[-1]
					ImSE_out_F[l] = ImSE_F[-1]
			WriteBlock(MatsFreq_out_F,ReSE_out_F,ImSE_out_F,f)
			k = k+1
f.close()

