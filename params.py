import scipy as sp
from ConfigParser import SafeConfigParser

config = SafeConfigParser()
config.read('xdmft.in')

###########################################################
# read the inputs from xdmft.in ###########################

# mpi part ################################################

np = 6
if config.has_option('mpi','np'):
	np = int(config.get('mpi','np'))

# params part #############################################

jobname = str(config.get('jobname','jname'))
offdiag = bool(int(config.get('params','offdiag')))

useHF    = int(config.get('params','useHF'))
if useHF: HF_SE = sp.array(tuple(config.get('params','HF_SelfEnergy').split(',')))

NStart   = int(config.get('params','NStart'))
NIter    = int(config.get('params','NIter'))
NBins    = int(config.get('params','NBins'))
PrevBins = int(config.get('params','PrevBins'))
beta     = float(config.get('params','beta'))
U        = float(config.get('params','U'))
J        = float(config.get('params','J'))
NMats    = int(config.get('ct-hyb','NMats'))
NLeg     = int(config.get('ct-hyb','NLeg'))
NTau     = int(config.get('ct-hyb','NTau'))
NAtom    = int(config.get('params','NAtom'))
NBand    = int(config.get('params','NBand'))
bands_T  = tuple(config.get('params','band_names').split(','))
equiv_T  = tuple(config.get('params','equivalence').split(','))
SymmS    = bool(int(config.get('params','SymmSigma')))
SymmG    = bool(int(config.get('params','SymmGreen')))

NFit = 8
if config.has_option('params','NFit'):
	NFit  = int(config.get('params','NFit'))
mix = 1.0
if config.has_option('params','mix'):
	mix  = float(config.get('params','mix'))

# ct-hyb part #############################################

SpinFlip    = bool(int(config.get('ct-hyb','SpinFlip')))
PairHopping = bool(int(config.get('ct-hyb','PairHopping')))

partition = 'autopartition'
if config.has_option('ct-hyb','partition'):
	partition = str(config.get('ct-hyb','partition'))
measure_Gtau  = bool(int(config.get('ct-hyb','measure_Gtau')))
measure_Gleg  = bool(int(config.get('ct-hyb','measure_Gleg')))
measure_pert = False
if config.has_option('ct-hyb','measure_pert'):
	measure_pert  = bool(int(config.get('ct-hyb','measure_pert')))
measure_Gat = False
if config.has_option('ct-hyb','measure_Gat'):
	measure_Gat   = bool(int(config.get('ct-hyb','measure_Gat')))
measure_tc = False
if config.has_option('ct-hyb','measure_tc'):
	measure_tc    = bool(int(config.get('ct-hyb','measure_tc')))

NWarmup       = int(float(config.get('ct-hyb','NWarmup')))
NCycles       = int(float(config.get('ct-hyb','NCycles')))
LengthCycle   = int(float(config.get('ct-hyb','LengthCycle')))

###########################################################
# process the inputs ######################################
equiv_F = sp.empty(len(equiv_T),dtype = sp.int16)	# changing strings to numbers
for i in range(len(equiv_T)): equiv_F[i] = int(equiv_T[i])

if NBins > 1: measure_Gtau = True	# just in case so we don't lose statistics on G(tau)

if NTau < 2*NMats:
	print 'NTau too small, setting it to ',2*NMats+1
	NTau = 2*NMats+1

