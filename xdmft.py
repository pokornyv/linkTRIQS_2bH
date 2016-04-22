# dmft calculation using triqs ct-hyb solver for a 2-band Hubbard model 
# with off-diagonal hybridizations
# linking the triqs1.2 to dmft code
# Vladislav Pokorny; 2014-2016; pokornyv@fzu.cz

from subprocess import call,check_output,STDOUT,CalledProcessError
from os import listdir,getcwd,environ,mkdir,uname,getcwd
from string import split
from sys import exit,argv
from time import ctime,time

from libqmc import *
import params as p

# log file
if p.NBins == 1:
	logfname = 'xdmft.log'
	scflog   = 'scf.log'
else:
	logfname = 'xdmft_bin.log'
	scflog   = 'scf_bin.log'

MaxMats=OmegaN(p.NMats-1,p.beta) # Matsubaras begin with n=0
RetCode = 0			# return code for programs

stars  = '**********************************************************'
hashes = '##########################################################'

PrintAndWrite('\n'+hashes+'\n# New calculation started.\n'+hashes+'\n',logfname)
PrintAndWrite('Master node name: '+str(uname()[1])+', OS '+str(uname()[0])+'-'+str(uname()[4]),logfname)
if environ.has_key('PBS_NUM_NODES') and environ.has_key('PBS_NUM_PPN'):
	PrintAndWrite('Number of CPU cores: '+str(environ['PBS_NUM_NODES'])\
	+' x '+str(environ['PBS_NUM_PPN'])+'\n',logfname)
PrintAndWrite('Working directory: '+str(getcwd()),logfname)

# on MetaCentrum machines the number of cores for mpi is set automatically
if environ.has_key('PBS_NUM_NODES'): meta = True
else: meta = False

# writing beta and NMats to the dmft.in file for dmftx
[mu,n] = WriteMuAndN(logfname)
WriteDmftFile(p.beta,p.NMats)
PrintAndWrite('Chemical potential from previous iteration: '+str(sp.around(mu,3)),logfname)
mu_old = mu

# creating directories for triqs, dmft and bin files
for ff in ['dmft_files','triqs_files']:
	if ff not in listdir('.'): mkdir(ff,0744)
if p.NBins > 1 and 'bin_files' not in listdir('.'): mkdir('bin_files',0744)

PrintAndWrite('U='+str(p.U)+', J='+str(p.J)+', beta='+str(p.beta)+'\n',logfname)

if p.NStart == 1:
	# starting a new calculation, removing sigma from previous run
	if 'sigma' in listdir('.'):
		PrintAndWrite('Moving old sigma file to sigma.prev.',logfname)
		call(['mv','sigma','sigma.prev'])
	if p.useHF:
		#PrintAndWrite('Using HF self-energy estimate in first iteration, HF_SE='+str(HF_SE),logfname)
		#WriteSigmaHF(NAtom,NBand,NMats,p.beta,eq_band,HF_SE)
		pass

for NLoop in range(p.NStart,p.NStart+p.NIter):
	PrintAndWrite(stars+'\n'+'Starting DMFT loop '+str(NLoop)+'/'+str(p.NStart+p.NIter-1)\
	+' on '+ctime()+'\n'+stars+'\n',logfname)
	if p.offdiag: PrintAndWrite('Solver running in off-diagonal mode.',logfname)
	if p.NBins>1: PrintAndWrite('Starting binning mode, bins '+str(PrevBins+1)+'-'+str(PrevBins+NBins)+'\n',logfname)
	for ff in ['mom1.tmp','mom2.tmp','ene.tmp','shift.tmp','gtau','delta.w','delta.tau','n.dmft']:
		if ff in listdir('.'): call(['rm',ff])
	# prepare sigma file from previous iteration for dmft code
	if 'sigma'+str(NLoop-1) in listdir('.'):
		PrintAndWrite('sigma'+str(NLoop-1)+' file present\n',logfname)
		AppendSigma(NLoop-1,p.NAtom)
	# run dmft code #######################################
	t = time()
	if p.offdiag: dmftcode = 'dmftx3' 
	else: dmftcode = 'dmftx'
	PrintAndWrite(dmftcode+' START',logfname)
	if meta: RetCode = call(['mpirun',dmftcode])
	else:    RetCode = call(['mpirun','-np',str(p.np),dmftcode])
	t = time()-t
	if RetCode != 0:
		PrintAndWrite('Warning: dmftx exited with non-zero return code '+str(RetCode)+'\n',logfname)
		exit()
	PrintAndWrite(dmftcode+' STOP after '+str(sp.around(t,2))+' seconds.\n',logfname)
	# prepare input for ct-qmc solver
	call(['mv','gw','gw'+str(NLoop)])
	call(['cp','dmft.out','dmft_files/dmft.out.'+str(NLoop)])
	[mu,n] = WriteMuAndN(logfname)
	PrintAndWrite('Chemical potential from dmft.out : {0: .3f}'.format(mu),logfname)
	PrintAndWrite('Total density      from dmft.out : {0: .3f}'.format(n/p.NAtom)+' / atom',logfname)
	PrintAndWrite('{0: 3d}\t{1: .6f}\t{2: .6f}'.format(NLoop,mu,1.0*n/p.NAtom),'mu.dat')
	# run ct-hyb solver
	PrintAndWrite('Starting ct-hyb solver on '+str(ctime()),logfname)
	t = time()
	if meta: RetCode = call(['mpirun','pytriqs','../cthyb_U.py',str(NLoop)])
	else:    RetCode = call(['mpirun','-np',str(p.np),'pytriqs','../cthyb_U.py',str(NLoop)])
	t = time()-t
	PrintAndWrite('ct-hyb STOP after '+str(sp.around(t/60.0,2))+' minutes.\n',logfname)
	if p.NBins == 1:
		call(['mv','gw'+str(NLoop),'dmft_files/gw'])
	if 'STOP' in listdir('.'): 
		PrintAndWrite('STOP file present, breaking the DMFT loop.\n',logfname)
		call(['rm','STOP'])
		break

for ff in ['aw.dmft','mom1.tmp','mom2.tmp','ene.tmp','shift.tmp','gw','gtau','delta.w','delta.tau','n.dmft']:
	if ff in listdir('.'): call(['mv',ff,'dmft_files/'+ff])

PrintAndWrite(argv[0]+' done, '+ctime(),logfname)

